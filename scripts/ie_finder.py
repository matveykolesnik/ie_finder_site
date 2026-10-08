#!/usr/bin/env python3
"""FASTA in, attL/attR, integrase and element span out, in assembly coordinates.

    ie_finder.py strain.fasta outdir
    ie_finder.py genomes_dir outdir --cores 8

Each genome goes through the steps below in order, in a temporary work
directory. Genomes run in parallel, up to ``--cores`` at a time. Every step
writes its own log in the genome's work directory; the report keeps the logs
that are not empty. ``run.sh`` builds the conda environment and starts this
script.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

# The pipeline modules log FutureWarnings from pandas that say nothing about
# the run; keep them out of the step logs.
warnings.filterwarnings("ignore", category=FutureWarning)

import yaml

from annotate_mge_region import main as blast_attl_windows
from annotate_trna_proximity import find_nearby_trnas, parse_integrases, parse_trna, write_results
from export_genome_features import export_genome_features
from extract_trna_region import extract_trna_sequence
from filter_confident_ie import filter_sample
from hmm_search import parse_faa, parse_gff, parse_tblout, run_hmmsearch, write_outputs
from merge_search_params import merge_search_params
from predict_orfs import predict_with_prodigal
from v3ps_filters import FilterThresholds, thresholds_to_params
from write_report import write_report

ROOT = Path(__file__).resolve().parents[1]
FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fsa")
PFAM_PROFILES = (ROOT / "pfam" / "PF00589.hmm", ROOT / "pfam" / "PF22022.hmm")


class StepFailed(Exception):
    """A step of one genome raised; its log holds the traceback."""

    def __init__(self, sample: str, step: str, log: Path):
        # All three go to Exception so the error survives the trip back from
        # a worker process, which rebuilds it from these arguments.
        super().__init__(sample, step, log)
        self.sample, self.step, self.log = sample, step, log

    def __str__(self) -> str:
        return f"{self.sample}: step {self.step} failed, see {self.log}"


@dataclass(frozen=True)
class Run:
    """Settings shared by every genome of one run."""

    results: Path
    hmm: Path
    thresholds: FilterThresholds
    params_used: Path
    all_candidates: bool


def collect_inputs(path: Path) -> dict[str, Path]:
    """Sample name -> FASTA. A directory contributes every FASTA it holds.

    The sample name is the file name without its last suffix.

    Raises:
        ValueError: on a missing or empty input, a directory without FASTA
            files, or two files that share a sample name.
    """
    if path.is_dir():
        files = sorted(p for p in path.iterdir() if p.suffix in FASTA_SUFFIXES and p.is_file())
        if not files:
            raise ValueError(f"No FASTA files in {path}")
    elif path.is_file():
        if path.stat().st_size == 0:
            raise ValueError(f"FASTA is empty: {path}")
        files = [path]
    else:
        raise ValueError(f"Input not found: {path}")
    samples: dict[str, Path] = {}
    for fasta in files:
        sample = fasta.name.rsplit(".", 1)[0] if "." in fasta.name else fasta.name
        if not sample or sample.startswith("."):
            raise ValueError(f"Could not derive a sample name from {fasta}")
        if sample in samples:
            raise ValueError(f"Two inputs share the sample name {sample}")
        samples[sample] = fasta
    return samples


@contextmanager
def output_to(log: Path):
    """Send this process's stdout and stderr, and its children's, to ``log``.

    The redirect is on the file descriptors, so Prodigal and Aragorn write
    into the step log just as the Python logging does.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    # Line-buffered, so a print and the output of a program started right
    # after it reach the log in the order they happened.
    sys.stdout.reconfigure(line_buffering=True)
    saved = os.dup(1), os.dup(2)
    with log.open("w") as handle:
        os.dup2(handle.fileno(), 1)
        os.dup2(handle.fileno(), 2)
        try:
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved[0], 1)
            os.dup2(saved[1], 2)
            os.close(saved[0])
            os.close(saved[1])


def step(sample: str, out: Path, name: str, func, /, *args, **kwargs) -> None:
    """Run one step with its output in ``<out>/<name>.log``."""
    log = out / f"{name}.log"
    with output_to(log):
        try:
            func(*args, **kwargs)
        except BaseException:
            traceback.print_exc()
            raise StepFailed(sample, name, log) from None


def find_integrases(faa: Path, gff: Path, hmm: Path, out: Path, cutoff: str) -> None:
    """hmmsearch the ORFs and write the integrase tables."""
    tblout = out / "integrase_hits.txt"
    run_hmmsearch(str(faa), str(hmm), str(tblout), cutoff)
    write_outputs(
        parse_tblout(tblout), parse_faa(faa), parse_gff(gff),
        out / "integrase_hits_summary.tsv", out / "integrase_orfs.tsv",
    )


def predict_trnas(fasta: Path, out: Path) -> None:
    """Aragorn on the assembly, tRNA genes only, batch output."""
    subprocess.run(["aragorn", "-w", "-t", "-o", str(out), str(fasta)], check=True)


def pair_trnas(integrases: Path, trnas: Path, out: Path, max_distance: int) -> None:
    """Integrases with an opposite-strand tRNA within ``max_distance``."""
    pairs = find_nearby_trnas(parse_integrases(integrases), parse_trna(trnas), max_distance)
    write_results(pairs, out)


def run_sample(sample: str, fasta: Path, run: Run) -> None:
    """Every step for one genome, writing into ``<results>/<sample>/``."""
    out = run.results / sample
    out.mkdir(parents=True, exist_ok=True)
    t = run.thresholds
    gff, faa = out / "orfs.gff", out / "orfs.faa"
    step(sample, out, "predict_orfs", predict_with_prodigal,
         str(fasta), str(gff), str(out / "orfs.ffn"), str(faa))
    step(sample, out, "hmm_search", find_integrases, faa, gff, run.hmm, out, t.integrase_hmm_cutoff)
    step(sample, out, "predict_trna", predict_trnas, fasta, out / "trna.tsv")
    step(sample, out, "trna_proximity", pair_trnas,
         out / "integrase_hits_summary.tsv", out / "trna.tsv", out / "integrase_trna.tsv",
         t.trna_max_distance_bp)
    step(sample, out, "extract_trna_region", extract_trna_sequence,
         str(fasta), str(out / "integrase_trna.tsv"), str(out / "mge_query.fa"))
    step(sample, out, "blast_mge", blast_attl_windows,
         str(fasta), str(out / "mge_query.fa"), str(out / "mge_blast_raw.tsv"), str(out),
         window_size=t.attl_window_bp)
    step(sample, out, "filter_confident_ie", filter_sample,
         sample=sample,
         trna_path=out / "integrase_trna.tsv",
         integrase_hits_path=out / "integrase_hits_summary.tsv",
         raw_blast_path=out / "mge_blast_raw.tsv",
         orfs_gff_path=gff,
         fasta_path=fasta,
         out_audit=out / "ie_filter_audit.tsv",
         thresholds=t)
    step(sample, out, "export_genome_features", export_genome_features,
         fasta=fasta,
         audit=out / "ie_filter_audit.tsv",
         integrases=out / "integrase_hits_summary.tsv",
         trna=out / "integrase_trna.tsv",
         out_gff3=out / f"{sample}.ie.gff3",
         out_gbk=out / f"{sample}.ie.gbk",
         out_tsv=out / "attachment_sites_genome.tsv",
         sample=sample,
         all_candidates=run.all_candidates)


def upstream_commit() -> str:
    """The MGE_finder commit this copy was taken from, as UPSTREAM records it."""
    for line in (ROOT / "UPSTREAM").read_text().splitlines():
        if line.startswith("commit:"):
            return line.split(":", 1)[1].strip()
    return ""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", type=Path, help="Assembly FASTA, or a directory of them")
    parser.add_argument("outdir", type=Path, nargs="?", default=ROOT / "outdir",
                        help="Where the GFF3, GenBank and report go (default: outdir/ in this repository)")
    parser.add_argument("--cores", "-j", type=int, default=os.cpu_count() or 1,
                        help="Genomes processed at the same time (default: all CPUs)")
    parser.add_argument("--search-params", type=Path, default=None,
                        help="File whose keys override search_params.yaml")
    parser.add_argument("--all-candidates", action="store_true",
                        help="Also write candidates that have attL but failed a later filter")
    parser.add_argument("--keep-work", action="store_true",
                        help="Keep the temporary work directory")
    args = parser.parse_args(argv)

    try:
        samples = collect_inputs(args.input.resolve())
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 1
    # Checked before any step starts, so a typo in a key stops the run at once.
    try:
        thresholds = merge_search_params(
            ROOT / "search_params.yaml",
            args.search_params.resolve() if args.search_params else None,
        )
    except (OSError, ValueError, yaml.YAMLError) as exc:
        print(f"Search parameters: {exc}", file=sys.stderr)
        return 1
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    work = Path(tempfile.mkdtemp(prefix="ie_finder_site."))
    try:
        results = work / "results"
        params_used = work / "search_params.yaml"
        params_used.write_text(yaml.safe_dump(thresholds_to_params(thresholds), sort_keys=False))
        hmm = results / "combined" / "pfam_combined.hmm"
        hmm.parent.mkdir(parents=True)
        hmm.write_text("".join(p.read_text() for p in PFAM_PROFILES))
        run = Run(results, hmm, thresholds, params_used, args.all_candidates)

        failures: list[StepFailed] = []
        workers = max(1, min(args.cores, len(samples)))
        if workers == 1:
            for sample, fasta in samples.items():
                try:
                    run_sample(sample, fasta, run)
                except StepFailed as exc:
                    failures.append(exc)
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                jobs = {sample: pool.submit(run_sample, sample, fasta, run) for sample, fasta in samples.items()}
                for job in jobs.values():
                    try:
                        job.result()
                    except StepFailed as exc:
                        failures.append(exc)

        failed = {exc.sample for exc in failures}
        for sample in samples:
            sample_dir = results / sample
            if sample not in failed:
                for suffix in (".ie.gff3", ".ie.gbk"):
                    shutil.copyfile(sample_dir / f"{sample}{suffix}", outdir / f"{sample}{suffix}")
                print(f"GFF3:    {outdir / f'{sample}.ie.gff3'}")
                print(f"GenBank: {outdir / f'{sample}.ie.gbk'}")
            # A report is written whenever there are tables to report on, so a
            # failed genome still says how far it got.
            if (sample_dir / "integrase_hits_summary.tsv").is_file() or (sample_dir / "ie_filter_audit.tsv").is_file():
                write_report(
                    sample_dir,
                    outdir / f"{sample}.ie.report.txt",
                    sample=sample,
                    mge_finder=str(ROOT),
                    upstream_commit=upstream_commit(),
                    params_path=params_used,
                )
                if sample not in failed:
                    print(f"Report:  {outdir / f'{sample}.ie.report.txt'}")

        sys.stdout.flush()
        for exc in failures:
            print(f"\n{exc}\n{exc.log.read_text(errors='replace')}", file=sys.stderr)
        if failures or args.keep_work:
            print(f"Work left at {work}", file=sys.stderr)
        else:
            shutil.rmtree(work)
        return 1 if failures else 0
    except BaseException:
        # An error in the driver itself: say where the partial results are.
        print(f"Work left at {work}", file=sys.stderr)
        raise


if __name__ == "__main__":
    sys.exit(main())
