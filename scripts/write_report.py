#!/usr/bin/env python3
"""Write one text report per genome from the finder working directory.

The report keeps the audit, the published coordinates and any step log that
actually said something. The tables and logs it reads are deleted with the
work directory.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

STEP_LOGS = (
    "predict_orfs.log",
    "hmm_search.log",
    "predict_trna.log",
    "trna_proximity.log",
    "extract_trna_region.log",
    "blast_mge.log",
    "extract_mge_region.log",
    "annotate_mge.log",
    "filter_confident_ie.log",
    "export_genome_features.log",
)
LOG_HEAD = 20
LOG_TAIL = 80


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _passed(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def _clip(text: str) -> str:
    lines = text.splitlines()
    if len(lines) <= LOG_HEAD + LOG_TAIL:
        return text.rstrip()
    omitted = len(lines) - LOG_HEAD - LOG_TAIL
    kept = lines[:LOG_HEAD] + [f"... {omitted} lines omitted ..."] + lines[-LOG_TAIL:]
    return "\n".join(kept)


def _candidate_line(row: dict[str, str]) -> str:
    status = "confident" if _passed(row.get("passed_confident", "")) else "rejected"
    reason = row.get("reject_reason", "").strip()
    reason_bit = f"  {reason}" if reason else ""
    return (
        f"{row.get('integrase_id', '')}\t{status}{reason_bit}\t"
        f"attL {row.get('attL_abs_lo', '')}-{row.get('attL_abs_hi', '')} "
        f"{row.get('attL_strand', '')}\t"
        f"tRNA {row.get('trna_start', '')}-{row.get('trna_end', '')} "
        f"{row.get('trna_strand', '')}\t"
        f"integrase {row.get('integrase_len_aa', '')} aa\t"
        f"exact {row.get('exact_anchored_bp', '')} bp\t"
        f"blast {row.get('best_attl_len_bp', '')} bp\t"
        f"{row.get('prodigal_hit_class', '')}"
    )


def _tsv_block(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "(none)"
    fields = list(rows[0].keys())
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(row.get(field, "") for field in fields))
    return "\n".join(lines)


def build_report(
    sample_dir: Path,
    *,
    sample: str,
    mge_finder: str = "",
    upstream_commit: str = "",
) -> str:
    """Assemble the report text for one sample directory."""
    integrases = _rows(sample_dir / "integrase_hits_summary.tsv")
    trna = _rows(sample_dir / "integrase_trna.tsv")
    blast = _rows(sample_dir / "mge_blast.tsv")
    audit = _rows(sample_dir / "ie_filter_audit.tsv")
    sites = _rows(sample_dir / "attachment_sites_genome.tsv")
    n_pass = sum(1 for row in audit if _passed(row.get("passed_confident", "")))
    n_fail = len(audit) - n_pass

    parts = [
        f"sample\t{sample}",
        f"code\t{mge_finder}",
        f"upstream_commit\t{upstream_commit}",
        "",
        "counts",
        f"integrase_hits\t{len(integrases)}",
        f"trna_pairs\t{len(trna)}",
        f"attL_blast_hits\t{len(blast)}",
        f"audit_rows\t{len(audit)}",
        f"confident\t{n_pass}",
        f"rejected\t{n_fail}",
        "",
        "candidates",
    ]
    if audit:
        parts.extend(_candidate_line(row) for row in audit)
    else:
        parts.append("(none)")
    parts.extend(["", "published_coordinates", _tsv_block(sites), "", "audit", _tsv_block(audit), "", "logs"])
    any_log = False
    for name in STEP_LOGS:
        path = sample_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            continue
        text = path.read_text(errors="replace").strip()
        if not text:
            continue
        any_log = True
        parts.extend(["", f"## {name}", _clip(text)])
    if not any_log:
        parts.append("(empty)")
    parts.append("")
    return "\n".join(parts)


def write_report(sample_dir: Path, out_path: Path, **kwargs: str) -> None:
    """Write the report for one sample."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(build_report(sample_dir, **kwargs))


def main() -> None:
    """Parse arguments and write one report."""
    parser = argparse.ArgumentParser(description="Write one IE finder report for a genome")
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mge-finder", default="")
    parser.add_argument("--upstream-commit", default="")
    args = parser.parse_args()
    write_report(
        args.sample_dir,
        args.out,
        sample=args.sample,
        mge_finder=args.mge_finder,
        upstream_commit=args.upstream_commit,
    )


if __name__ == "__main__":
    main()
