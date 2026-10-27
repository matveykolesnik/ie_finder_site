#!/usr/bin/env python3
"""BLAST each paired tRNA against the genomic window next to it.

The raw hits are the attL candidates; ``filter_confident_ie.py`` picks attL
from them.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from genome import contig_record, read_genome
from logger import get_logger
from v3ps_filters import DEFAULT_ATTL_WINDOW_BP, build_qseqid, closest_trna_rows, load_thresholds

logger = get_logger("annotate_mge_region")

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
"""Column order for BLAST tabular output (format 6)."""


def run_blast_on_region(
    query_rec: SeqRecord | list[SeqRecord],
    subject_seq: Seq,
    tmp_dir: str = ".",
) -> pd.DataFrame:
    """Run nucleotide BLAST of tRNA queries against a genomic window.

    Args:
        query_rec: tRNA query sequence record, or several in one search.
        subject_seq: Genomic window sequence.
        tmp_dir: Directory for the temporary query, subject and output files.

    Returns:
        BLAST hits as a DataFrame with columns listed in ``BLAST_COLS``; ``qseqid``
        is the id of the query record, whatever characters it holds.
        Returns an empty DataFrame when BLAST produces no alignments.
    """
    queries = query_rec if isinstance(query_rec, list) else [query_rec]
    # BLAST sees plain ids and the hits get the records' own back, so an id
    # with ':' or '|' in it reaches the caller unchanged.
    plain = {f"q{i}": q.id for i, q in enumerate(queries)}
    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp:
        qfa = os.path.join(tmp, "query.fa")
        sfa = os.path.join(tmp, "subject.fa")
        bout = os.path.join(tmp, "blast.tsv")

        SeqIO.write([SeqRecord(q.seq, id=key, description="") for key, q in zip(plain, queries)], qfa, "fasta")
        SeqIO.write([SeqRecord(subject_seq, id="subject", description="")], sfa, "fasta")

        subprocess.run(
            [
                "blastn", "-query", qfa, "-subject", sfa,
                "-outfmt", f"6 {' '.join(BLAST_COLS)}",
                "-word_size", "4", "-dust", "no",
                "-out", bout,
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        try:
            hits = pd.read_csv(bout, sep="\t", header=None, names=BLAST_COLS, dtype={"qseqid": str})
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=BLAST_COLS)
    hits["qseqid"] = hits["qseqid"].map(plain)
    return hits


def read_pairs(pairs_tsv: str) -> dict[str, tuple[str, str, int, int, str]]:
    """Query id -> (integrase id, contig, tRNA start, tRNA end, tRNA strand).

    The pairs are those ``extract_trna_region`` writes as queries, under the
    same ids. The id is only looked up here, never split: a contig name may
    hold the ':' that separates its fields.
    """
    pairs = pd.read_csv(pairs_tsv, sep="\t", dtype={"contig": str})
    out = {}
    for _, row in closest_trna_rows(pairs).iterrows():
        fields = (
            str(row["integrase_id"]), str(row["contig"]).strip(),
            int(row["trna_start"]), int(row["trna_end"]), str(row["trna_strand"]).strip(),
        )
        out[build_qseqid(*fields)] = fields
    return out


def main(
    genome_fasta: str,
    query_fa: str,
    pairs_tsv: str,
    out_tsv: str,
    tmp_dir: str = ".",
    window_size: int = DEFAULT_ATTL_WINDOW_BP,
) -> None:
    """BLAST every tRNA query in ``query_fa`` against its window and write all hits.

    The window is ``window_size`` bp on the 3' side of the tRNA. Each hit row
    carries the query id, ``integrase_id``, ``contig`` and ``wstart``, the
    1-based contig coordinate of the window's first base.

    Args:
        genome_fasta: Path to the sample genome FASTA file.
        query_fa: Path to ``mge_query.fa`` (tRNA sequences used as BLAST queries).
        pairs_tsv: Path to ``integrase_trna.tsv``, the pairs behind the queries.
        out_tsv: Output path for the raw hits (``mge_blast_raw.tsv``).
        tmp_dir: Directory for temporary BLAST files.
        window_size: Length of the searched window in bp.

    Raises:
        ValueError: when a query is not one of the pairs, rather than leaving
            its pair without hits.
    """
    genome = read_genome(genome_fasta)
    pairs = read_pairs(pairs_tsv)
    frames: list[pd.DataFrame] = []

    for rec in SeqIO.parse(query_fa, "fasta"):
        if rec.id not in pairs:
            raise ValueError(f"Query {rec.id} in {query_fa} is not a pair of {pairs_tsv}")
        integrase_id, contig, trna_start, trna_end, strand = pairs[rec.id]

        if strand == "+":
            window_start, window_end = trna_end, trna_end + window_size
        else:
            window_start, window_end = max(1, trna_start - window_size), trna_start
        contig_seq = contig_record(genome, contig).seq
        subject = contig_seq[window_start - 1 : min(len(contig_seq), window_end)]

        hits = run_blast_on_region(
            SeqRecord(rec.seq, id=rec.id, description=""),
            subject,
            tmp_dir,
        )
        if hits.empty:
            logger.info(f"No BLAST hits for {rec.id}")
            continue
        hits["integrase_id"] = integrase_id
        hits["contig"] = contig
        hits["wstart"] = window_start
        frames.append(hits)

    if frames:
        pd.concat(frames).to_csv(out_tsv, sep="\t", index=False)
        logger.info(f"{sum(len(df) for df in frames)} raw hits saved to {out_tsv}")
    else:
        open(out_tsv, "w").close()
        logger.info("No raw BLAST hits collected.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="BLAST tRNA queries against the genomic window next to each tRNA."
    )
    parser.add_argument("--ffn", required=True, help="Genome FASTA path.")
    parser.add_argument("--query", required=True, help="mge_query.fa path.")
    parser.add_argument("--pairs", required=True, help="integrase_trna.tsv path.")
    parser.add_argument("--out_tsv", required=True, help="Output mge_blast_raw.tsv path.")
    parser.add_argument("--tmp_dir", default=".", help="Temporary directory for BLAST.")
    parser.add_argument(
        "--params", "--config", dest="params", default="search_params.yaml",
        help="Search parameters (search_params.yaml).",
    )
    cli_args = parser.parse_args()
    thresholds = load_thresholds(cli_args.params)
    main(
        cli_args.ffn, cli_args.query, cli_args.pairs, cli_args.out_tsv, cli_args.tmp_dir,
        window_size=thresholds.attl_window_bp,
    )
