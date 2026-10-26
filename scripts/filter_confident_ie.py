#!/usr/bin/env python3
"""Apply confidence filters to integrative element candidates for one sample."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

from genome import contig_record, read_genome
from logger import get_logger
from tandems import find_tandems, read_integrases
from v3ps_filters import (
    FilterThresholds,
    build_qseqid,
    closest_trna_rows,
    evaluate_ie_candidate,
    load_thresholds,
    parse_orfs_gff,
)

logger = get_logger("filter_confident_ie")


def integrase_coords(hits_path: Path) -> dict[str, tuple[int, int, int]]:
    """Load integrase coordinates and protein length keyed by ORF identifier.

    Args:
        hits_path: Path to ``integrase_hits_summary.tsv``.

    Returns:
        Dictionary mapping ``orf_id`` to ``(start, end, protein_len_aa)``;
        coordinates are 1-based inclusive.
    """
    df = pd.read_csv(hits_path, sep="\t")
    out: dict[str, tuple[int, int]] = {}
    for _, row in df.iterrows():
        out[str(row["orf_id"])] = (int(row["start"]), int(row["end"]), int(row["protein_len_aa"]))
    return out


def load_raw_blast(raw_path: Path) -> pd.DataFrame:
    """Load raw BLAST output, returning an empty frame when the file is missing.

    Args:
        raw_path: Path to ``mge_blast_raw.tsv``.

    Returns:
        Parsed BLAST DataFrame, or an empty DataFrame when the file is absent.
    """
    if not raw_path.is_file() or raw_path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(raw_path, sep="\t", dtype={"contig": str})


def resolve_shared_trna(
    rows: list[dict],
    int_coords: dict[str, tuple[int, int, int]],
    contig_seqs: dict[str, str],
    thresholds: FilterThresholds,
) -> None:
    """Confident elements that reached the same tRNA: duplicates and tandems.

    Two candidates can reach one tRNA. An element with an integrase at each end
    is found from both (``integrase_at`` attR and attL), and elements integrated
    one after another into the same tRNA each find it, since every integration
    leaves the tRNA whole at the attR end. The same attL twice is one element:
    the first row keeps it, an integrase next to the tRNA before one next to
    attL, and the others are rejected as ``duplicate``. Different attLs are a
    tandem: the nearest element runs from the tRNA to its attL, and each further
    one from the attL copy before it -- its attR -- to its own attL. Where a
    further element's integrase lies before that copy, the two are nested
    rather than in tandem and are left as they are.

    Rows are updated in place: ``attR_lo``/``attR_hi``, ``attR_after`` (the
    integrase of the element before), the element span and the checks that
    depend on it.
    """
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("passed_confident"):
            key = (row["contig"], row["trna_start"], row["trna_end"], row["trna_strand"])
            groups[key].append(row)
    for (contig, trna_start, trna_end, strand), members in groups.items():
        if len(members) < 2:
            continue
        kept: list[dict] = []
        first: dict[tuple[int, int], dict] = {}
        for row in sorted(members, key=lambda r: r.get("integrase_at") != "attR"):
            attl = (int(row["attL_abs_lo"]), int(row["attL_abs_hi"]))
            if attl in first:
                row["passed_confident"] = False
                row["reject_reason"] = "duplicate"
                row["duplicate_of"] = first[attl]["integrase_id"]
                continue
            first[attl] = row
            kept.append(row)

        def from_trna(row: dict) -> int:
            if strand == "+":
                return int(row["attL_abs_lo"]) - int(trna_end)
            return int(trna_start) - int(row["attL_abs_hi"])

        kept.sort(key=from_trna)
        seq = contig_seqs[contig]
        for inner, outer in zip(kept, kept[1:]):
            int_lo, int_hi = sorted(int_coords[outer["integrase_id"]][:2])
            if strand == "+" and int_lo <= int(inner["attL_abs_hi"]):
                continue
            if strand == "-" and int_hi >= int(inner["attL_abs_lo"]):
                continue
            outer["attR_lo"], outer["attR_hi"] = int(inner["attL_abs_lo"]), int(inner["attL_abs_hi"])
            outer["attR_after"] = inner["integrase_id"]
            lo = min(int(outer["attL_abs_lo"]), outer["attR_lo"], int_lo)
            hi = max(int(outer["attL_abs_hi"]), outer["attR_hi"], int_hi)
            outer["ie_id"] = f"{outer['integrase_id']}:{contig}:{lo}-{hi}"
            outer["ie_len_nt"] = hi - lo + 1
            outer["ie_has_ambiguous_n"] = "N" in seq[lo - 1:hi]
            if thresholds.reject_ambiguous_n_ie and outer["ie_has_ambiguous_n"]:
                outer["passed_no_ambiguous_n"] = outer["passed_confident"] = False
                outer["reject_reason"] = "ambiguous_n"
            if thresholds.ie_min_nt > 0 and outer["ie_len_nt"] < thresholds.ie_min_nt:
                outer["passed_ie_len"] = outer["passed_confident"] = False
                outer["reject_reason"] = outer["reject_reason"] or "ie_too_short"


def filter_sample(
    *,
    sample: str,
    trna_path: Path,
    integrase_hits_path: Path,
    raw_blast_path: Path,
    orfs_gff_path: Path,
    fasta_path: Path,
    out_audit: Path,
    thresholds: FilterThresholds,
) -> int:
    """Evaluate every integrase-tRNA candidate of one sample and write the audit.

    Args:
        sample: Sample identifier (genome name).
        trna_path: Path to ``integrase_trna.tsv``.
        integrase_hits_path: Path to ``integrase_hits_summary.tsv``.
        raw_blast_path: Path to ``mge_blast_raw.tsv``.
        orfs_gff_path: Path to Prodigal ``orfs.gff``.
        fasta_path: Path to the assembly FASTA. attL and the element are read from it.
        out_audit: Output path for ``ie_filter_audit.tsv``.
        thresholds: Confidence filter thresholds.

    Returns:
        Number of integrase identifiers that passed all confidence filters.
    """
    trna_df = closest_trna_rows(pd.read_csv(trna_path, sep="\t", dtype={"contig": str}))
    out_audit.parent.mkdir(parents=True, exist_ok=True)
    if trna_df.empty:
        pd.DataFrame(columns=["integrase_id", "passed_confident"]).to_csv(
            out_audit, sep="\t", index=False
        )
        return 0

    int_coords = integrase_coords(integrase_hits_path)
    raw_blast = load_raw_blast(raw_blast_path)
    raw_by_q = raw_blast.groupby("qseqid") if not raw_blast.empty else None
    cds_by_contig = parse_orfs_gff(orfs_gff_path)
    genome = read_genome(fasta_path)

    audit_rows: list[dict] = []
    for _, trna_row in trna_df.iterrows():
        integrase_id = str(trna_row["integrase_id"])
        contig = str(trna_row["contig"]).strip()
        trna_start = int(trna_row["trna_start"])
        trna_end = int(trna_row["trna_end"])
        trna_strand = str(trna_row["trna_strand"]).strip()
        trna_len = trna_end - trna_start + 1

        if integrase_id not in int_coords:
            audit_rows.append({
                "sample": sample,
                "integrase_id": integrase_id,
                "passed_confident": False,
                "reject_reason": "integrase_not_in_hits",
            })
            continue

        int_start, int_end, int_aa = int_coords[integrase_id]
        qseqid = build_qseqid(integrase_id, contig, trna_start, trna_end, trna_strand)

        raw_hits = None
        if raw_by_q is not None:
            try:
                raw_hits = raw_by_q.get_group(qseqid)
            except KeyError:
                raw_hits = pd.DataFrame()

        audit = evaluate_ie_candidate(
            integrase_id=integrase_id,
            contig=contig,
            trna_start=trna_start,
            trna_end=trna_end,
            trna_strand=trna_strand,
            trna_len=trna_len,
            integrase_start=int_start,
            integrase_end=int_end,
            integrase_protein_aa=int_aa,
            raw_hits=raw_hits,
            contig_seq=str(contig_record(genome, contig).seq).upper(),
            cds_by_contig=cds_by_contig,
            thresholds=thresholds,
        )
        audit["sample"] = sample
        audit["qseqid"] = qseqid
        # Pairs made through attL (attl_side.py) say so; the others face the tRNA.
        at = trna_row.get("integrase_at", "attR")
        audit["integrase_at"] = at if isinstance(at, str) and at else "attR"
        audit["attR_lo"], audit["attR_hi"] = sorted((trna_start, trna_end))
        audit["attR_after"] = ""
        audit_rows.append(audit)

    contig_seqs = {c: str(r.seq).upper() for c, r in genome.items()}
    resolve_shared_trna(audit_rows, int_coords, contig_seqs, thresholds)
    if thresholds.tandem_search:
        # Elements in tandem at one tRNA: split at inner copies, add the ones past attL.
        audit_rows += find_tandems(
            audit_rows, read_integrases(integrase_hits_path, thresholds), int_coords,
            contig_seqs, cds_by_contig, thresholds, str(out_audit.parent),
        )
    n_passed = sum(bool(row.get("passed_confident")) for row in audit_rows)
    pd.DataFrame(audit_rows).to_csv(out_audit, sep="\t", index=False)
    logger.info(f"{sample}: {n_passed} confident IE(s) from {len(trna_df)} integrase-tRNA pair(s) -> {out_audit.name}"
                if thresholds.tandem_search else
                f"{sample}: {n_passed}/{len(trna_df)} confident IE(s) -> {out_audit.name}")
    return n_passed


def main() -> None:
    """Command-line entry point for per-sample confident IE filtering."""
    parser = argparse.ArgumentParser(
        description="Filter integrative element candidates to a confident set."
    )
    parser.add_argument("--sample", required=True, help="Sample identifier.")
    parser.add_argument("--trna", required=True, help="Path to integrase_trna.tsv.")
    parser.add_argument(
        "--integrase-hits",
        required=True,
        help="Path to integrase_hits_summary.tsv.",
    )
    parser.add_argument("--blast-raw", required=True, help="Path to mge_blast_raw.tsv.")
    parser.add_argument("--orfs-gff", required=True, help="Path to orfs.gff.")
    parser.add_argument("--fasta", required=True, help="Path to the assembly FASTA.")
    parser.add_argument(
        "--params", "--config", dest="params", default="search_params.yaml",
        help="Search parameters (search_params.yaml).",
    )
    parser.add_argument("--out-audit", required=True, help="Output ie_filter_audit.tsv path.")
    args = parser.parse_args()

    thresholds = load_thresholds(Path(args.params))
    filter_sample(
        sample=args.sample,
        trna_path=Path(args.trna),
        integrase_hits_path=Path(args.integrase_hits),
        raw_blast_path=Path(args.blast_raw),
        orfs_gff_path=Path(args.orfs_gff),
        fasta_path=Path(args.fasta),
        out_audit=Path(args.out_audit),
        thresholds=thresholds,
    )


if __name__ == "__main__":
    main()
