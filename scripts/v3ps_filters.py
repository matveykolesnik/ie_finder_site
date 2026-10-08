"""Confident integrative-element cascade, in the order of the Results section.

1. Integrase paired with an oppositely oriented tRNA, 3' ends within
   ``trna_max_distance_bp`` (``annotate_trna_proximity.py``).
2. Exact duplication of the tRNA 3' end of at least ``attl_candidate_min_bp``
   (perfect run from the 3' anchor).
3. Integrase longer than ``integrase_min_aa`` amino acids.
4. attL in a non-coding region (no overlap with a Prodigal CDS on the contig).
5. attL length: exact run of at least ``attl_exact_min_bp``, or the less
   specific search — BLAST length, mismatches allowed — of at least
   ``attl_relaxed_min_bp``.
6. No alignment gaps: the alignment length equals both the query span and the
   subject span.
7. Deduplication (``dedup_ie_representatives``).

Every candidate gets a full audit row. ``reject_reason`` is the first failing
step in the order above.

The reported attL is one of the 3'-anchored, strand-consistent BLAST hits in
the search window. By default it is the hit with the highest bitscore, so a
long but poor alignment far from the tRNA (often a second copy of the tRNA
gene) does not win over a shorter, near-perfect repeat. ``attl_select_by:
length`` restores the published rule, which takes the longest hit.

The exact run is the longest perfect match of the reported attL against the
tRNA region, starting at the anchored end, with an offset of up to ``shift``
bases to absorb indels.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

DEFAULT_SHIFT = 3
DEFAULT_TRNA_MAX_DISTANCE_BP = 500
DEFAULT_ATTL_WINDOW_BP = 300_000
DEFAULT_ATTL_CANDIDATE_MIN_BP = 8
DEFAULT_ATTL_EXACT_MIN_BP = 14
DEFAULT_ATTL_RELAXED_MIN_BP = 17
DEFAULT_INTEGRASE_MIN_AA = 300
DEFAULT_IE_MIN_NT = 0
DEFAULT_ATTL_SELECT_BY = "bitscore"
ATTL_SELECT_CHOICES = ("bitscore", "length")


@dataclass(frozen=True)
class FilterThresholds:
    """Numeric thresholds applied when selecting confident integrative elements.

    Attributes:
        shift: Max distance (bp) between the BLAST hit end and the tRNA 3' end.
        trna_max_distance_bp: Max distance between integrase and tRNA 3' ends.
        attl_window_bp: Genomic window searched for the attL duplication.
        attl_candidate_min_bp: Minimum length of the exact anchored run.
        attl_exact_min_bp: Exact-run length that passes the length gate.
        attl_relaxed_min_bp: BLAST length that passes the length gate. Mismatches
            are allowed; there is no identity cap.
        attl_reject_gapped: Reject duplications whose alignment length differs
            from the query span or the subject span.
        integrase_min_aa: Integrases of this length and shorter are rejected.
        ie_min_nt: Minimum element length in nucleotides; 0 disables the filter.
        reject_ambiguous_n_ie: Reject elements containing ambiguous N bases.
        attl_select_by: How the reported attL is chosen among anchored hits:
            ``bitscore`` (best alignment) or ``length`` (longest hit, the
            published rule).
    """

    shift: int = DEFAULT_SHIFT
    trna_max_distance_bp: int = DEFAULT_TRNA_MAX_DISTANCE_BP
    attl_window_bp: int = DEFAULT_ATTL_WINDOW_BP
    attl_candidate_min_bp: int = DEFAULT_ATTL_CANDIDATE_MIN_BP
    attl_exact_min_bp: int = DEFAULT_ATTL_EXACT_MIN_BP
    attl_relaxed_min_bp: int = DEFAULT_ATTL_RELAXED_MIN_BP
    attl_reject_gapped: bool = True
    integrase_min_aa: int = DEFAULT_INTEGRASE_MIN_AA
    ie_min_nt: int = DEFAULT_IE_MIN_NT
    reject_ambiguous_n_ie: bool = False
    attl_select_by: str = DEFAULT_ATTL_SELECT_BY


def load_thresholds(config_path: Path | str | None) -> FilterThresholds:
    """Load thresholds from ``ie_finder_config.yaml``; defaults when absent."""
    if config_path is None:
        return FilterThresholds()
    config_path = Path(config_path)
    if not config_path.is_file():
        return FilterThresholds()
    with config_path.open() as fh:
        cfg = yaml.safe_load(fh) or {}
    f = cfg.get("filters", {}) or {}
    select_by = str(f.get("attl_select_by", DEFAULT_ATTL_SELECT_BY)).strip().lower()
    if select_by not in ATTL_SELECT_CHOICES:
        raise ValueError(
            f"filters.attl_select_by must be one of {ATTL_SELECT_CHOICES}, got {select_by!r}"
        )
    return FilterThresholds(
        shift=int(f.get("v3ps_shift", DEFAULT_SHIFT)),
        trna_max_distance_bp=int(f.get("trna_max_distance_bp", DEFAULT_TRNA_MAX_DISTANCE_BP)),
        attl_window_bp=int(f.get("attl_window_bp", DEFAULT_ATTL_WINDOW_BP)),
        attl_candidate_min_bp=int(f.get("attl_candidate_min_bp", DEFAULT_ATTL_CANDIDATE_MIN_BP)),
        attl_exact_min_bp=int(f.get("attl_exact_min_bp", DEFAULT_ATTL_EXACT_MIN_BP)),
        attl_relaxed_min_bp=int(f.get("attl_relaxed_min_bp", DEFAULT_ATTL_RELAXED_MIN_BP)),
        attl_reject_gapped=bool(f.get("attl_reject_gapped", True)),
        integrase_min_aa=int(f.get("integrase_min_aa", DEFAULT_INTEGRASE_MIN_AA)),
        ie_min_nt=int(f.get("ie_min_nt", DEFAULT_IE_MIN_NT)),
        reject_ambiguous_n_ie=bool(cfg.get("reject_ambiguous_n_ie", False)),
        attl_select_by=select_by,
    )


def parse_orfs_gff(path: Path | str) -> dict[str, list[tuple[int, int, str, str]]]:
    """Parse Prodigal CDS features: contig -> [(start, end, strand, cds_id)], 1-based."""
    cds_by_contig: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    path = Path(path)
    if not path.is_file():
        return cds_by_contig
    with path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            cds_id = ""
            for attr in parts[8].split(";"):
                if attr.startswith("ID="):
                    cds_id = attr[3:]
                    break
            cds_by_contig[parts[0]].append((start, end, strand, cds_id))
    return cds_by_contig


def _best_overlap(
    attl_lo: int,
    attl_hi: int,
    cds_list: list[tuple[int, int, str, str]],
) -> tuple[int, tuple[int, int, str, str] | None]:
    best_ov = 0
    best_cds = None
    for c0, c1, cstrand, cid in cds_list:
        ov = max(0, min(attl_hi, c1) - max(attl_lo, c0) + 1)
        if ov > best_ov:
            best_ov = ov
            best_cds = (c0, c1, cstrand, cid)
    return best_ov, best_cds


def classify_prodigal_overlap(
    attl_lo: int,
    attl_hi: int,
    cds_list: list[tuple[int, int, str, str]],
) -> dict[str, Any]:
    """Classify attL vs Prodigal CDS: intergenic / partial_CDS / fully_inside_CDS / no_cds_on_contig."""
    attl_len = attl_hi - attl_lo + 1
    if not cds_list:
        return {
            "prodigal_overlap_bp": 0,
            "prodigal_overlap_frac": 0.0,
            "prodigal_hit_class": "no_cds_on_contig",
            "prodigal_hit_cds_id": "",
            "prodigal_hit_cds_strand": "",
        }
    best_ov, best_cds = _best_overlap(attl_lo, attl_hi, cds_list)
    if best_ov == 0:
        hit_class = "intergenic"
    elif best_ov >= attl_len:
        hit_class = "fully_inside_CDS"
    else:
        hit_class = "partial_CDS"
    return {
        "prodigal_overlap_bp": int(best_ov),
        "prodigal_overlap_frac": float(best_ov) / attl_len if attl_len > 0 else 0.0,
        "prodigal_hit_class": hit_class,
        "prodigal_hit_cds_id": best_cds[3] if best_cds else "",
        "prodigal_hit_cds_strand": best_cds[2] if best_cds else "",
    }


def _oriented_anchored(raw_hits: pd.DataFrame, trna_len: int, trna_strand: str, shift: int) -> pd.DataFrame:
    """Hits anchored at the tRNA 3' end (within ``shift``) and running with the tRNA strand."""
    sub = raw_hits.copy()
    for col in ("qstart", "qend", "sstart", "send", "length", "wstart"):
        if col in sub.columns:
            sub[col] = sub[col].astype(int)
    kept = sub[sub["qend"] >= trna_len - shift + 1]
    if trna_strand == "+":
        kept = kept[kept["sstart"] < kept["send"]]
    else:
        kept = kept[kept["sstart"] > kept["send"]]
    if kept.empty:
        return kept
    kept = kept.copy()
    kept["s_lo"] = kept[["sstart", "send"]].min(axis=1)
    kept["s_hi"] = kept[["sstart", "send"]].max(axis=1)
    kept["abs_lo"] = kept["wstart"] + kept["s_lo"] - 1
    kept["abs_hi"] = kept["wstart"] + kept["s_hi"] - 1
    kept["mismatch_n"] = kept["mismatch"].astype(int) if "mismatch" in kept.columns else 0
    kept["gapopen_n"] = kept["gapopen"].astype(int) if "gapopen" in kept.columns else 0
    return kept


def _hit_record(best: pd.Series, n_hits: int, tier: str) -> dict[str, Any]:
    length = int(best["length"])
    qspan = abs(int(best["qend"]) - int(best["qstart"])) + 1
    sspan = abs(int(best["send"]) - int(best["sstart"])) + 1
    return {
        "best_attl_len_bp": length,
        "attL_abs_lo": int(best["abs_lo"]),
        "attL_abs_hi": int(best["abs_hi"]),
        "attL_strand": "+" if int(best["sstart"]) < int(best["send"]) else "-",
        "attL_mismatch": int(best["mismatch_n"]),
        "attL_gapopen": int(best["gapopen_n"]),
        "attL_has_gap": length != qspan or length != sspan,
        "attL_pident": float(best["pident"]) if "pident" in best.index else float("nan"),
        "attL_bitscore": float(best["bitscore"]) if "bitscore" in best.index else float("nan"),
        "attL_evalue": float(best["evalue"]) if "evalue" in best.index else float("nan"),
        "attL_tier": tier,
        "n_v3ps_hits": int(n_hits),
        "qstart": int(best["qstart"]),
        "qend": int(best["qend"]),
    }


def select_v3ps_strict_hit(
    raw_hits: pd.DataFrame,
    trna_len: int,
    trna_strand: str,
    shift: int = DEFAULT_SHIFT,
    min_len_bp: int = 0,
    select_by: str = DEFAULT_ATTL_SELECT_BY,
) -> dict[str, Any] | None:
    """Best 3'-anchored, strand-consistent hit of at least ``min_len_bp`` (any identity).

    ``select_by="bitscore"`` takes the highest bitscore, then the longer hit on a
    tie. ``select_by="length"`` takes the longest hit, as the published finder did.
    """
    if select_by not in ATTL_SELECT_CHOICES:
        raise ValueError(f"select_by must be one of {ATTL_SELECT_CHOICES}, got {select_by!r}")
    if raw_hits is None or raw_hits.empty:
        return None
    kept = _oriented_anchored(raw_hits, trna_len, trna_strand, shift)
    if min_len_bp:
        kept = kept[kept["length"] >= min_len_bp]
    if kept.empty:
        return None
    if select_by == "length":
        best = kept.loc[kept["length"].idxmax()]
    else:
        ranked = kept.assign(_bits=kept["bitscore"].astype(float)).sort_values(
            ["_bits", "length"], ascending=False, kind="mergesort"
        )
        best = ranked.iloc[0]
    return _hit_record(best, len(kept), "candidate")


def exact_anchored_run(attl: str, trna: str, trna_strand: str, shift: int = DEFAULT_SHIFT) -> int:
    """Longest perfect run of ``attl`` against ``trna`` starting at the anchored end.

    The anchored end is the last base for a ``+`` tRNA and the first base for a
    ``-`` tRNA. An offset of up to ``shift`` bases on either sequence absorbs
    indels before the run starts; the run stops at the first mismatch.
    """
    a, t = attl.upper(), trna.upper()
    if trna_strand == "+":
        a, t = a[::-1], t[::-1]
    best = 0
    for off in range(0, shift + 1):
        for x, y in ((a, t[off:]), (a[off:], t)):
            n = 0
            for ca, ct in zip(x, y):
                if ca != ct:
                    break
                n += 1
            best = max(best, n)
    return best


def select_attl_hit(
    raw_hits: pd.DataFrame,
    trna_len: int,
    trna_strand: str,
    thresholds: FilterThresholds,
) -> dict[str, Any] | None:
    """Best 3'-anchored BLAST hit of at least ``attl_candidate_min_bp``.

    This is the reported attL, chosen by ``thresholds.attl_select_by``. The exact
    run, the length gate and the gap check are all applied to it.
    """
    return select_v3ps_strict_hit(
        raw_hits,
        trna_len,
        trna_strand,
        thresholds.shift,
        thresholds.attl_candidate_min_bp,
        thresholds.attl_select_by,
    )


def integrase_aa_length(start: int, end: int) -> int:
    """Protein length from ORF nucleotide span (1-based inclusive): span // 3.

    The span includes the stop codon, so this is the protein length + 1.
    Kept as is: the published counts were made with this value.
    """
    return (int(end) - int(start) + 1) // 3


def build_qseqid(
    integrase_id: str,
    contig: str,
    trna_start: int,
    trna_end: int,
    trna_strand: str,
) -> str:
    """BLAST query id as written in ``mge_query.fa``: ``integrase_id:contig:start-end:strand``."""
    return f"{integrase_id}:{contig}:{trna_start}-{trna_end}:{trna_strand}"


def evaluate_ie_candidate(
    *,
    integrase_id: str,
    contig: str,
    trna_start: int,
    trna_end: int,
    trna_strand: str,
    trna_len: int,
    integrase_start: int,
    integrase_end: int,
    ie_id: str | None,
    ie_len_nt: int | None,
    ie_has_n: bool,
    raw_hits: pd.DataFrame | None,
    attl_seq: str = "",
    trna_seq: str = "",
    contig_seq: str = "",
    cds_by_contig: dict[str, list[tuple[int, int, str, str]]],
    thresholds: FilterThresholds,
) -> dict[str, Any]:
    """Evaluate one candidate against every filter; ``reject_reason`` is the first failure."""
    aa_len = integrase_aa_length(integrase_start, integrase_end)
    row: dict[str, Any] = {
        "integrase_id": integrase_id,
        "ie_id": ie_id or "",
        "contig": contig,
        "trna_start": trna_start,
        "trna_end": trna_end,
        "trna_strand": trna_strand,
        "trna_len": trna_len,
        "integrase_len_aa": aa_len,
        "ie_len_nt": ie_len_nt if ie_len_nt is not None else 0,
        "ie_has_ambiguous_n": ie_has_n,
        "candidate_attl_len_bp": 0,
        "exact_anchored_bp": 0,
        "best_attl_len_bp": 0,
        "attL_abs_lo": 0,
        "attL_abs_hi": 0,
        "attL_strand": "",
        "attL_mismatch": 0,
        "attL_gapopen": 0,
        "attL_has_gap": False,
        "attL_pident": float("nan"),
        "attL_bitscore": float("nan"),
        "attL_evalue": float("nan"),
        "attL_tier": "",
        "n_v3ps_hits": 0,
        "prodigal_hit_class": "",
        "prodigal_overlap_bp": 0,
        "reject_reason": "",
        "passed_candidate": False,
        "passed_integrase_len": False,
        "passed_intergenic": False,
        "passed_attl_len": False,
        "passed_no_gap": False,
        "passed_no_ambiguous_n": False,
        "passed_ie_len": False,
        "passed_confident": False,
    }
    reasons: list[str] = []

    if raw_hits is None or raw_hits.empty:
        row["reject_reason"] = "v3ps_no_raw_blast"
        return row

    # Reported attL: best 3'-anchored BLAST hit of at least 8 bp (see attl_select_by).
    chosen = select_attl_hit(raw_hits, trna_len, trna_strand, thresholds)
    if chosen is None:
        row["reject_reason"] = "v3ps_no_strict_hit"
        return row
    row.update(chosen)
    row["candidate_attl_len_bp"] = chosen["best_attl_len_bp"]
    # The cut-out island does not always contain the BLAST hit: on a minus-strand
    # tRNA the extract starts at hit_end, so the first bases are not attL.
    # Measure the repeat on the assembly instead.
    if contig_seq:
        lo, hi = int(chosen["attL_abs_lo"]), int(chosen["attL_abs_hi"])
        if 1 <= lo <= hi <= len(contig_seq):
            attl_seq = contig_seq[lo - 1:hi]
        ts, te = int(trna_start), int(trna_end)
        if ts > te:
            ts, te = te, ts
        if 1 <= ts <= te <= len(contig_seq):
            trna_seq = contig_seq[ts - 1:te]
    exact = exact_anchored_run(attl_seq, trna_seq, trna_strand, thresholds.shift) if attl_seq and trna_seq else 0
    row["exact_anchored_bp"] = exact

    # 2. Exact duplication of at least 8 bp.
    row["passed_candidate"] = exact >= thresholds.attl_candidate_min_bp
    if not row["passed_candidate"]:
        reasons.append("exact_run_too_short")

    # 3. Integrase longer than 300 aa.
    row["passed_integrase_len"] = aa_len > thresholds.integrase_min_aa
    if not row["passed_integrase_len"]:
        reasons.append("integrase_too_short")

    # 4. attL in a non-coding region.
    cds_list = cds_by_contig.get(contig, [])
    if not cds_list:
        for key, val in cds_by_contig.items():
            if contig in key or key in contig:
                cds_list = val
                break
    prodigal = classify_prodigal_overlap(chosen["attL_abs_lo"], chosen["attL_abs_hi"], cds_list)
    row.update(prodigal)
    row["passed_intergenic"] = prodigal["prodigal_hit_class"] == "intergenic"
    if not row["passed_intergenic"]:
        reasons.append(f"cds_{prodigal['prodigal_hit_class']}")

    # 5. Length: exact run >= 14, or BLAST length >= 17 with mismatches allowed.
    exact_ok = exact >= thresholds.attl_exact_min_bp
    relaxed_ok = chosen["best_attl_len_bp"] >= thresholds.attl_relaxed_min_bp
    row["passed_attl_len"] = exact_ok or relaxed_ok
    row["attL_tier"] = "exact" if exact_ok else ("relaxed" if relaxed_ok else "")
    if not row["passed_attl_len"]:
        reasons.append("attl_too_short")

    # 6. No alignment gaps.
    row["passed_no_gap"] = (not thresholds.attl_reject_gapped) or (not chosen["attL_has_gap"])
    if not row["passed_no_gap"]:
        reasons.append("attl_gapped")
    row["passed_no_ambiguous_n"] = (not thresholds.reject_ambiguous_n_ie) or (not ie_has_n)
    if not row["passed_no_ambiguous_n"]:
        reasons.append("ambiguous_n")

    # optional element length
    if thresholds.ie_min_nt > 0:
        row["passed_ie_len"] = ie_len_nt is not None and ie_len_nt >= thresholds.ie_min_nt
    else:
        row["passed_ie_len"] = True
    if not row["passed_ie_len"]:
        reasons.append("ie_too_short")

    row["passed_confident"] = not reasons
    row["reject_reason"] = reasons[0] if reasons else ""
    return row
