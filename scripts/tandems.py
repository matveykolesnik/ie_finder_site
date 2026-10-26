"""Tandems: elements integrated one after another at the same tRNA.

An element integrates at the 3' end of a tRNA gene and leaves a copy of that end
at its far end (attL). The next element can integrate at the same site, so a
tandem reads, outwards from the tRNA:

    tRNA - element 1 - copy - element 2 - copy ...

The cascade takes attL to be the best-scoring copy, so a tandem comes out one of
two ways, both wrong:

- the outer copy won, and elements 1 and 2 are called as one;
- the inner copy won, and element 2 is left out: its integrase sits next to a
  copy, not a tRNA gene, so it is never paired.

(``filter_confident_ie.resolve_shared_trna`` already resolves the tandems whose
elements were each found from an integrase of their own.)

Here, for each tRNA with confident elements, the copies of its 3' end are found
again, from the tRNA to ``tandem_window_bp`` past the outermost attL, with the
cascade's attL gates (3'-anchored and strand-consistent, the exact-run and length
gates, intergenic, the gap rule) plus ``tandem_min_att_bp`` and
``tandem_min_att_pident``. A copy within ``MIN_SEGMENT_BP`` of another junction is
a duplicated att, not a junction. The junctions -- the tRNA, the elements'
attLs and the copies -- cut the stretch into segments. A segment is an element
when it holds one of the elements' integrases, or an integrase that faces one of
its ends as the cascade requires:

- facing the near copy as it would face a tRNA: opposite strand, 3' ends within
  ``trna_max_distance_bp`` (``integrase_at`` attR);
- with ``integrase_at: either``, facing the far copy from the tRNA's strand, its
  3' end within ``trna_max_distance_bp`` of it (``integrase_at`` attL).

Integrases are the hits the cascade uses, longer than ``integrase_min_aa``. Up to
the outermost attL, a segment that is no element is left with the element before
it (a remnant, or an att-like stretch), so an element is split only where both
parts are elements. Past it, segments are taken while each is an element.

Each member of a tandem is bounded as the cascade bounds an element: from the
copy before it (attR; the tRNA for the first) to its own copy (attL), both
included, so neighbours share a copy; ``attR_after`` names the member before.
Members the cascade did not call are added to the audit with ``found_by =
tandem``; elements cut back to an inner copy get that copy as attL. Every member
carries ``tandem`` (a key per tRNA) and ``tandem_member`` (k/n).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from annotate_mge_region import run_blast_on_region
from logger import get_logger
from v3ps_filters import (
    FilterThresholds,
    _hit_record,
    _oriented_anchored,
    build_qseqid,
    classify_prodigal_overlap,
    exact_anchored_run,
    integrase_aa_length,
)

logger = get_logger("tandems")

MIN_SEGMENT_BP = 1000
# Fields of an audit row that describe its attL, copied with a junction.
ATTL_FIELDS = (
    "best_attl_len_bp", "attL_abs_lo", "attL_abs_hi", "attL_strand", "attL_mismatch",
    "attL_gapopen", "attL_has_gap", "attL_pident", "attL_bitscore", "attL_evalue",
    "attL_tier", "qstart", "qend", "candidate_attl_len_bp", "exact_anchored_bp",
    "prodigal_hit_class", "prodigal_overlap_bp",
)


def _sign(value: object) -> str:
    text = str(value).strip()
    return "+" if text in {"1", "+", "+1"} else "-" if text in {"-1", "-"} else ""


def read_integrases(path: Path, thresholds: FilterThresholds) -> dict[str, list[dict]]:
    """The integrase hits long enough for the cascade, by contig."""
    frame = pd.read_csv(path, sep="\t", dtype={"contig_id": str})
    out: dict[str, list[dict]] = defaultdict(list)
    for _, r in frame.drop_duplicates("orf_id").iterrows():
        aa = integrase_aa_length(int(r["start"]), int(r["end"]), int(r["protein_len_aa"]),
                                 thresholds.integrase_len_from)
        if aa > thresholds.integrase_min_aa:
            out[str(r["contig_id"]).strip()].append({
                "orf_id": str(r["orf_id"]), "start": int(r["start"]), "end": int(r["end"]),
                "strand": _sign(r["strand"]), "aa": aa,
            })
    return out


def find_copies(
    contig_seq: str,
    trna: tuple[int, int],
    strand: str,
    reach: int,
    cds: list,
    thresholds: FilterThresholds,
    tmp_dir: str,
) -> list[dict]:
    """Copies of the tRNA 3' end on its 3' side, out to contig coordinate ``reach``.

    Each copy is an attL-like record (the fields of ``ATTL_FIELDS``) that passes the
    cascade's attL gates and the tandem ones; overlapping hits keep the best.
    """
    lo, hi = trna
    if strand == "+":
        start, end = hi + 1, min(len(contig_seq), reach)
    else:
        start, end = max(1, reach), lo - 1
    if end - start + 1 < thresholds.tandem_min_att_bp:
        return []
    trna_seq = Seq(contig_seq[lo - 1:hi])
    query = SeqRecord(trna_seq if strand == "+" else trna_seq.reverse_complement(), id="trna", description="")
    hits = run_blast_on_region(query, Seq(contig_seq[start - 1:end]), tmp_dir)
    if hits.empty:
        return []
    hits["wstart"] = start
    kept = _oriented_anchored(hits, hi - lo + 1, strand, thresholds.shift)
    copies = []
    for _, h in kept.iterrows():
        rec = _hit_record(h, len(kept), "")
        if rec["best_attl_len_bp"] < thresholds.tandem_min_att_bp or rec["attL_pident"] < thresholds.tandem_min_att_pident:
            continue
        if thresholds.attl_reject_gapped and rec["attL_has_gap"]:
            continue
        seq = contig_seq[rec["attL_abs_lo"] - 1:rec["attL_abs_hi"]]
        exact = exact_anchored_run(seq, contig_seq[lo - 1:hi], strand, thresholds.shift)
        if exact < thresholds.attl_candidate_min_bp:
            continue
        exact_ok = exact >= thresholds.attl_exact_min_bp
        if not (exact_ok or rec["best_attl_len_bp"] >= thresholds.attl_relaxed_min_bp):
            continue
        overlap = classify_prodigal_overlap(rec["attL_abs_lo"], rec["attL_abs_hi"], cds)
        if overlap["prodigal_hit_class"] != "intergenic":
            continue
        rec.update(overlap, exact_anchored_bp=exact, candidate_attl_len_bp=rec["best_attl_len_bp"],
                   attL_tier="exact" if exact_ok else "relaxed")
        copies.append(rec)
    copies.sort(key=lambda r: -r["attL_bitscore"])
    out: list[dict] = []
    for rec in copies:
        if not any(rec["attL_abs_lo"] <= o["attL_abs_hi"] and o["attL_abs_lo"] <= rec["attL_abs_hi"] for o in out):
            out.append(rec)
    return out


def _gap(a: tuple[int, int], b: tuple[int, int]) -> int:
    """Bases between two intervals (0 when they touch or overlap)."""
    return max(0, max(a[0], b[0]) - min(a[1], b[1]) - 1)


def faces_near(integ: dict, junction: tuple[int, int], strand: str, max_distance: int) -> bool:
    """The integrase faces this copy as it would a tRNA: opposite strand, beyond it, 3' ends close."""
    if strand == "+":
        return integ["strand"] == "-" and 0 <= integ["start"] - junction[1] <= max_distance
    return integ["strand"] == "+" and 0 <= junction[0] - integ["end"] <= max_distance


def faces_far(integ: dict, junction: tuple[int, int], strand: str, max_distance: int) -> bool:
    """The integrase faces this copy from the tRNA's strand, before it, 3' end close (attL layout)."""
    if strand == "+":
        return integ["strand"] == "+" and 0 <= junction[0] - integ["end"] - 1 <= max_distance
    return integ["strand"] == "-" and 0 <= integ["start"] - junction[1] - 1 <= max_distance


def _inside(integ: dict, lo: int, hi: int) -> bool:
    return lo <= integ["start"] and integ["end"] <= hi


def _set_attl(row: dict, rec: dict) -> None:
    for key in ATTL_FIELDS:
        if key in rec:
            row[key] = rec[key]


def _set_span(row: dict, integ: tuple[int, int], contig_seq: str, thresholds: FilterThresholds) -> None:
    """The element span from attR, attL and the integrase, and the checks that depend on it."""
    lo = min(int(row["attL_abs_lo"]), int(row["attR_lo"]), *integ)
    hi = max(int(row["attL_abs_hi"]), int(row["attR_hi"]), *integ)
    row["ie_id"] = f"{row['integrase_id']}:{row['contig']}:{lo}-{hi}"
    row["ie_len_nt"] = hi - lo + 1
    row["ie_has_ambiguous_n"] = "N" in contig_seq[lo - 1:hi]
    if thresholds.reject_ambiguous_n_ie and row["ie_has_ambiguous_n"]:
        row["passed_no_ambiguous_n"] = row["passed_confident"] = False
        row["reject_reason"] = "ambiguous_n"
    if thresholds.ie_min_nt > 0 and row["ie_len_nt"] < thresholds.ie_min_nt:
        row["passed_ie_len"] = row["passed_confident"] = False
        row["reject_reason"] = row.get("reject_reason") or "ie_too_short"


def find_tandems(
    rows: list[dict],
    integrases: dict[str, list[dict]],
    int_coords: dict[str, tuple[int, int, int]],
    contig_seqs: dict[str, str],
    cds_by_contig: dict[str, list],
    thresholds: FilterThresholds,
    tmp_dir: str,
) -> list[dict]:
    """Split and add tandem members in the audit rows; returns the added rows.

    ``rows`` are updated in place. Every confident row gets ``tandem``,
    ``tandem_member`` and ``found_by`` (empty unless it is in a tandem or added).
    """
    for row in rows:
        row.setdefault("tandem", "")
        row.setdefault("tandem_member", "")
        row.setdefault("found_by", "")
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("passed_confident"):
            groups[(row["contig"], int(row["trna_start"]), int(row["trna_end"]), row["trna_strand"])].append(row)
    added: list[dict] = []
    for (contig, t_lo, t_hi, strand), chain in groups.items():
        seq = contig_seqs[contig]
        plus = strand == "+"

        def dist(iv: tuple[int, int]) -> int:      # distance from the tRNA, outwards
            return iv[0] - t_hi if plus else t_lo - iv[1]

        chain.sort(key=lambda r: dist((int(r["attL_abs_lo"]), int(r["attL_abs_hi"]))))
        if any(not r.get("attR_after") for r in chain[1:]):
            logger.info(f"{contig}:{t_lo}-{t_hi}: elements nested, not in tandem; left as they are")
            continue
        outer = chain[-1]
        reach = int(outer["attL_abs_hi"]) + thresholds.tandem_window_bp if plus \
            else int(outer["attL_abs_lo"]) - thresholds.tandem_window_bp
        copies = find_copies(seq, (t_lo, t_hi), strand, reach, cds_by_contig.get(contig, []), thresholds, tmp_dir)

        # junctions: the tRNA, the elements' attLs (as they stand now), then the
        # copies that are not a duplicated att of one of them
        junctions = [{"iv": (t_lo, t_hi), "rec": None, "row": None}]
        junctions += [{"iv": (int(r["attL_abs_lo"]), int(r["attL_abs_hi"])), "row": r,
                       "rec": {f: r[f] for f in ATTL_FIELDS if f in r}} for r in chain]
        for rec in sorted(copies, key=lambda r: dist((r["attL_abs_lo"], r["attL_abs_hi"]))):
            iv = (rec["attL_abs_lo"], rec["attL_abs_hi"])
            if all(_gap(iv, j["iv"]) >= MIN_SEGMENT_BP for j in junctions):
                junctions.append({"iv": iv, "rec": rec, "row": None})
        junctions.sort(key=lambda j: dist(j["iv"]))
        if len(junctions) == 2:
            continue                                   # one element, no further copy
        n_inner = max(i for i, j in enumerate(junctions) if j["row"] is not None)

        own = {r["integrase_id"]: r for r in chain}
        integ_here = integrases.get(contig, [])
        segments = []                                  # segments[i - 1] lies between junctions i - 1 and i
        for i in range(1, len(junctions)):
            a, b = junctions[i - 1]["iv"], junctions[i]["iv"]
            lo, hi = (a[0], b[1]) if plus else (b[0], a[1])
            inside = [g for g in integ_here if _inside(g, lo, hi)]
            segments.append({
                "rows": [own[g["orf_id"]] for g in inside if g["orf_id"] in own],
                "near": [g for g in inside if faces_near(g, a, strand, thresholds.trna_max_distance_bp)],
                "far": [g for g in inside if thresholds.integrase_at == "either"
                        and faces_far(g, b, strand, thresholds.trna_max_distance_bp)],
            })

        def owned(i: int) -> bool:
            s = segments[i - 1]
            return bool(s["rows"] or s["near"] or s["far"])

        # members: runs of consecutive segments. Up to the outermost attL a new
        # member starts after an element's attL, or where an integrase faces the
        # copy (or faced it from the segment before) and the member so far is an
        # element; past it, each segment is a member while it is an element.
        members: list[list[int]] = [[1]]
        for i in range(2, len(junctions)):
            if i <= n_inner:
                boundary = junctions[i - 1]["row"] is not None or (
                    any(owned(x) for x in members[-1])
                    and (segments[i - 1]["near"] or segments[i - 2]["far"]))
                if boundary:
                    members.append([i])
                else:
                    members[-1].append(i)
            elif owned(i) and not segments[i - 1]["rows"]:
                members.append([i])
            else:
                break
        # a member that is no element (a remnant) joins the one before it,
        # unless an element's attL separates them
        merged: list[list[int]] = []
        for m in members:
            if merged and not any(owned(i) for i in m) and junctions[m[0] - 1]["row"] is None:
                merged[-1] += m
            else:
                merged.append(m)
        members = merged
        if any(sum(len(segments[i - 1]["rows"]) for i in m) > 1 for m in members) \
                or any(not any(owned(i) for i in m) for m in members):
            logger.info(f"{contig}:{t_lo}-{t_hi}: copies do not delimit one element each; left as they are")
            continue
        if len(members) < 2:
            continue

        # each member's element: an existing row, or an integrase facing its first
        # copy (attR layout) or its last one (attL layout); checked for all members
        # before anything is changed
        plan = []
        for m in members:
            segs = [segments[i - 1] for i in m]
            row_list = [r for s in segs for r in s["rows"]]
            near = segs[0]["near"]
            far = segs[-1]["far"]
            if row_list:
                plan.append((m, row_list[0], None, ""))
            elif near or far:
                plan.append((m, None, near[0] if near else far[0], "attR" if near else "attL"))
            else:
                break
        if len(plan) < len(members):
            logger.info(f"{contig}:{t_lo}-{t_hi}: a member without an integrase at either end; left as it is")
            continue

        key = f"{contig}:{t_lo}-{t_hi}"
        previous = None
        for k, (m, row, integ, at) in enumerate(plan, start=1):
            j_r, j_l = junctions[m[0] - 1], junctions[m[-1]]
            if row is None:
                template = chain[0]
                row = {f: template.get(f) for f in ("sample", "contig", "trna_start", "trna_end", "trna_strand",
                                                    "trna_len")}
                row.update({
                    "integrase_id": integ["orf_id"], "integrase_len_aa": integ["aa"],
                    "qseqid": build_qseqid(integ["orf_id"], contig, t_lo, t_hi, strand),
                    "integrase_at": at, "reject_reason": "", "found_by": "tandem",
                    **{f: True for f in ("passed_candidate", "passed_integrase_len", "passed_intergenic",
                                         "passed_attl_len", "passed_no_gap", "passed_no_ambiguous_n",
                                         "passed_ie_len", "passed_confident")},
                })
                int_coords.setdefault(integ["orf_id"], (integ["start"], integ["end"], integ["aa"]))
                added.append(row)
            if j_l["row"] is not row:              # the attL moves to another copy
                _set_attl(row, j_l["rec"])
            row["attR_lo"], row["attR_hi"] = j_r["iv"]
            row["attR_after"] = previous["integrase_id"] if previous is not None else ""
            _set_span(row, sorted(int_coords[row["integrase_id"]][:2]), seq, thresholds)
            row["tandem"] = key
            row["tandem_member"] = f"{k}/{len(plan)}"
            previous = row
        logger.info(f"{key} ({strand}): tandem of {len(plan)}: "
                    + ", ".join(f"{r['integrase_id']}{' (added)' if r.get('found_by') else ''}"
                                for r in [r for r in rows + added if r.get("tandem") == key]))
    return added
