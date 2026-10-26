"""Integrases at the attL end of an element, paired with their tRNA through attL.

The published finder pairs an integrase with the tRNA it faces -- opposite
strand, 3' ends within ``trna_max_distance_bp`` -- and looks for attL from that
tRNA. Some elements carry the integrase at the other end instead: on the tRNA's
strand, its 3' end facing attL, the copy of the tRNA 3' end that closes the
element. The tailless Thermus phages of the IN93/P23-77 kind are built so:

    tRNA (attR) ---> ... ---> integrase ---> attL

Nothing next to such an integrase is a tRNA gene, so the first pairing leaves
it alone and the element is never a candidate. Here each integrase left alone
is offered to every tRNA that has it on its 3' side, on the same strand, within
``attl_window_bp``. Those tRNAs are BLASTed, as ``blast_mge`` does, against the
``trna_max_distance_bp`` next to the integrase 3' end; the nearest tRNA with a
3'-anchored, strand-consistent hit there that passes the attL length gates
(the exact run from the anchor, ``attl_exact_min_bp`` or ``attl_relaxed_min_bp``)
becomes the integrase's pair. The gates are applied here, not only later,
because the pair is defined by its attL: a chance match of a few bases to a
nearer tRNA would otherwise take the integrase from the tRNA it copies. The pair
goes into ``integrase_trna.tsv`` with ``integrase_at = attL`` and its hits into
``mge_blast_raw.tsv``, so the cascade in ``filter_confident_ie.py`` treats it
as any other candidate. Pairs made the usual way get ``integrase_at = attR``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio.SeqRecord import SeqRecord

from annotate_mge_region import run_blast_on_region
from annotate_trna_proximity import parse_integrases, parse_trna
from genome import contig_record, read_genome
from logger import get_logger
from v3ps_filters import FilterThresholds, build_qseqid, exact_anchored_run, select_attl_hit

logger = get_logger("attl_side")

PAIR_COLUMNS = [
    "integrase_id", "model", "integrase_start", "integrase_end", "integrase_strand",
    "contig", "contig_length", "trna_start", "trna_end", "trna_strand", "tRNA_type",
    "distance", "integrase_at",
]


def strand_sign(value: object) -> str:
    """'+' or '-' for Prodigal's 1/-1 and for +/-; '' otherwise."""
    text = str(value).strip()
    if text in {"1", "+", "+1"}:
        return "+"
    if text in {"-1", "-"}:
        return "-"
    return ""


def attl_side_candidates(
    integrases: list[dict],
    trnas: list[dict],
    paired: set[str],
    window: int,
) -> dict[str, list[dict]]:
    """For each unpaired integrase, the tRNAs whose 3' side it lies on.

    A tRNA qualifies when it is on the integrase's contig and strand and the
    integrase lies wholly on its 3' side, its 3' end within ``window`` of the
    tRNA 3' end -- the stretch ``blast_mge`` searches for attL. Nearest first.
    """
    by_contig: dict[str, list[dict]] = {}
    for t in trnas:
        by_contig.setdefault(t["contig"].strip(), []).append(t)
    out: dict[str, list[dict]] = {}
    for integ in integrases:
        if integ["orf_id"] in paired:
            continue
        strand = strand_sign(integ["strand"])
        found = []
        for t in by_contig.get(integ["contig_id"].strip(), []):
            if t["strand"] != strand:
                continue
            if strand == "+":
                if t["end"] >= integ["start"]:
                    continue
                dist = integ["end"] - t["end"]
            else:
                if t["start"] <= integ["end"]:
                    continue
                dist = t["start"] - integ["start"]
            if dist <= window:
                found.append((dist, t))
        if found:
            out[integ["orf_id"]] = [t for _, t in sorted(found, key=lambda x: x[0])]
    return out


def attl_window(
    integ_start: int, integ_end: int, strand: str, contig_len: int, max_distance: int, pad: int
) -> tuple[int, int]:
    """1-based inclusive stretch past the integrase 3' end where attL may lie.

    attL must begin within ``max_distance`` of the 3' end; ``pad`` (the longest
    query plus the anchor shift) leaves room for the whole copy.
    """
    if strand == "+":
        return integ_end + 1, min(contig_len, integ_end + max_distance + pad)
    return max(1, integ_start - max_distance - pad), integ_start - 1


def gap_to_attl(chosen: dict, integ_start: int, integ_end: int, strand: str) -> int:
    """Bases between the integrase 3' end and attL; negative if they overlap."""
    if strand == "+":
        return int(chosen["attL_abs_lo"]) - integ_end - 1
    return integ_start - int(chosen["attL_abs_hi"]) - 1


def passes_attl_gates(chosen: dict, contig_seq: str, t: dict, thresholds: FilterThresholds) -> bool:
    """The attL length gates of the cascade, on the forward-strand sequences as there."""
    attl = contig_seq[int(chosen["attL_abs_lo"]) - 1:int(chosen["attL_abs_hi"])]
    trna = contig_seq[t["start"] - 1:t["end"]]
    exact = exact_anchored_run(attl, trna, t["strand"], thresholds.shift)
    if exact < thresholds.attl_candidate_min_bp:
        return False
    return exact >= thresholds.attl_exact_min_bp or chosen["best_attl_len_bp"] >= thresholds.attl_relaxed_min_bp


def trna_record(contig_seq, t: dict, qseqid: str) -> SeqRecord:
    """The tRNA on its own strand, as ``extract_trna_region`` writes it."""
    seq = contig_seq[t["start"] - 1:t["end"]]
    if t["strand"] == "-":
        seq = seq.reverse_complement()
    return SeqRecord(seq, id=qseqid, description="")


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype={"contig": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def find_attl_side_pairs(
    genome_fasta: Path,
    integrases_path: Path,
    trna_path: Path,
    pairs_path: Path,
    raw_path: Path,
    tmp_dir: Path,
    thresholds: FilterThresholds,
) -> int:
    """Pair unpaired integrases through attL; append to the pair and hit tables.

    Returns the number of pairs added.
    """
    pairs = _read(pairs_path)
    if not pairs.empty:
        pairs["integrase_at"] = "attR"
    paired = set(pairs["integrase_id"].astype(str)) if not pairs.empty else set()
    integrases = parse_integrases(str(integrases_path))
    candidates = attl_side_candidates(
        integrases, parse_trna(str(trna_path)), paired, thresholds.attl_window_bp
    )
    genome = read_genome(genome_fasta)
    by_id = {i["orf_id"]: i for i in integrases}
    new_pairs: list[dict] = []
    new_hits: list[pd.DataFrame] = []
    for orf_id, trnas in candidates.items():
        integ = by_id[orf_id]
        contig = integ["contig_id"].strip()
        strand = strand_sign(integ["strand"])
        contig_seq = contig_record(genome, contig).seq
        queries = [
            trna_record(contig_seq, t, build_qseqid(orf_id, contig, t["start"], t["end"], t["strand"]))
            for t in trnas
        ]
        pad = max(len(q.seq) for q in queries) + thresholds.shift
        w0, w1 = attl_window(
            integ["start"], integ["end"], strand, len(contig_seq), thresholds.trna_max_distance_bp, pad
        )
        if w1 < w0:
            continue
        hits = run_blast_on_region(queries, contig_seq[w0 - 1:w1], str(tmp_dir))
        if hits.empty:
            continue
        hits["integrase_id"] = orf_id
        hits["contig"] = contig
        hits["wstart"] = w0
        for t, q in zip(trnas, queries):
            own = hits[hits["qseqid"] == q.id]
            if own.empty:
                continue
            chosen = select_attl_hit(own, t["end"] - t["start"] + 1, t["strand"], thresholds)
            if chosen is None:
                continue
            gap = gap_to_attl(chosen, integ["start"], integ["end"], strand)
            if not 0 <= gap <= thresholds.trna_max_distance_bp:
                continue
            if not passes_attl_gates(chosen, str(contig_seq).upper(), t, thresholds):
                continue
            new_pairs.append({
                "integrase_id": orf_id,
                "model": integ["model_accession"],
                "integrase_start": integ["start"],
                "integrase_end": integ["end"],
                "integrase_strand": integ["strand"],
                "contig": contig,
                "contig_length": integ["contig_length"],
                "trna_start": t["start"],
                "trna_end": t["end"],
                "trna_strand": t["strand"],
                "tRNA_type": t.get("tRNA_type", ""),
                "distance": gap,
                "integrase_at": "attL",
            })
            new_hits.append(own)
            logger.info(
                f"{orf_id}: attL {chosen['attL_abs_lo']}-{chosen['attL_abs_hi']} "
                f"{gap} bp past its 3' end copies {t.get('tRNA_type', '')} "
                f"{t['start']}-{t['end']} {t['strand']}"
            )
            break  # the nearest tRNA with a copy there

    if not pairs.empty or new_pairs:
        out = pd.concat([pairs, pd.DataFrame(new_pairs)], ignore_index=True) if new_pairs else pairs
        out.reindex(columns=PAIR_COLUMNS).to_csv(pairs_path, sep="\t", index=False)
    if new_hits:
        raw = _read(raw_path)
        frames = [raw, *new_hits] if not raw.empty else new_hits
        pd.concat(frames, ignore_index=True).to_csv(raw_path, sep="\t", index=False)
    logger.info(f"{len(new_pairs)} integrase(s) paired through attL")
    return len(new_pairs)
