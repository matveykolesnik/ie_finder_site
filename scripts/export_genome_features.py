#!/usr/bin/env python3
"""Write attL/attR on the original assembly, as GFF3.

This script reads the audit table, where attL and the tRNA are in 1-based
inclusive contig coordinates, and writes features on the assembly itself.

Default output is the confident set (the manuscript cascade). ``--all-candidates``
also writes candidates that have an attL coordinate but failed a later filter.

With the tandem search on (``tandems.py``), the members of a tandem carry
``tandem=TA<n>`` and ``tandem_member=k/n`` (and ``found_by=tandem`` when the
cascade alone did not call them), and every element is also written as a
sequence of its own, FASTA and GenBank, sorted into single elements, tandems
whole and tandem members (``write_sequences``). The sequences are read outwards
from the tRNA, so an element at a minus-strand tRNA is reverse-complemented; the
headers give the contig coordinates.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from logger import get_logger

logger = get_logger("export_genome_features")

SOURCE = "IE_finder"
TSV_COLUMNS = [
    "sample",
    "integrase_id",
    "contig",
    "ie_start",
    "ie_end",
    "ie_strand",
    "attL_start",
    "attL_end",
    "attL_strand",
    "attR_start",
    "attR_end",
    "attR_strand",
    "attR_type",
    "trna_strand",
    "integrase_start",
    "integrase_end",
    "integrase_strand",
    "confidence",
    "integrase_at",
    "attR_after",
]
TANDEM_COLUMNS = ["element", "tandem", "tandem_member", "found_by"]


def as_bool(value: object) -> bool:
    """Interpret audit flags written by pandas or as text."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def gff_strand(value: object) -> str:
    """Map Prodigal and audit strand labels onto a GFF3 strand column."""
    text = str(value).strip()
    if text in {"+", "1", "+1"}:
        return "+"
    if text in {"-", "-1"}:
        return "-"
    return "."


def opposite_strand(strand: str) -> str:
    """The other GFF3 strand; ``.`` stays ``.``."""
    return {"+": "-", "-": "+"}.get(strand, ".")


def strand_int(value: object) -> int | None:
    """Map a strand label onto a Biopython strand integer."""
    text = gff_strand(value)
    if text == "+":
        return 1
    if text == "-":
        return -1
    return None


def gff_escape(text: object) -> str:
    """Escape a GFF3 attribute value."""
    return (
        str(text)
        .replace("%", "%25")
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace(",", "%2C")
        .replace("\t", "%09")
        .replace("\n", "%0A")
        .replace("\r", "%0D")
    )


def ordered_span(start: int, end: int) -> tuple[int, int]:
    """Return a 1-based inclusive interval with the lower coordinate first."""
    return (start, end) if start <= end else (end, start)


def _read_table(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype={"contig": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _int_field(row: pd.Series, column: str) -> int | None:
    if column not in row.index or pd.isna(row[column]):
        return None
    try:
        value = int(row[column])
    except (TypeError, ValueError):
        return None
    return value if value > 0 else None


def load_integrases(path: Path) -> dict[str, pd.Series]:
    """Index integrase nucleotide coordinates by ORF id."""
    frame = _read_table(path)
    if frame.empty or "orf_id" not in frame.columns:
        return {}
    frame["orf_id"] = frame["orf_id"].astype(str)
    return {str(row["orf_id"]): row for _, row in frame.drop_duplicates("orf_id").iterrows()}


def load_trna_types(path: Path | None) -> dict[tuple[str, int, int], str]:
    """Map (integrase id, tRNA start, tRNA end) to the Aragorn amino-acid label."""
    if path is None:
        return {}
    frame = _read_table(path)
    if frame.empty or "integrase_id" not in frame.columns:
        return {}
    type_col = "tRNA_type" if "tRNA_type" in frame.columns else None
    if type_col is None:
        return {}
    out: dict[tuple[str, int, int], str] = {}
    for _, row in frame.iterrows():
        start, end = _int_field(row, "trna_start"), _int_field(row, "trna_end")
        if start is None or end is None:
            continue
        lo, hi = ordered_span(start, end)
        out[(str(row["integrase_id"]), lo, hi)] = str(row[type_col])
    return out


def select_rows(audit: pd.DataFrame, all_candidates: bool) -> pd.DataFrame:
    """Keep confident elements, plus attL-bearing candidates when requested."""
    if audit.empty or "integrase_id" not in audit.columns:
        return pd.DataFrame()
    kept: list[pd.Series] = []
    for _, row in audit.iterrows():
        passed = "passed_confident" in row.index and as_bool(row["passed_confident"])
        has_attl = _int_field(row, "attL_abs_lo") is not None and _int_field(row, "attL_abs_hi") is not None
        if passed or (all_candidates and has_attl):
            kept.append(row)
    if not kept:
        return pd.DataFrame()
    return pd.DataFrame(kept)


def build_elements(
    audit: pd.DataFrame,
    integrases: dict[str, pd.Series],
    trna_types: dict[tuple[str, int, int], str],
    records: dict[str, SeqRecord],
    *,
    sample: str,
    all_candidates: bool,
) -> list[dict]:
    """Collect one feature bundle per element, in contig coordinates."""
    chosen = select_rows(audit, all_candidates)
    elements: list[dict] = []
    for _, row in chosen.iterrows():
        contig = str(row.get("contig", "")).strip()
        if contig not in records:
            raise SystemExit(f"Contig {contig!r} from the audit is not in the FASTA")
        attl_lo, attl_hi = _int_field(row, "attL_abs_lo"), _int_field(row, "attL_abs_hi")
        attr_lo, attr_hi = _int_field(row, "attR_abs_lo"), _int_field(row, "attR_abs_hi")
        trna_lo, trna_hi = _int_field(row, "trna_start"), _int_field(row, "trna_end")
        if None in (attl_lo, attl_hi, attr_lo, attr_hi, trna_lo, trna_hi):
            logger.warning(f"Skipping {row.get('integrase_id')}: incomplete att coordinates")
            continue
        attl_start, attl_end = ordered_span(attl_lo, attl_hi)
        # attR is the tRNA bases aligned to attL or, in a tandem, the attL copy
        # of the element before (resolve_shared_trna, tandems.py); the tRNA type
        # labels both.
        attr_start, attr_end = ordered_span(attr_lo, attr_hi)
        trna_start, trna_end = ordered_span(trna_lo, trna_hi)
        integrase_id = str(row["integrase_id"])
        hit = integrases.get(integrase_id)
        int_start = int_end = None
        int_strand = "."
        if hit is not None:
            int_start_raw, int_end_raw = _int_field(hit, "start"), _int_field(hit, "end")
            if int_start_raw is not None and int_end_raw is not None:
                int_start, int_end = ordered_span(int_start_raw, int_end_raw)
            if "strand" in hit.index:
                int_strand = gff_strand(hit["strand"])
        # The element and both att sites take the integrase strand. Without a
        # hit the tRNA gives the same answer: an integrase next to the tRNA sits
        # opposite it, one at attL on its strand.
        ie_strand = int_strand
        if ie_strand == ".":
            ie_strand = gff_strand(row.get("trna_strand", "."))
            if _text_field(row, "integrase_at") != "attL":
                ie_strand = opposite_strand(ie_strand)
        bounds = [attl_start, attl_end, attr_start, attr_end]
        if int_start is not None and int_end is not None:
            bounds.extend([int_start, int_end])
        span = (min(bounds), max(bounds))
        passed = "passed_confident" in row.index and as_bool(row["passed_confident"])
        elements.append({
            "sample": sample,
            "integrase_id": integrase_id,
            "contig": contig,
            "ie_start": span[0],
            "ie_end": span[1],
            "ie_strand": ie_strand,
            "attL_start": attl_start,
            "attL_end": attl_end,
            "attL_strand": ie_strand,
            "attR_start": attr_start,
            "attR_end": attr_end,
            "attR_strand": ie_strand,
            "attR_type": trna_types.get((integrase_id, trna_start, trna_end), ""),
            "trna_strand": gff_strand(row.get("trna_strand", ".")),
            "integrase_start": int_start or 0,
            "integrase_end": int_end or 0,
            "integrase_strand": int_strand,
            "confidence": "confident" if passed else "candidate",
            "integrase_at": _text_field(row, "integrase_at") or "attR",
            "attR_after": _text_field(row, "attR_after"),
            "tandem": _text_field(row, "tandem"),
            "tandem_member": _text_field(row, "tandem_member"),
            "found_by": _text_field(row, "found_by"),
        })
    elements.sort(key=lambda item: (item["contig"], item["ie_start"], item["integrase_id"]))
    for index, item in enumerate(elements, start=1):
        item["name"] = item["element"] = f"IE{index}"
    # A member the tandem search added has no integrase-tRNA pair to name the
    # tRNA; it is the tRNA of the tandem's other members.
    site_type = {item["tandem"]: item["attR_type"] for item in elements if item["tandem"] and item["attR_type"]}
    for item in elements:
        if item["tandem"] and not item["attR_type"]:
            item["attR_type"] = site_type.get(item["tandem"], "")
    # tandems.py keys a tandem by its tRNA; name them TA1, TA2, ... in element order.
    tandem_names: dict[str, str] = {}
    for item in elements:
        if item["tandem"]:
            item["tandem"] = tandem_names.setdefault(item["tandem"], f"TA{len(tandem_names) + 1}")
    # attR_after names the element before by its integrase; give its name instead.
    names = {item["integrase_id"]: item["name"] for item in elements}
    for item in elements:
        if item["attR_after"]:
            item["attR_after"] = names.get(item["attR_after"], "")
    return elements


def _text_field(row: pd.Series, column: str) -> str:
    if column not in row.index or pd.isna(row[column]):
        return ""
    return str(row[column]).strip()


def _extra_notes(item: dict) -> list[str]:
    """What sets an element apart from the usual layout, if anything."""
    notes = []
    if item.get("integrase_at") == "attL":
        notes.append("integrase at attL")
    if item.get("attR_after"):
        notes.append(f"attR is the attL of {item['attR_after']}")
    if item.get("tandem"):
        notes.append(f"member {item['tandem_member']} of tandem {item['tandem']}")
    if item.get("found_by"):
        notes.append(f"found by the {item['found_by']} search")
    return notes


def _gff_line(contig: str, feature_type: str, start: int, end: int, strand: str, attributes: str) -> str:
    return "\t".join([
        contig, SOURCE, feature_type, str(start), str(end), ".", strand, ".", attributes,
    ])


def write_gff3(path: Path, elements: list[dict], records: dict[str, SeqRecord]) -> None:
    """Write a GFF3 file whose sequence ids are the FASTA record ids."""
    lines = ["##gff-version 3"]
    seen: list[str] = []
    for item in elements:
        if item["contig"] not in seen:
            seen.append(item["contig"])
    for contig in seen:
        lines.append(f"##sequence-region {contig} 1 {len(records[contig].seq)}")
    for item in elements:
        name = item["name"]
        parent = (
            f"ID={name};Name={name};integrase_id={gff_escape(item['integrase_id'])};"
            f"confidence={item['confidence']}"
        )
        # Only where they say something, so the usual element reads as before.
        if item.get("integrase_at") == "attL":
            parent += ";integrase_at=attL"
        if item.get("attR_after"):
            parent += f";attR_after={gff_escape(item['attR_after'])}"
        if item.get("tandem"):
            parent += f";tandem={item['tandem']};tandem_member={item['tandem_member']}"
        if item.get("found_by"):
            parent += f";found_by={gff_escape(item['found_by'])}"
        lines.append(_gff_line(
            item["contig"], "mobile_genetic_element",
            item["ie_start"], item["ie_end"], item["ie_strand"], parent,
        ))
        lines.append(_gff_line(
            item["contig"], "attachment_site",
            item["attL_start"], item["attL_end"], item["attL_strand"],
            f"ID={name}_attL;Parent={name};Note=attL",
        ))
        attr_attr = f"ID={name}_attR;Parent={name};Note=attR"
        if item["attR_type"]:
            attr_attr += f";product={gff_escape(item['attR_type'])}"
        lines.append(_gff_line(
            item["contig"], "attachment_site",
            item["attR_start"], item["attR_end"], item["attR_strand"], attr_attr,
        ))
        if item["integrase_start"] and item["integrase_end"]:
            lines.append(_gff_line(
                item["contig"], "CDS",
                item["integrase_start"], item["integrase_end"], item["integrase_strand"],
                f"ID={name}_integrase;Parent={name};Note=integrase;product=tyrosine integrase",
            ))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _gb_feature(feature_type: str, start: int, end: int, strand: str, qualifiers: dict[str, list[str]]) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start - 1, end, strand=strand_int(strand)),
        type=feature_type,
        qualifiers=qualifiers,
    )


def write_table(path: Path, elements: list[dict], tandems: bool = False) -> None:
    """Write a coordinate table, one element per row; with ``tandems``, the tandem columns too."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = TSV_COLUMNS + (TANDEM_COLUMNS if tandems else [])
    frame = pd.DataFrame(elements, columns=columns) if elements else pd.DataFrame(columns=columns)
    frame.to_csv(path, sep="\t", index=False)


def tandem_groups(elements: list[dict]) -> dict[str, list[dict]]:
    """Members of each tandem, in order outwards from the tRNA."""
    groups: dict[str, list[dict]] = {}
    for item in elements:
        if item.get("tandem"):
            groups.setdefault(item["tandem"], []).append(item)
    for members in groups.values():
        members.sort(key=lambda item: int(item["tandem_member"].split("/")[0]))
    return groups


SINGLE_DIR = Path("single_elements")
TANDEM_FULL_DIR = Path("tandem_elements") / "full"
TANDEM_SPLIT_DIR = Path("tandem_elements") / "split"
# Every single element and every tandem member, as GenBank only: one file per element.
ALL_DIR = Path("all")


def _cut(
    name: str,
    record_id: str,
    lo: int,
    hi: int,
    plus: bool,
    parts: list[dict],
    source: SeqRecord,
    description: str,
    part_label: bool,
) -> SeqRecord:
    """The contig span lo..hi read outwards from the tRNA, with each part's features.

    On a minus-strand tRNA the span is reverse-complemented and the features
    turned with it. ``parts`` are elements; each gets a ``mobile_element`` (when
    ``part_label``, i.e. a tandem), its attR and attL and its integrase.
    """

    def local(start: int, end: int) -> tuple[int, int]:
        return (start - lo + 1, end - lo + 1) if plus else (hi - end + 1, hi - start + 1)

    def turned(strand: str) -> str:
        return strand if plus or strand not in "+-" else ("-" if strand == "+" else "+")

    seq = source.seq[lo - 1:hi] if plus else source.seq[lo - 1:hi].reverse_complement()
    features: list[SeqFeature] = []
    for m in parts:
        a, b = local(m["ie_start"], m["ie_end"])
        notes = [m["name"], m["confidence"], f"{m['contig']}:{m['ie_start']}-{m['ie_end']}", *_extra_notes(m)]
        features.append(_gb_feature("mobile_element", a, b, ".", {
            "mobile_element_type": ["integrative element"], "note": notes, "label": [m["name"]]}))
        for site in ("attR", "attL"):
            a, b = local(m[f"{site}_start"], m[f"{site}_end"])
            site_notes = [site] + ([m["attR_type"]] if site == "attR" and m["attR_type"] else [])
            features.append(_gb_feature("misc_feature", a, b, turned(m[f"{site}_strand"]),
                                        {"note": site_notes, "label": [f"{m['name']} {site}"]}))
        if m["integrase_start"] and m["integrase_end"]:
            a, b = local(m["integrase_start"], m["integrase_end"])
            features.append(_gb_feature("CDS", a, b, turned(m["integrase_strand"]), {
                "note": ["integrase"], "product": ["tyrosine integrase"], "label": [f"{m['name']} integrase"]}))
    features.sort(key=lambda f: (int(f.location.start), -int(f.location.end)))
    # GenBank's LOCUS name holds 16 characters: the element or tandem name; the id carries the sample.
    record = SeqRecord(seq, id=record_id, name=name[:16], description=description, features=features)
    record.annotations["molecule_type"] = "DNA"
    return record


def _write_gbk(directory: Path, record: SeqRecord) -> None:
    """One record as <id>.gbk."""
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / f"{record.id}.gbk").open("w") as handle:
        SeqIO.write([record], handle, "genbank")


def _write_pair(directory: Path, record: SeqRecord, header: str) -> None:
    """One record as <id>.gbk and <id>.fna."""
    _write_gbk(directory, record)
    seq = str(record.seq)
    (directory / f"{record.id}.fna").write_text(
        f">{header}\n" + "".join(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80)))


def write_sequences(out_dir: Path, elements: list[dict], records: dict[str, SeqRecord], sample: str) -> tuple[int, int]:
    """Every element as a sequence of its own, sorted by tandem membership.

    Writes, under ``out_dir``, one GenBank and one FASTA file per sequence:

    - ``single_elements/<sample>_<IE>``: the elements that are not in a tandem;
    - ``tandem_elements/full/<sample>_<TA>``: each tandem whole, from its first
      member's attR (in the tRNA) to its last member's attL, members marked;
    - ``tandem_elements/split/<sample>_<IE>``: each member of a tandem alone;
    - ``all/<sample>_<IE>.gbk``: the single elements and the tandem members
      together, GenBank only, the same files as in the two directories above.

    All are read outwards from the tRNA. Returns (single elements, tandems).
    """
    n_single = 0
    for m in elements:
        if m.get("tandem"):
            continue
        plus = m["trna_strand"] != "-"
        where = f"{m['contig']}:{m['ie_start']}-{m['ie_end']}({'+' if plus else '-'})"
        record = _cut(m["name"], f"{sample}_{m['name']}", m["ie_start"], m["ie_end"], plus, [m],
                      records[m["contig"]],
                      f"{sample} integrative element {m['name']} at {m['attR_type'] or 'a tRNA'}, {where}, "
                      f"read outwards from the tRNA", False)
        _write_pair(out_dir / SINGLE_DIR, record,
                    f"{record.id} {where} {m['attR_type']} confidence={m['confidence']}"
                    + "".join(f"; {n}" for n in _extra_notes(m)))
        _write_gbk(out_dir / ALL_DIR, record)
        n_single += 1
    groups = tandem_groups(elements)
    for name, members in groups.items():
        contig = members[0]["contig"]
        plus = members[0]["trna_strand"] != "-"
        lo = min(m["ie_start"] for m in members)
        hi = max(m["ie_end"] for m in members)
        where = f"{contig}:{lo}-{hi}({'+' if plus else '-'})"
        trna = members[0]["attR_type"]
        record = _cut(name, f"{sample}_{name}", lo, hi, plus, members, records[contig],
                      f"{sample} tandem {name} at {trna or 'a tRNA'}, {where}, read outwards from the tRNA: "
                      + ", ".join(m["name"] for m in members), True)
        spans = []
        for m in members:
            a, b = ((m["ie_start"] - lo + 1, m["ie_end"] - lo + 1) if plus
                    else (hi - m["ie_end"] + 1, hi - m["ie_start"] + 1))
            spans.append(f"{m['name']}:{a}-{b}")
        _write_pair(out_dir / TANDEM_FULL_DIR, record, f"{record.id} {where} {trna} members {','.join(spans)}")
        for m in members:
            m_where = f"{m['contig']}:{m['ie_start']}-{m['ie_end']}({'+' if plus else '-'})"
            part = _cut(m["name"], f"{sample}_{m['name']}", m["ie_start"], m["ie_end"], plus, [m],
                        records[contig],
                        f"{sample} integrative element {m['name']}, member {m['tandem_member']} of tandem {name} "
                        f"at {trna or 'a tRNA'}, {m_where}, read outwards from the tRNA", False)
            _write_pair(out_dir / TANDEM_SPLIT_DIR, part,
                        f"{part.id} tandem={sample}_{name} member={m['tandem_member']} {m_where} {trna}"
                        + (f" found_by={m['found_by']}" if m.get("found_by") else ""))
            _write_gbk(out_dir / ALL_DIR, part)
    return n_single, len(groups)


def export_genome_features(
    *,
    fasta: Path,
    audit: Path,
    integrases: Path,
    trna: Path | None,
    out_gff3: Path,
    out_tsv: Path,
    sample: str,
    all_candidates: bool,
    tandems: bool = False,
    out_sequences: Path | None = None,
) -> int:
    """Write GFF3 and a coordinate table; with ``tandems``, the element sequences too
    (``write_sequences``, under ``out_sequences``).

    Returns the element count.
    """
    records = {rec.id: rec for rec in SeqIO.parse(fasta, "fasta")}
    elements = build_elements(
        _read_table(audit),
        load_integrases(integrases),
        load_trna_types(trna),
        records,
        sample=sample,
        all_candidates=all_candidates,
    )
    write_gff3(out_gff3, elements, records)
    write_table(out_tsv, elements, tandems)
    logger.info(f"{sample}: {len(elements)} element(s) -> {out_gff3.name}")
    if tandems and out_sequences is not None:
        n_single, n_tandems = write_sequences(out_sequences, elements, records, sample)
        logger.info(f"{sample}: {n_single} single element(s), {n_tandems} tandem(s) -> {out_sequences}")
    return len(elements)


def main() -> None:
    """Parse arguments and write the output files."""
    parser = argparse.ArgumentParser(
        description="Write attL/attR features in the coordinates of the input assembly"
    )
    parser.add_argument("--fasta", required=True, type=Path, help="Input assembly FASTA")
    parser.add_argument("--audit", required=True, type=Path, help="ie_filter_audit.tsv")
    parser.add_argument("--integrases", required=True, type=Path, help="integrase_hits_summary.tsv")
    parser.add_argument("--trna", type=Path, default=None, help="integrase_trna.tsv")
    parser.add_argument("--out-gff3", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--all-candidates",
        action="store_true",
        help="Also write candidates that have attL coordinates but failed a later filter",
    )
    args = parser.parse_args()
    export_genome_features(
        fasta=args.fasta,
        audit=args.audit,
        integrases=args.integrases,
        trna=args.trna,
        out_gff3=args.out_gff3,
        out_tsv=args.out_tsv,
        sample=args.sample,
        all_candidates=args.all_candidates,
    )


if __name__ == "__main__":
    main()
