"""Tandems: elements integrated one after another at one tRNA (tandems.py).

A synthetic contig carries a tRNA and two copies of its 3' end, 10 and 20 kb
away; integrase A faces the tRNA and integrase B faces the first copy. Needs
blastn on PATH, as the pipeline does.
"""
import random
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from export_genome_features import build_elements, write_gff3, write_sequences
from tandems import faces_far, faces_near, find_tandems
from v3ps_filters import FilterThresholds

TRNA = "GGGCGATTAGCTCAGCTGGTCAGAGCGCACGCCTGATAAGCGTGAGGTCGGTGGTTCAAGTCCACCATCGCCCACCA"
ATT = TRNA[-45:]
L = 30000


def build_contig() -> str:
    rng = random.Random(7)
    seq = [rng.choice("ACGTGC") for _ in range(L)]
    for start, piece in ((1001, TRNA), (10001, ATT), (20001, ATT)):
        seq[start - 1:start - 1 + len(piece)] = piece
    return "".join(seq)


PLUS = build_contig()
MINUS = str(Seq(PLUS).reverse_complement())


def flip(lo: int, hi: int) -> tuple[int, int]:
    return L - hi + 1, L - lo + 1


def integrases(plus: bool, with_b: bool = True) -> dict:
    found = [("A", 1150, 2350), ("B", 10120, 11320)][: 2 if with_b else 1]
    out = []
    for orf, a, b in found:
        lo, hi = (a, b) if plus else flip(a, b)
        out.append({"orf_id": orf, "start": lo, "end": hi, "strand": "-" if plus else "+", "aa": 400})
    return {"c": out}


def coords(table: dict) -> dict:
    return {g["orf_id"]: (g["start"], g["end"], g["aa"]) for g in table["c"]}


def cds(table: dict) -> dict:
    return {"c": [(g["start"], g["end"], g["strand"], g["orf_id"]) for g in table["c"]]}


def row(plus: bool, attl: tuple[int, int]) -> dict:
    trna = (1001, 1077) if plus else flip(1001, 1077)
    attr = (1033, 1077) if plus else flip(1033, 1077)     # the tRNA bases aligned to attL (ATT)
    lo, hi = attl if plus else flip(*attl)
    return {
        "sample": "s", "integrase_id": "A", "contig": "c", "trna_start": trna[0], "trna_end": trna[1],
        "trna_strand": "+" if plus else "-", "trna_len": 77, "integrase_len_aa": 400,
        "attL_abs_lo": lo, "attL_abs_hi": hi, "attL_strand": "+" if plus else "-", "best_attl_len_bp": 45,
        "qstart": 33, "qend": 77,
        "passed_confident": True, "reject_reason": "", "integrase_at": "attR",
        "attR_abs_lo": attr[0], "attR_abs_hi": attr[1], "attR_after": "", "ie_id": "", "ie_len_nt": 0,
    }


def run(rows: list[dict], plus: bool, with_b: bool = True) -> list[dict]:
    table = integrases(plus, with_b)
    with tempfile.TemporaryDirectory() as tmp:
        return find_tandems(rows, table, coords(table), {"c": PLUS if plus else MINUS}, cds(table),
                            FilterThresholds(), tmp)


class FacingTest(unittest.TestCase):
    def test_near_copy_is_faced_from_the_opposite_strand(self):
        b = {"start": 10120, "end": 11320, "strand": "-"}
        self.assertTrue(faces_near(b, (10001, 10045), "+", 700))
        self.assertFalse(faces_near(dict(b, strand="+"), (10001, 10045), "+", 700))
        self.assertFalse(faces_near(b, (10001, 10045), "+", 50))

    def test_far_copy_is_faced_from_the_trnas_strand(self):
        g = {"start": 18000, "end": 19800, "strand": "+"}
        self.assertTrue(faces_far(g, (20001, 20045), "+", 700))
        self.assertFalse(faces_far(g, (20001, 20045), "-", 700))


class TandemTest(unittest.TestCase):
    def near(self, got, want, slack=3):
        self.assertTrue(all(abs(x - y) <= slack for x, y in zip(got, want)), f"{got} vs {want}")

    def check_split(self, rows, added, plus):
        a = rows[0]
        self.assertEqual(len(added), 1)
        b = added[0]
        span = (lambda lo, hi: (lo, hi)) if plus else flip
        self.assertEqual(b["integrase_id"], "B")
        self.assertEqual(b["found_by"], "tandem")
        # BLAST may stop an alignment a base or two short of the anchor (v3ps_shift)
        self.near((a["attL_abs_lo"], a["attL_abs_hi"]), span(10001, 10045))
        self.near((b["attL_abs_lo"], b["attL_abs_hi"]), span(20001, 20045))
        # neighbours share the copy between them exactly
        self.assertEqual((b["attR_abs_lo"], b["attR_abs_hi"]), (a["attL_abs_lo"], a["attL_abs_hi"]))
        self.assertEqual(b["attR_after"], "A")
        self.assertEqual((a["tandem_member"], b["tandem_member"]), ("1/2", "2/2"))
        self.assertEqual(a["tandem"], b["tandem"])
        self.assertEqual(b["integrase_at"], "attR")

    def test_two_elements_called_as_one_are_split(self):
        rows = [row(True, (20001, 20045))]
        self.check_split(rows, run(rows, True), True)

    def test_element_past_attl_is_added(self):
        rows = [row(True, (10001, 10045))]
        self.check_split(rows, run(rows, True), True)

    def test_minus_strand_mirrors_it(self):
        rows = [row(False, (20001, 20045))]
        self.check_split(rows, run(rows, False), False)
        rows = [row(False, (10001, 10045))]
        self.check_split(rows, run(rows, False), False)

    def test_no_integrase_at_the_copy_no_tandem(self):
        rows = [row(True, (10001, 10045))]
        self.assertEqual(run(rows, True, with_b=False), [])
        self.assertEqual((rows[0]["attL_abs_lo"], rows[0]["tandem"]), (10001, ""))

    def test_inner_copy_without_integrase_does_not_split(self):
        rows = [row(True, (20001, 20045))]
        self.assertEqual(run(rows, True, with_b=False), [])
        self.assertEqual((rows[0]["attL_abs_lo"], rows[0]["tandem"]), (20001, ""))


class WriteTest(unittest.TestCase):
    """Sequences sorted into single elements, tandems whole and members, read from the tRNA;
    the elements' GenBank files together in all/."""

    def elements(self, rows):
        table = integrases(False)
        hits = {g["orf_id"]: pd.Series({"start": g["start"], "end": g["end"], "strand": g["strand"]})
                for g in table["c"]}
        records = {"c": SeqRecord(Seq(MINUS), id="c")}
        return build_elements(pd.DataFrame(rows), hits, {}, records, sample="s", all_candidates=False), records

    def test_tandem_written_whole_and_split(self):
        rows = [row(False, (20001, 20045))]
        rows += run(rows, False)
        elements, records = self.elements(rows)
        self.assertEqual([e["tandem"] for e in elements], ["TA1", "TA1"])
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertEqual(write_sequences(out, elements, records, "s"), (0, 1))
            self.assertFalse((out / "single_elements").exists())
            full = out / "tandem_elements" / "full"
            self.assertEqual(sorted(p.name for p in full.iterdir()), ["s_TA1.fna", "s_TA1.gbk"])
            tandem = SeqIO.read(full / "s_TA1.fna", "fasta")
            self.assertTrue(str(tandem.seq).startswith(ATT))
            marked = [f for f in SeqIO.read(full / "s_TA1.gbk", "genbank").features if f.type == "mobile_element"]
            self.assertEqual(len(marked), 2)
            self.assertEqual(int(marked[0].location.start), 0)
            split = out / "tandem_elements" / "split"
            first = min(elements, key=lambda e: e["tandem_member"])
            members = {p.stem for p in split.glob("*.fna")}
            self.assertEqual(members, {f"s_{e['name']}" for e in elements})
            seq1 = SeqIO.read(split / f"s_{first['name']}.fna", "fasta")
            self.assertTrue(str(seq1.seq).startswith(ATT))
            shared = rows[0]["attL_abs_hi"] - rows[0]["attL_abs_lo"] + 1
            total = sum(len(SeqIO.read(p, "fasta")) for p in split.glob("*.fna"))
            self.assertEqual(total - len(tandem), shared)
            # all/ holds the members' GenBank files, not the whole tandem
            every = out / "all"
            self.assertEqual(sorted(p.name for p in every.iterdir()), sorted(p.name for p in split.glob("*.gbk")))
            for p in every.iterdir():
                self.assertEqual(p.read_text(), (split / p.name).read_text())

    def test_single_element(self):
        rows = [row(False, (10001, 10045))]
        self.assertEqual(run(rows, False, with_b=False), [])
        elements, records = self.elements(rows)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertEqual(write_sequences(out, elements, records, "s"), (1, 0))
            self.assertFalse((out / "tandem_elements").exists())
            rec = SeqIO.read(out / "single_elements" / "s_IE1.gbk", "genbank")
            self.assertTrue(str(rec.seq).startswith(ATT))
            self.assertEqual(rec.name, "IE1")
            self.assertEqual([p.name for p in (out / "all").iterdir()], ["s_IE1.gbk"])
            self.assertEqual((out / "all" / "s_IE1.gbk").read_text(),
                             (out / "single_elements" / "s_IE1.gbk").read_text())

    # The minus-strand tRNA sits 1 kb from the contig end it faces, so its attL
    # search window was cut there, as the cascade would record.
    CUT = {"window_truncated": True, "window_searched_bp": 28924}

    def test_cut_window_goes_to_truncated_elements(self):
        rows = [{**row(False, (10001, 10045)), **self.CUT}]
        self.assertEqual(run(rows, False, with_b=False), [])
        elements, records = self.elements(rows)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertEqual(write_sequences(out, elements, records, "s"), (1, 0))
            cut = out / "truncated_elements" / "single_elements"
            self.assertEqual(sorted(p.name for p in cut.iterdir()), ["s_IE1.fna", "s_IE1.gbk"])
            self.assertIn("attL search window cut by the contig end", (cut / "s_IE1.fna").read_text())
            self.assertFalse((out / "single_elements").exists())
            self.assertFalse((out / "all").exists())
            write_gff3(out / "s.gff3", elements, records)
            self.assertIn(";window_truncated=true", (out / "s.gff3").read_text())

    def test_cut_window_moves_a_tandem_whole(self):
        rows = [{**row(False, (20001, 20045)), **self.CUT}]
        rows += run(rows, False)
        elements, records = self.elements(rows)
        # the member the tandem search added takes the flag from its tRNA's row
        self.assertEqual([e["window_truncated"] for e in elements], [True, True])
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self.assertEqual(write_sequences(out, elements, records, "s"), (0, 1))
            cut = out / "truncated_elements" / "tandem_elements"
            self.assertEqual(sorted(p.name for p in (cut / "full").iterdir()), ["s_TA1.fna", "s_TA1.gbk"])
            self.assertEqual(len(list((cut / "split").glob("*.fna"))), 2)
            self.assertFalse((out / "tandem_elements").exists())
            self.assertFalse((out / "all").exists())


if __name__ == "__main__":
    unittest.main()
