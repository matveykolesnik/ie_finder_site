"""Integrases at the attL end: pairing through attL, and elements sharing a tRNA."""
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from attl_side import attl_side_candidates, attl_window, gap_to_attl
from export_genome_features import export_genome_features
from filter_confident_ie import resolve_shared_trna
from v3ps_filters import FilterThresholds, load_thresholds, thresholds_from_params


def integrase(orf, start, end, strand, contig="c"):
    return {"orf_id": orf, "start": start, "end": end, "strand": strand,
            "contig_id": contig, "contig_length": 1_000_000, "model_accession": "PF00589"}


def trna(start, end, strand, contig="c"):
    return {"contig": contig, "start": start, "end": end, "strand": strand, "tRNA_type": "tRNA-Ile"}


class CandidatesTest(unittest.TestCase):
    """tRNA (attR) ---> ... ---> integrase ---> attL, on one strand."""

    def test_plus_strand_trna_before_the_integrase(self):
        found = attl_side_candidates(
            [integrase("i1", 20000, 21000, "1")],
            [trna(100, 176, "+"), trna(30000, 30076, "+"), trna(5000, 5076, "-")],
            set(), 300_000,
        )
        # the tRNA after the integrase and the one on the other strand do not qualify
        self.assertEqual([t["start"] for t in found["i1"]], [100])

    def test_minus_strand_mirrors_it(self):
        found = attl_side_candidates(
            [integrase("i1", 30000, 31000, "-1")],
            [trna(50000, 50076, "-"), trna(10000, 10076, "-")],
            set(), 300_000,
        )
        self.assertEqual([t["start"] for t in found["i1"]], [50000])

    def test_nearest_trna_first(self):
        found = attl_side_candidates(
            [integrase("i1", 20000, 21000, "+")],
            [trna(100, 176, "+"), trna(5000, 5076, "+")],
            set(), 300_000,
        )
        self.assertEqual([t["start"] for t in found["i1"]], [5000, 100])

    def test_window_and_contig_limit_the_search(self):
        found = attl_side_candidates(
            [integrase("i1", 400000, 401000, "+"), integrase("i2", 20000, 21000, "+", contig="d")],
            [trna(100, 176, "+")],
            set(), 300_000,
        )
        self.assertEqual(found, {})

    def test_integrases_paired_with_a_trna_are_left_out(self):
        found = attl_side_candidates(
            [integrase("i1", 20000, 21000, "+")], [trna(100, 176, "+")], {"i1"}, 300_000,
        )
        self.assertEqual(found, {})


class WindowTest(unittest.TestCase):
    def test_window_starts_past_the_integrase_3_end(self):
        self.assertEqual(attl_window(20000, 21000, "+", 100000, 500, 80), (21001, 21580))
        self.assertEqual(attl_window(30000, 31000, "-", 100000, 500, 80), (29420, 29999))

    def test_window_stops_at_the_contig_ends(self):
        self.assertEqual(attl_window(99000, 99900, "+", 100000, 500, 80), (99901, 100000))
        self.assertEqual(attl_window(300, 1300, "-", 100000, 500, 80), (1, 299))

    def test_gap_between_integrase_and_attl(self):
        self.assertEqual(gap_to_attl({"attL_abs_lo": 21020, "attL_abs_hi": 21064}, 20000, 21000, "+"), 19)
        self.assertEqual(gap_to_attl({"attL_abs_lo": 29936, "attL_abs_hi": 29980}, 30000, 31000, "-"), 19)


def confident(integrase_id, attl, trna_span, strand, at="attR"):
    lo, hi = trna_span
    return {
        "integrase_id": integrase_id, "contig": "c", "trna_start": lo, "trna_end": hi,
        "trna_strand": strand, "attL_abs_lo": attl[0], "attL_abs_hi": attl[1],
        "passed_confident": True, "reject_reason": "", "integrase_at": at,
        "attR_lo": lo, "attR_hi": hi, "attR_after": "",
        "ie_len_nt": 0, "ie_has_ambiguous_n": False,
    }


class SharedTrnaTest(unittest.TestCase):
    """THY_K25-5: an island at tRNA-Ile, then a G2 prophage in the copy it left."""

    seqs = {"c": "A" * 200_000}

    def test_tandem_on_the_minus_strand(self):
        inner = confident("island", (99296, 99342), (133618, 133694), "-")
        outer = confident("g2", (79802, 79846), (133618, 133694), "-", at="attL")
        coords = {"island": (130000, 131000, 400), "g2": (79865, 81052, 395)}
        resolve_shared_trna([outer, inner], coords, self.seqs, FilterThresholds())
        self.assertEqual((outer["attR_lo"], outer["attR_hi"]), (99296, 99342))
        self.assertEqual(outer["attR_after"], "island")
        self.assertEqual(outer["ie_len_nt"], 99342 - 79802 + 1)
        self.assertTrue(outer["passed_confident"])
        self.assertEqual((inner["attR_lo"], inner["attR_hi"]), (133618, 133694))

    def test_tandem_on_the_plus_strand(self):
        inner = confident("a", (20000, 20044), (100, 176), "+")
        outer = confident("b", (45000, 45044), (100, 176), "+", at="attL")
        coords = {"a": (300, 1400, 380), "b": (43000, 44950, 390)}
        resolve_shared_trna([inner, outer], coords, self.seqs, FilterThresholds())
        self.assertEqual((outer["attR_lo"], outer["attR_hi"]), (20000, 20044))
        self.assertEqual(outer["ie_len_nt"], 45044 - 20000 + 1)

    def test_same_attl_twice_is_one_element(self):
        by_trna = confident("int_at_trna", (20000, 20044), (100, 176), "+")
        by_attl = confident("int_at_attl", (20000, 20044), (100, 176), "+", at="attL")
        coords = {"int_at_trna": (300, 1400, 380), "int_at_attl": (18000, 19950, 390)}
        resolve_shared_trna([by_attl, by_trna], coords, self.seqs, FilterThresholds())
        self.assertTrue(by_trna["passed_confident"])
        self.assertFalse(by_attl["passed_confident"])
        self.assertEqual(by_attl["reject_reason"], "duplicate")
        self.assertEqual(by_attl["duplicate_of"], "int_at_trna")

    def test_nested_elements_are_left_alone(self):
        # the further element's integrase sits before the nearer attL
        inner = confident("a", (20000, 20044), (100, 176), "+")
        outer = confident("b", (45000, 45044), (100, 176), "+")
        coords = {"a": (300, 1400, 380), "b": (500, 1600, 380)}
        resolve_shared_trna([inner, outer], coords, self.seqs, FilterThresholds())
        self.assertEqual((outer["attR_lo"], outer["attR_hi"]), (100, 176))
        self.assertEqual(outer["attR_after"], "")

    def test_rejected_rows_are_not_touched(self):
        inner = confident("a", (20000, 20044), (100, 176), "+")
        outer = confident("b", (45000, 45044), (100, 176), "+", at="attL")
        outer["passed_confident"] = False
        coords = {"a": (300, 1400, 380), "b": (43000, 44950, 390)}
        resolve_shared_trna([inner, outer], coords, self.seqs, FilterThresholds())
        self.assertEqual(outer["attR_after"], "")


class ExportTest(unittest.TestCase):
    def test_tandem_and_attl_side_are_marked_and_the_rest_reads_as_before(self):
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            SeqIO.write([SeqRecord(Seq("A" * 200_000), id="c", description="")], tmp / "g.fna", "fasta")
            island = confident("island", (99296, 99342), (133618, 133694), "-")
            g2 = confident("g2", (79802, 79846), (133618, 133694), "-", at="attL")
            g2.update({"attR_lo": 99296, "attR_hi": 99342, "attR_after": "island"})
            for row in (island, g2):
                row["attL_strand"] = "-"
            pd.DataFrame([island, g2]).to_csv(tmp / "audit.tsv", sep="\t", index=False)
            pd.DataFrame([
                {"orf_id": "island", "start": 130000, "end": 131000, "strand": "1"},
                {"orf_id": "g2", "start": 79865, "end": 81052, "strand": "-1"},
            ]).to_csv(tmp / "integrases.tsv", sep="\t", index=False)
            pd.DataFrame([
                {"integrase_id": "island", "trna_start": 133618, "trna_end": 133694, "tRNA_type": "tRNA-Ile"},
                {"integrase_id": "g2", "trna_start": 133618, "trna_end": 133694, "tRNA_type": "tRNA-Ile"},
            ]).to_csv(tmp / "trna.tsv", sep="\t", index=False)
            export_genome_features(
                fasta=tmp / "g.fna", audit=tmp / "audit.tsv", integrases=tmp / "integrases.tsv",
                trna=tmp / "trna.tsv", out_gff3=tmp / "g.gff3", out_gbk=tmp / "g.gbk",
                out_tsv=tmp / "g.tsv", sample="g", all_candidates=False,
            )
            lines = [l.split("\t") for l in (tmp / "g.gff3").read_text().splitlines() if not l.startswith("#")]
        parents = {l[8].split(";")[0]: l for l in lines if l[2] == "mobile_genetic_element"}
        # by position: the G2 prophage first, the island second
        self.assertEqual(parents["ID=IE1"][3:5], ["79802", "99342"])
        self.assertTrue(parents["ID=IE1"][8].endswith(";integrase_at=attL;attR_after=IE2"))
        self.assertEqual(parents["ID=IE2"][8], "ID=IE2;Name=IE2;integrase_id=island;confidence=confident")
        attr = [l for l in lines if l[8].startswith("ID=IE1_attR")][0]
        self.assertEqual(attr[3:5], ["99296", "99342"])
        self.assertIn("product=tRNA-Ile", attr[8])


class ParamsTest(unittest.TestCase):
    def test_default_takes_either_end(self):
        self.assertEqual(FilterThresholds().integrase_at, "either")
        root = Path(__file__).resolve().parents[1]
        self.assertEqual(load_thresholds(root / "search_params.yaml").integrase_at, "either")

    def test_published_rule_and_bad_values(self):
        self.assertEqual(thresholds_from_params({"integrase_at": "tRNA"}).integrase_at, "trna")
        with self.assertRaisesRegex(ValueError, "integrase_at must be one of"):
            thresholds_from_params({"integrase_at": "attl"})


if __name__ == "__main__":
    unittest.main()
