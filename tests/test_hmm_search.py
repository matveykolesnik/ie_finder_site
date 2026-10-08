"""Integrase hits from an hmmsearch table."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from hmm_search import parse_faa, parse_tblout, write_outputs

TBLOUT = """\
# target name  accession  query name       accession   E-value  score  bias
1_52           -          Phage_integrase  PF00589.27  2.1e-33  102.6  0.1
1_713          -          Phage_integrase  PF00589.27  2.4e-25   76.4  0.0
1_52           -          Phage_int_M      PF22022.2   4.7e-10   26.9  0.0
"""


class ParseTbloutTest(unittest.TestCase):
    def test_each_orf_keeps_its_best_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hits.txt"
            path.write_text(TBLOUT)
            hits = parse_tblout(path)
        self.assertEqual(hits["1_52"], ("PF00589.27", 102.6, 2.1e-33))
        self.assertEqual(hits["1_713"][0], "PF00589.27")

    def test_summary_follows_orf_order(self):
        hits = {"1_713": ("PF00589.27", 76.4, 2.4e-25), "1_52": ("PF00589.27", 102.6, 2.1e-33)}
        coords = {"1_52": (100, 1200, "-1", 366), "1_99": (5, 50, "1", 15), "1_713": (7000, 8200, "1", 400)}
        with tempfile.TemporaryDirectory() as tmp:
            summary, orfs = Path(tmp) / "summary.tsv", Path(tmp) / "orfs.tsv"
            write_outputs(hits, coords, {"1": 9000}, summary, orfs)
            lines = summary.read_text().splitlines()
        self.assertEqual(lines[0].split("\t")[-3:], ["bitscore", "evalue", "protein_len_aa"])
        self.assertEqual([line.split("\t")[0] for line in lines[1:]], ["1_52", "1_713"])


class ParseFaaTest(unittest.TestCase):
    def test_protein_length_leaves_out_the_stop(self):
        faa = (
            ">c_1 # 1 # 18 # 1 # ID=1_1;partial=00\nMKV\nLA*\n"
            ">c_2 # 20 # 31 # -1 # ID=1_2;partial=10\nMKVL\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "orfs.faa"
            path.write_text(faa)
            orfs = parse_faa(path)
        self.assertEqual(orfs["c_1"], (1, 18, "1", 5))
        # A gene cut off at the contig edge has no stop codon.
        self.assertEqual(orfs["c_2"], (20, 31, "-1", 4))


if __name__ == "__main__":
    unittest.main()
