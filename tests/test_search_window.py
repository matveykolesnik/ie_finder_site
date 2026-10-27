"""The attL search window, and the flag for a window that ran into a contig end."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from v3ps_filters import FilterThresholds, attl_search_window, evaluate_ie_candidate

TRNA = "GGGCCCTTTA"
CONTIG = "C" * 100 + TRNA + "C" * 60 + TRNA + "C" * 50       # 230 bp; tRNA 101-110 (+), copy 171-180


class WindowTest(unittest.TestCase):
    def test_plus_strand_runs_past_the_trna_end(self):
        self.assertEqual(attl_search_window(101, 110, "+", 1000, 300), (110, 410, False))
        self.assertEqual(attl_search_window(101, 110, "+", 410, 300), (110, 410, False))
        self.assertEqual(attl_search_window(101, 110, "+", 300, 300), (110, 300, True))

    def test_minus_strand_runs_before_the_trna_start(self):
        self.assertEqual(attl_search_window(701, 710, "-", 1000, 300), (401, 701, False))
        self.assertEqual(attl_search_window(701, 710, "-", 1000, 700), (1, 701, False))
        self.assertEqual(attl_search_window(701, 710, "-", 1000, 800), (1, 701, True))


def evaluate(window_bp: int) -> dict:
    hits = pd.DataFrame([{
        "qseqid": "1_5:1:101-110:+", "sseqid": "subject", "pident": 100.0,
        "length": 10, "mismatch": 0, "gapopen": 0, "qstart": 1, "qend": 10,
        "sstart": 62, "send": 71, "evalue": 1e-3, "bitscore": 20.0, "wstart": 110,
    }])
    return evaluate_ie_candidate(
        integrase_id="1_5", contig="1", trna_start=101, trna_end=110, trna_strand="+",
        trna_len=10, integrase_start=120, integrase_end=165, integrase_protein_aa=15, raw_hits=hits,
        contig_seq=CONTIG, cds_by_contig={}, thresholds=FilterThresholds(attl_window_bp=window_bp),
    )


class AuditFlagTest(unittest.TestCase):
    def test_window_cut_by_the_contig_end_is_flagged(self):
        row = evaluate(300_000)
        self.assertTrue(row["window_truncated"])
        self.assertEqual(row["window_searched_bp"], 121)

    def test_window_inside_the_contig_is_not(self):
        row = evaluate(100)
        self.assertFalse(row["window_truncated"])
        self.assertEqual(row["window_searched_bp"], 101)

    def test_candidate_without_hits_is_flagged_too(self):
        row = evaluate_ie_candidate(
            integrase_id="1_5", contig="1", trna_start=101, trna_end=110, trna_strand="+",
            trna_len=10, integrase_start=120, integrase_end=165, integrase_protein_aa=15,
            raw_hits=pd.DataFrame(), contig_seq=CONTIG, cds_by_contig={}, thresholds=FilterThresholds(),
        )
        self.assertEqual(row["reject_reason"], "v3ps_no_raw_blast")
        self.assertTrue(row["window_truncated"])


if __name__ == "__main__":
    unittest.main()
