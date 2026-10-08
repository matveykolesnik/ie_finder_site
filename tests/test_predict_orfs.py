"""Prodigal mode by assembly size."""
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from predict_orfs import prodigal_mode


def _fasta(tmp: Path, *lengths: int) -> Path:
    path = tmp / "genome.fna"
    records = [SeqRecord(Seq("A" * n), id=f"c{i}", description="") for i, n in enumerate(lengths)]
    SeqIO.write(records, path, "fasta")
    return path


class ProdigalModeTest(unittest.TestCase):
    def test_small_input_uses_meta(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(prodigal_mode(_fasta(Path(tmp), 13_000)), "meta")

    def test_total_length_decides(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(prodigal_mode(_fasta(Path(tmp), 12_000, 8_000)), "single")


if __name__ == "__main__":
    unittest.main()
