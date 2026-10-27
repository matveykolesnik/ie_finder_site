"""Genome-coordinate GFF3 for attL/attR."""
import sys
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from export_genome_features import export_genome_features


def _write_inputs(tmp: Path) -> None:
    record = SeqRecord(Seq("A" * 400), id="contig1", description="")
    SeqIO.write([record], tmp / "genome.fna", "fasta")
    pd.DataFrame([
        {
            "integrase_id": "contig1_1",
            "ie_id": "contig1_1:contig1:10-370",
            "contig": "contig1",
            "trna_start": 300,
            "trna_end": 370,
            "trna_strand": "-",
            "attL_abs_lo": 10,
            "attL_abs_hi": 29,
            "attL_strand": "+",
            "attR_abs_lo": 300,
            "attR_abs_hi": 319,
            "passed_confident": True,
        },
        {
            "integrase_id": "contig1_2",
            "ie_id": "contig1_2:contig1:1-50",
            "contig": "contig1",
            "trna_start": 40,
            "trna_end": 50,
            "trna_strand": "+",
            "attL_abs_lo": 1,
            "attL_abs_hi": 20,
            "attL_strand": "+",
            "attR_abs_lo": 41,
            "attR_abs_hi": 50,
            "passed_confident": False,
        },
    ]).to_csv(tmp / "audit.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"orf_id": "contig1_1", "start": 80, "end": 250, "strand": "1"},
    ]).to_csv(tmp / "integrases.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"integrase_id": "contig1_1", "trna_start": 300, "trna_end": 370, "tRNA_type": "tRNA-Ile"},
    ]).to_csv(tmp / "trna.tsv", sep="\t", index=False)


class ExportGenomeFeaturesTest(unittest.TestCase):
    def test_confident_features_use_assembly_coordinates(self):
        import tempfile
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            _write_inputs(tmp)
            n = export_genome_features(
                fasta=tmp / "genome.fna",
                audit=tmp / "audit.tsv",
                integrases=tmp / "integrases.tsv",
                trna=tmp / "trna.tsv",
                out_gff3=tmp / "out.gff3",
                out_tsv=tmp / "out.tsv",
                sample="strain",
                all_candidates=False,
            )
            self.assertEqual(n, 1)
            gff = (tmp / "out.gff3").read_text().splitlines()
            self.assertEqual(gff[0], "##gff-version 3")
            features = [line.split("\t") for line in gff if not line.startswith("#")]
            self.assertEqual(
                [(row[2], int(row[3]), int(row[4]), row[6]) for row in features],
                [
                    ("mobile_genetic_element", 10, 319, "+"),
                    ("attachment_site", 10, 29, "+"),
                    ("attachment_site", 300, 319, "+"),
                    ("CDS", 80, 250, "+"),
                ],
            )
            self.assertIn("Note=attL", features[1][8])
            self.assertIn("Note=attR", features[2][8])
            self.assertIn("product=tRNA-Ile", features[2][8])
            table = pd.read_csv(tmp / "out.tsv", sep="\t")
            self.assertEqual(list(table["integrase_id"]), ["contig1_1"])
            self.assertEqual(int(table.loc[0, "attL_start"]), 10)
            self.assertEqual(table.loc[0, "confidence"], "confident")

    def test_all_candidates_keeps_a_failed_attl(self):
        import tempfile
        with tempfile.TemporaryDirectory() as raw:
            tmp = Path(raw)
            _write_inputs(tmp)
            n = export_genome_features(
                fasta=tmp / "genome.fna",
                audit=tmp / "audit.tsv",
                integrases=tmp / "integrases.tsv",
                trna=None,
                out_gff3=tmp / "out.gff3",
                out_tsv=tmp / "out.tsv",
                sample="strain",
                all_candidates=True,
            )
            self.assertEqual(n, 2)
            table = pd.read_csv(tmp / "out.tsv", sep="\t")
            self.assertEqual(set(table["confidence"]), {"confident", "candidate"})
            # contig1_2 has no integrase hit: its strand is the one opposite the + tRNA.
            strands = dict(zip(table["integrase_id"], table["ie_strand"]))
            self.assertEqual(strands, {"contig1_1": "+", "contig1_2": "-"})


if __name__ == "__main__":
    unittest.main()
