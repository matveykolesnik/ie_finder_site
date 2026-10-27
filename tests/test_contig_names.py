"""Contig names with ':' or '|', the characters the BLAST query id and BLAST itself use.

Needs blastn on PATH, as the pipeline does.
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

from annotate_mge_region import main as blast_attl_windows
from annotate_mge_region import run_blast_on_region
from extract_trna_region import extract_trna_sequence
from v3ps_filters import build_qseqid

TRNA = "GGGCGATTAGCTCAGCTGGTCAGAGCGCACGCCTGATAAGCGTGAGGTCGGTGGTTCAAGTCCACCATCGCCCACCA"
ATT = TRNA[-45:]
NAMES = ("contig:1", "gi|12345|ref|NC_000001.1|", "lcl|contig_1")


def contig_seq() -> str:
    """A tRNA at 501-577 (+) and a copy of its 3' end at 1501-1545."""
    rng = random.Random(3)
    seq = [rng.choice("ACGT") for _ in range(3000)]
    for start, piece in ((501, TRNA), (1501, ATT)):
        seq[start - 1:start - 1 + len(piece)] = piece
    return "".join(seq)


class ContigNameTest(unittest.TestCase):
    def test_attl_copy_is_found_whatever_the_contig_is_called(self):
        for name in NAMES:
            with self.subTest(contig=name), tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                genome = tmp / "genome.fna"
                SeqIO.write([SeqRecord(Seq(contig_seq()), id=name, description="")], genome, "fasta")
                pairs = tmp / "integrase_trna.tsv"
                pd.DataFrame([{
                    "integrase_id": f"{name}_1", "model": "PF00589.27", "integrase_start": 600,
                    "integrase_end": 1400, "integrase_strand": -1, "contig": name, "contig_length": 3000,
                    "trna_start": 501, "trna_end": 577, "trna_strand": "+", "tRNA_type": "tRNA-Ile",
                    "distance": 23,
                }]).to_csv(pairs, sep="\t", index=False)
                query = tmp / "mge_query.fa"
                extract_trna_sequence(genome, pairs, query)
                raw = tmp / "mge_blast_raw.tsv"
                blast_attl_windows(str(genome), str(query), str(pairs), str(raw), str(tmp))
                hits = pd.read_csv(raw, sep="\t", dtype={"contig": str, "qseqid": str})
                self.assertEqual(set(hits["qseqid"]), {build_qseqid(f"{name}_1", name, 501, 577, "+")})
                self.assertEqual(set(hits["contig"]), {name})
                self.assertEqual(set(hits["integrase_id"]), {f"{name}_1"})
                copy = hits[(hits["qend"] == len(TRNA)) & (hits["length"] == len(ATT))]
                self.assertEqual(list(copy["wstart"] + copy["sstart"] - 1), [1501])

    def test_unknown_query_stops_the_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = tmp / "genome.fna"
            SeqIO.write([SeqRecord(Seq(contig_seq()), id="c", description="")], genome, "fasta")
            pairs = tmp / "integrase_trna.tsv"
            pd.DataFrame(columns=["integrase_id", "contig", "trna_start", "trna_end", "trna_strand",
                                  "distance"]).to_csv(pairs, sep="\t", index=False)
            query = tmp / "mge_query.fa"
            SeqIO.write([SeqRecord(Seq(TRNA), id="x:c:501-577:+", description="")], query, "fasta")
            with self.assertRaises(ValueError):
                blast_attl_windows(str(genome), str(query), str(pairs), str(tmp / "raw.tsv"), str(tmp))

    def test_several_queries_keep_their_own_ids(self):
        subject = Seq(contig_seq()[1000:2000])
        queries = [SeqRecord(Seq(TRNA), id=build_qseqid(f"{n}_1", n, 501, 577, "+"), description="")
                   for n in NAMES]
        with tempfile.TemporaryDirectory() as tmp:
            hits = run_blast_on_region(queries, subject, tmp)
        self.assertEqual(set(hits["qseqid"]), {q.id for q in queries})


if __name__ == "__main__":
    unittest.main()
