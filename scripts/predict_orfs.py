"""Predict open reading frames with pyrodigal for one genome assembly.

pyrodigal runs Prodigal's gene finder inside Python. With the version pinned
in envs/IE_finder.yaml it calls the same genes, with the same ORF ids and
proteins, as the Prodigal 2.6.3 program the pipeline used before.
"""

import argparse
import os

import pyrodigal
from Bio import SeqIO

from logger import get_logger

log = get_logger("predict_orfs")

# Prodigal's single mode trains on the input itself and refuses fewer bases
# than this in total.
SINGLE_MODE_MIN_BP = 20_000


def prodigal_mode(fna_path):
    """``single`` for an assembly Prodigal can train on, ``meta`` below 20 kb.

    Meta mode uses Prodigal's pre-trained models, so it works on a phage
    genome or a lone plasmid.
    """
    total = sum(len(rec.seq) for rec in SeqIO.parse(fna_path, "fasta"))
    return "single" if total >= SINGLE_MODE_MIN_BP else "meta"


def predict_orfs(fna_path, gff_path, ffn_path, faa_path):
    """Predict genes on every contig and write GFF, gene and protein FASTA files.

    Single mode trains one model on all contigs together, as the Prodigal
    program does; inputs under 20 kb use meta mode. The files have Prodigal's
    layout: ORF ids are ``<contig>_<n>``, and each contig's GFF section starts
    with the ``# Sequence Data`` line that ``hmm_search.py`` reads.

    Args:
        fna_path (str): Path to the input genome FASTA file.
        gff_path (str): Path to the output GFF file.
        ffn_path (str): Path to the output nucleotide sequence file.
        faa_path (str): Path to the output protein sequence file.
    """
    for path in (gff_path, ffn_path, faa_path):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    records = list(SeqIO.parse(fna_path, "fasta"))
    mode = prodigal_mode(fna_path)
    if mode == "meta":
        log.warning(
            f"{fna_path} has fewer than {SINGLE_MODE_MIN_BP} bp, too little for Prodigal "
            "to train on; using its metagenomic mode"
        )
        finder = pyrodigal.GeneFinder(meta=True)
    else:
        finder = pyrodigal.GeneFinder()
        finder.train(*(str(rec.seq) for rec in records))
    log.info(f"Predicting genes with pyrodigal {pyrodigal.__version__} ({mode} mode) on {fna_path}")

    n_genes = 0
    with open(gff_path, "w") as gff, open(ffn_path, "w") as ffn, open(faa_path, "w") as faa:
        for index, rec in enumerate(records):
            genes = finder.find_genes(str(rec.seq))
            # One "##gff-version" line per file, as Prodigal writes it.
            genes.write_gff(gff, sequence_id=rec.id, header=index == 0)
            genes.write_genes(ffn, sequence_id=rec.id)
            genes.write_translations(faa, sequence_id=rec.id)
            n_genes += len(genes)
    log.info(f"{n_genes} genes on {len(records)} contig(s): {gff_path}, {ffn_path}, {faa_path}")


def main() -> None:
    """Parse command-line arguments and predict ORFs."""
    parser = argparse.ArgumentParser(description="Predict ORFs with pyrodigal")
    parser.add_argument("--fna", required=True, help="Path to the input genome FASTA file.")
    parser.add_argument("--gff", required=True, help="Path to the output GFF file.")
    parser.add_argument("--ffn", required=True, help="Path to the output nucleotide sequence file")
    parser.add_argument("--faa", required=True, help="Path to the output protein sequence file")
    args = parser.parse_args()

    predict_orfs(args.fna, args.gff, args.ffn, args.faa)


if __name__ == "__main__":
    main()
