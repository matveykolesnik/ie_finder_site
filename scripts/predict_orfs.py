"""Predict open reading frames with Prodigal for one genome assembly."""

import argparse
import subprocess
import os

from Bio import SeqIO

from logger import get_logger

log = get_logger("predict_orfs")

# Prodigal's single mode trains on the input itself and stops with an error
# below this many bases in total.
SINGLE_MODE_MIN_BP = 20_000


def prodigal_mode(fna_path):
    """``single`` for an assembly Prodigal can train on, ``meta`` below 20 kb.

    Meta mode uses Prodigal's pre-trained models, so it works on a phage
    genome or a lone plasmid.
    """
    total = sum(len(rec.seq) for rec in SeqIO.parse(fna_path, "fasta"))
    return "single" if total >= SINGLE_MODE_MIN_BP else "meta"


def predict_with_prodigal(fna_path, gff_path, ffn_path, faa_path):
    """Run Prodigal to predict ORFs from a given genome FASTA file.

    This function creates the necessary output directories, constructs the command to run
    Prodigal in single mode (meta mode for inputs under 20 kb), and executes it. It generates three output files:
    a GFF file for gene annotations, an FFN file for nucleotide sequences, and an FAA file for protein sequences.

    Args:
        fna_path (str): Path to the input genome FASTA file.
        gff_path (str): Path to the output GFF file.
        ffn_path (str): Path to the output nucleotide sequence file.
        faa_path (str): Path to the output protein sequence file.

    Raises:
        subprocess.CalledProcessError: If Prodigal execution fails.
    """
    os.makedirs(os.path.dirname(gff_path), exist_ok=True)
    os.makedirs(os.path.dirname(ffn_path), exist_ok=True)
    os.makedirs(os.path.dirname(faa_path), exist_ok=True)

    mode = prodigal_mode(fna_path)
    if mode == "meta":
        log.warning(
            f"{fna_path} has fewer than {SINGLE_MODE_MIN_BP} bp, too little for Prodigal "
            "to train on; using its metagenomic mode"
        )
    cmd = [
        "prodigal",
        "-i", fna_path,
        "-o", gff_path,
        "-d", ffn_path,
        "-a", faa_path,
        "-f", "gff",
        "-p", mode,
    ]

    log.info(f"Running Prodigal:\n{' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        log.info(f"ORF annotation complete for: {fna_path}")
        log.info(f"GFF file created: {gff_path}")
        log.info(f"FFN file created: {ffn_path}")
        log.info(f"FAA file created: {faa_path}")
    except subprocess.CalledProcessError as e:
        log.error(f" Error in Prodigal execution: {e}")
        raise


def main() -> None:
    """Parse command-line arguments and run Prodigal ORF prediction."""
    parser = argparse.ArgumentParser(description="Predict ORFs using Prodigal CLI")
    parser.add_argument("--fna", required=True, help="Path to the input genome FASTA file.")
    parser.add_argument("--gff", required=True, help="Path to the output GFF file.")
    parser.add_argument("--ffn", required=True, help="Path to the output nucleotide sequence file")
    parser.add_argument("--faa", required=True, help="Path to the output protein sequence file")
    args = parser.parse_args()

    predict_with_prodigal(args.fna, args.gff, args.ffn, args.faa)


if __name__ == "__main__":
    main()
