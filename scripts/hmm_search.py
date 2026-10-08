"""Search the ORFs for tyrosine integrases with HMMER and add their coordinates."""

import argparse
import re
import subprocess

from logger import get_logger

logger = get_logger("hmm_search")

CUTOFF_FLAGS = {"ga": ["--cut_ga"], "none": []}
"""hmmsearch flags per ``integrase_hmm_cutoff``. ``none`` keeps HMMER's
reporting threshold (E-value 10), as the published finder did."""


def run_hmmsearch(faa_path, hmm_path, output_tbl, cutoff="ga"):
    """Search every ORF against the integrase profiles, writing ``--tblout``.

    Args:
        faa_path: Path to the protein FASTA file (``orfs.faa``).
        hmm_path: Path to the HMM file; it may hold several profiles.
        output_tbl: Path for the per-sequence table.
        cutoff: ``ga`` applies each profile's Pfam gathering threshold;
            ``none`` reports every hit with E-value up to 10.

    Raises:
        subprocess.CalledProcessError: When hmmsearch fails.
    """
    # E-values scale with the search space. hmmscan, which the published finder
    # ran, counts the profiles; hmmsearch counts the ORFs. -Z keeps the old
    # E-values, so cutoff "none" reports the same hits as before.
    n_profiles = sum(1 for line in open(hmm_path) if line.startswith("NAME "))
    cmd = [
        "hmmsearch", *CUTOFF_FLAGS[cutoff], "-Z", str(n_profiles),
        "--noali", "--tblout", output_tbl, hmm_path, faa_path,
    ]
    logger.info(f"Running {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def parse_tblout(tbl_path):
    """Best-scoring profile for each ORF in an ``hmmsearch --tblout`` table.

    Columns are: target (ORF id), target accession, query (profile name),
    query accession, full-sequence E-value, score, ...

    Returns:
        dict: ORF id -> (profile accession, score, E-value) of its best hit.
    """
    best = {}
    with open(tbl_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            orf_id, model_acc = parts[0], parts[3]
            evalue, score = float(parts[4]), float(parts[5])
            if orf_id not in best or score > best[orf_id][1]:
                best[orf_id] = (model_acc, score, evalue)
    return best


def parse_faa(faa_path):
    """Parse the FAA file to extract ORF coordinates.

    The FAA file header is expected to have the format:
      >JBKBIM010000001.1_1 # 3 # 1430 # -1 # ID=1_1;partial=10;...
    The function returns a dictionary mapping the ORF ID to a tuple of (start, end, strand).

    Args:
        faa_path (str): Path to the FAA file.

    Returns:
        dict: A dictionary mapping ORF IDs (str) to tuples (start (int), end (int), strand (str)).
    """
    orf_coords = {}
    with open(faa_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                parts = header.split(' # ')
                if len(parts) < 4:
                    logger.error(f"Incorrect FAA header format: {line.strip()}")
                    continue
                orf_id = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    logger.error(f"Error converting coordinates in FAA: {line.strip()}")
                    continue
                strand = parts[3]
                orf_coords[orf_id] = (start, end, strand)
    return orf_coords


def parse_gff(gff_path):
    """Parse the GFF file to extract contig lengths.

    The GFF file is expected to contain a line with the following format:
      # Sequence Data: seqnum=1;seqlen=128375;seqhdr="JBKBIM010000027.1 MAG: Thermus sp. isolate ..."
    This function extracts the sequence length (seqlen) and the contig ID (the first token in seqhdr).

    Args:
        gff_path (str): Path to the GFF file.

    Returns:
        dict: A dictionary mapping contig IDs (str) to their lengths (int).
    """
    contig_lengths = {}
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                m = re.search(r'seqlen=(\d+);seqhdr="([^"]+)"', line)
                if m:
                    seqlen = int(m.group(1))
                    seqhdr = m.group(2)
                    contig_id = seqhdr.split()[0]
                    contig_lengths[contig_id] = seqlen
                else:
                    logger.error(f"Incorrect GFF line: {line.strip()}")
    return contig_lengths


def write_outputs(hits, faa_coords, contig_lengths, summary_file, orfs_file):
    """Write the integrase summary and ORF coordinate files, in ORF order.

    The summary file has the columns: ORF ID, model accession, start, end,
    strand, contig ID, contig length, bitscore, E-value. The ORFs file has:
    ORF ID, start, end.

    Args:
        hits (dict): ORF id -> (model accession, score, E-value).
        faa_coords (dict): ORF id -> (start, end, strand), in FAA order.
        contig_lengths (dict): Contig id -> contig length.
        summary_file (str): Path to the output summary file.
        orfs_file (str): Path to the output ORFs file.
    """
    for orf_id in hits:
        if orf_id not in faa_coords:
            logger.error(f"ORF {orf_id} not found in FAA")
    with open(summary_file, 'w') as summ_f, open(orfs_file, 'w') as orfs_f:
        summ_f.write("orf_id\tmodel_accession\tstart\tend\tstrand\tcontig_id\tcontig_length\tbitscore\tevalue\n")
        orfs_f.write("orf_id\tstart\tend\n")
        for orf_id, (start, end, strand) in faa_coords.items():
            if orf_id not in hits:
                continue
            model, score, evalue = hits[orf_id]
            # Contig ID is the ORF ID with the trailing _<orf_index> suffix removed.
            contig_id = '_'.join(orf_id.split('_')[:-1])
            if contig_id not in contig_lengths:
                logger.error(f"Contig {contig_id} not found in GFF")
                continue
            contig_len = contig_lengths[contig_id]
            summ_f.write(f"{orf_id}\t{model}\t{start}\t{end}\t{strand}\t{contig_id}\t{contig_len}\t{score}\t{evalue:g}\n")
            orfs_f.write(f"{orf_id}\t{start}\t{end}\n")
    logger.info(f"{len(hits)} integrase hit(s) written to {summary_file}")


def main() -> None:
    """Search one sample's ORFs and write the integrase summary tables."""
    parser = argparse.ArgumentParser(
        description="Find integrase ORFs with hmmsearch and add their coordinates"
    )
    parser.add_argument("--faa", required=True, help="Path to FAA file (orfs.faa)")
    parser.add_argument("--gff", required=True, help="Path to GFF file (orfs.gff)")
    parser.add_argument("--out", required=True, help="Path to tblout file (integrase_hits.txt)")
    parser.add_argument("--summary", required=True, help="Output file for integrase summary table")
    parser.add_argument("--orfs", required=True, help="Output file for ORF coordinates")
    parser.add_argument("--hmm", required=True, help="HMM file with the integrase profiles")
    parser.add_argument(
        "--cutoff", choices=sorted(CUTOFF_FLAGS), default="ga",
        help="ga: Pfam gathering thresholds; none: every hit with E-value up to 10",
    )
    args = parser.parse_args()

    run_hmmsearch(args.faa, args.hmm, args.out, args.cutoff)
    write_outputs(
        parse_tblout(args.out), parse_faa(args.faa), parse_gff(args.gff), args.summary, args.orfs
    )


if __name__ == "__main__":
    main()
