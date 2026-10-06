"""Annotate extracted integrative elements with attachment sites and integrase features."""

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from logger import Logger

logger = Logger(name="annotate_mge", level=Logger.Level.INFO).get_logger()


def reverse_coords(seq_len: int, start: int, end: int) -> tuple[int, int]:
    """Recalculate coordinates after reverse-complementing a sequence.

    Args:
        seq_len: Total sequence length.
        start: Original start coordinate (0-based).
        end: Original end coordinate (0-based, exclusive).

    Returns:
        A tuple (new_start, new_end) after reverse complement.
    """
    return seq_len - end, seq_len - start


def load_input_files(trna_fa_path: str,
                     integrase_path: str,
                     blast_path: str
                     ) -> tuple[dict[str, SeqIO.SeqRecord], pd.DataFrame, pd.DataFrame]:
    """Load input data: tRNA FASTA, integrase hits, and BLAST hits.

    Args:
        trna_fa_path: Path to FASTA file of extracted tRNA regions.
        integrase_path: Path to integrase_hits_summary.tsv.
        blast_path: Path to mge_blast.tsv.

    Returns:
        A tuple of:
            trna_records: dict mapping core_id to SeqRecord,
            integrase_df: DataFrame of integrase hits,
            blast_df: DataFrame of BLAST hits.

    Notes:
        If any input table is empty, returns empty dict/DataFrames.
    """
    try:
        trna_records = {
            rec.id.split(":")[0]: rec
            for rec in SeqIO.parse(trna_fa_path, "fasta")
        }
        integrase_df = pd.read_csv(integrase_path, sep="\t")
        blast_df = pd.read_csv(blast_path, sep="\t")
        return trna_records, integrase_df, blast_df
    except pd.errors.EmptyDataError:
        logger.warning("One or more input tables are empty.")
        return {}, pd.DataFrame(), pd.DataFrame()


def get_annotation_data(record_id: str,
                        trna_records: dict[str, SeqIO.SeqRecord],
                        integrase_df: pd.DataFrame,
                        blast_df: pd.DataFrame
                        ) -> tuple[str, SeqIO.SeqRecord, pd.Series, pd.Series]:
    """Extract tRNA record, integrase hit, and BLAST hit for this region.

    Args:
        record_id: FASTA record ID of the MGE region.
        trna_records: Dict of tRNA SeqRecords keyed by core_id.
        integrase_df: DataFrame of integrase hits.
        blast_df: DataFrame of BLAST hits.

    Returns:
        A tuple (core_id, trna_record, integrase_row, blast_row). If any
        piece is missing, returns (None, None, None, None).
    """
    core_id = record_id.split(":")[0]
    if core_id not in trna_records or core_id not in integrase_df["orf_id"].values:
        return None, None, None, None

    trna_rec = trna_records[core_id]
    int_row = integrase_df[integrase_df["orf_id"] == core_id].iloc[0]
    blast_row = (
        blast_df[blast_df["integrase_id"] == core_id]
        .sort_values("evalue")
        .iloc[0]
    )
    return core_id, trna_rec, int_row, blast_row


def annotate_record(record,
                    trna_rec,
                    int_row: pd.Series,
                    blast_row: pd.Series
                    ) -> tuple:
    """Annotate a single MGE region with att sites and integrase CDS.

    This will reverse-complement the region if the tRNA is on the '+' strand,
    recalculate integrase coordinates, and always label:
      - attL: the BLAST hit region
      - attR: the tRNA region

    Args:
        record: SeqRecord of the MGE region.
        trna_rec: SeqRecord of the extracted tRNA (ID ends with strand).
        int_row: DataFrame row of integrase hit summary.
        blast_row: DataFrame row of BLAST hit info.

    Returns:
        A tuple (updated_record, features_list, att_data_dict).
    """
    features = []

    # Determine strand from tRNA record ID
    trna_strand = trna_rec.id.split(":")[-1]

    # Extract genomic start of MGE region from record ID
    mge_start = int(record.id.split(":")[2].split("-")[0])
    seq_len = len(record.seq)

    # Integrase genomic coordinates
    int_start = int(int_row["start"])
    int_end = int(int_row["end"])
    att_len = int(blast_row["length"])
    trna_len = len(trna_rec.seq)

    # Local integrase coordinates in MGE region
    local_int_start = int_start - mge_start
    local_int_end = int_end - mge_start

    # Reverse-complement if tRNA on '+'
    if trna_strand == "+":
        record.seq = record.seq.reverse_complement()
        new_int_start, new_int_end = (
            seq_len - local_int_end,
            seq_len - local_int_start
        )
        int_strand = -1
    else:
        new_int_start, new_int_end = local_int_start, local_int_end
        int_strand = 1

    # Compute att sites: attL = BLAST hit, attR = tRNA
    if trna_strand == "+":
        # After RC: BLAST at end, tRNA at start
        bl_start = seq_len - att_len
        bl_end = seq_len
        tr_start = 0
        tr_end = trna_len
    else:
        # Normal orientation: BLAST at start, tRNA at end
        bl_start = 0
        bl_end = att_len
        tr_start = seq_len - trna_len
        tr_end = seq_len

    # Add features: attL, attR, integrase CDS
    features.append(
        SeqFeature(FeatureLocation(bl_start, bl_end),
                   type="misc_feature",
                   qualifiers={"note": "attL"})
    )
    features.append(
        SeqFeature(FeatureLocation(tr_start, tr_end),
                   type="misc_feature",
                   qualifiers={"note": "attR"})
    )
    features.append(
        SeqFeature(FeatureLocation(new_int_start, new_int_end,
                                   strand=int_strand),
                   type="CDS",
                   qualifiers={"note": "integrase"})
    )

    # Prepare att table data
    att_data = {
        "attL_start": bl_start,
        "attL_end": bl_end,
        "attR_start": tr_start,
        "attR_end": tr_end,
        "att_length": bl_end - bl_start
    }

    return record, features, att_data


def annotate_mge(fasta_path: str,
                 trna_fa_path: str,
                 integrase_path: str,
                 blast_path: str,
                 output_gbk: str,
                 output_att: str) -> None:
    """Annotate all MGE regions in a FASTA with features and write outputs.

    Args:
        fasta_path: Path to FASTA of MGE regions.
        trna_fa_path: Path to tRNA FASTA (mge_query.fa).
        integrase_path: Path to integrase_hits_summary.tsv.
        blast_path: Path to mge_blast.tsv.
        output_gbk: Path for output GenBank file.
        output_att: Path for output attachment_sites.tsv.
    """
    logger.info("Loading input files")
    trna_records, integrase_df, blast_df = load_input_files(
        trna_fa_path, integrase_path, blast_path
    )

    if not trna_records or integrase_df.empty or blast_df.empty:
        logger.warning("Missing input data; writing empty outputs.")
        open(output_gbk, "w").close()
        open(output_att, "w").close()
        return

    annotated, att_table = [], []
    for record in SeqIO.parse(fasta_path, "fasta"):
        logger.info(f"Annotating {record.id}")
        core_id, trna_rec, int_row, blast_row = get_annotation_data(
            record.id, trna_records, integrase_df, blast_df
        )
        if not core_id:
            logger.warning(f"No data for {record.id.split(':')[0]}; skipping.")
            continue

        record, features, att_data = annotate_record(
            record, trna_rec, int_row, blast_row
        )
        record.features = features
        record.annotations["molecule_type"] = "DNA"
        annotated.append(record)
        att_data["integrase_id"] = core_id
        att_table.append(att_data)

    if annotated:
        logger.info(f"Writing GenBank to {output_gbk}")
        SeqIO.write(annotated, output_gbk, "genbank")
        logger.info(f"Writing attachment table to {output_att}")
        pd.DataFrame(att_table).to_csv(output_att, sep="\t", index=False)
    else:
        logger.warning("No valid annotations; writing empty outputs.")
        open(output_gbk, "w").close()
        open(output_att, "w").close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate MGE regions with att sites and integrase"
    )
    parser.add_argument("--fasta", required=True,
                        help="Input MGE region FASTA")
    parser.add_argument("--trna_fa", required=True,
                        help="Extracted tRNA regions FASTA")
    parser.add_argument("--integrase", required=True,
                        help="Integrase hits summary TSV")
    parser.add_argument("--blast", required=True,
                        help="BLAST hits TSV")
    parser.add_argument("--out_gbk", required=True,
                        help="Output GenBank file")
    parser.add_argument("--out_att", required=True,
                        help="Output attachment sites TSV")
    args = parser.parse_args()

    annotate_mge(
        args.fasta, args.trna_fa, args.integrase, args.blast,
        args.out_gbk, args.out_att
    )