"""Quality-control helpers for genome and integrative element FASTA files."""

from __future__ import annotations

from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:  # pragma: no cover
    SeqIO = None  # type: ignore


def fasta_contains_ambiguous_n(path: Path | str) -> bool:
    """Return whether any sequence record contains ambiguous IUPAC N bases.

    Args:
        path: Path to a FASTA file containing one or more records.

    Returns:
        True if at least one record contains ``N`` or ``n``; False otherwise.
        Returns True when the file cannot be parsed.
    """
    path = Path(path)
    if SeqIO is None:
        text = path.read_text(encoding="utf-8", errors="replace")
        return "N" in text.upper()
    try:
        for rec in SeqIO.parse(path, "fasta"):
            seq = str(rec.seq)
            if "N" in seq or "n" in seq:
                return True
    except Exception:
        return True
    return False


def parse_fasta_lengths_and_n_flags(
    path: Path | str,
) -> tuple[dict[str, int], dict[str, bool]]:
    """Parse per-record sequence length and ambiguous-base flags from a FASTA file.

    Args:
        path: Path to a FASTA file.

    Returns:
        A tuple ``(lengths, has_ambiguous_n)`` where each dictionary is keyed by
        record identifier. ``lengths[id]`` is the sequence length in nucleotides;
        ``has_ambiguous_n[id]`` is True when the record contains ``N``.
    """
    path = Path(path)
    lengths: dict[str, int] = {}
    has_n: dict[str, bool] = {}
    current_id: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(parts)
                    lengths[current_id] = len(seq)
                    has_n[current_id] = "N" in seq.upper()
                current_id = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if current_id is not None:
        seq = "".join(parts)
        lengths[current_id] = len(seq)
        has_n[current_id] = "N" in seq.upper()
    return lengths, has_n
