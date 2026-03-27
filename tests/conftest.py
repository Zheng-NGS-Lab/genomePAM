"""Shared fixtures for the genomePAM umi test suite."""

import sys
import gzip
from pathlib import Path
from typing import List

import pytest

GUIDESEQ_DIR = str(Path(__file__).resolve().parent.parent / "modules" / "guideseq")


@pytest.fixture(autouse=True)
def guideseq_on_path() -> None:
    """Add modules/guideseq to sys.path so `from umi import ...` works."""
    if GUIDESEQ_DIR not in sys.path:
        sys.path.insert(0, GUIDESEQ_DIR)


def write_fastq_record(
    header: str, sequence: str, quality: str
) -> str:
    """Format a single FASTQ record as a string.

    Args:
        header: Read identifier (without leading @).
        sequence: Nucleotide sequence.
        quality: Phred+33 quality string (same length as sequence).

    Returns:
        Four-line FASTQ record with trailing newlines.
    """
    return f"@{header}\n{sequence}\n+\n{quality}\n"


def write_fastq_file(
    path: Path, records: List[str], gzipped: bool = False
) -> Path:
    """Write FASTQ records to a file, optionally gzipped.

    Args:
        path: Output file path.
        records: List of formatted FASTQ record strings.
        gzipped: If True, write as gzip-compressed file.

    Returns:
        The path written to.
    """
    content = "".join(records)
    if gzipped:
        with gzip.open(path, "wb") as fh:
            fh.write(content.encode("UTF-8"))
    else:
        path.write_text(content)
    return path


def make_index_sequence(barcode_half: str, pad: str = "N") -> str:
    """Build a 16-base index sequence with barcode at positions 1-7.

    The demultiplex module extracts seq[1:8], so position 0 is padding.

    Args:
        barcode_half: 7-base barcode substring.
        pad: Character for the leading padding position.

    Returns:
        16-character index sequence string.
    """
    return pad + barcode_half + "N" * 8


def make_index2_sequence(molecular_barcode: str = "TTTTGGGG") -> str:
    """Build index2 sequence with molecular barcode at positions 8-15.

    get_umi() extracts i2[1][8:16] as the molecular barcode.

    Args:
        molecular_barcode: 8-base molecular barcode.

    Returns:
        16-character index2 sequence.
    """
    return "N" * 8 + molecular_barcode
