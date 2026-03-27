"""Shared fixtures for integration tests in the genomePAM pipeline."""

from pathlib import Path
from typing import Dict

import pytest


def make_fasta_padding(length: int) -> str:
    """Create deterministic padding sequence of repeated ACGT.

    Args:
        length: Number of bases to generate.

    Returns:
        String of ACGT repeated (deterministically, not random).
    """
    pattern = "ACGT"
    return (pattern * ((length // len(pattern)) + 1))[:length]


@pytest.fixture
def synthetic_reference_genome(tmp_path: Path) -> Path:
    """Create a tiny FASTA file with a known Cas9 target at a known position.

    The target sequence GTGAGCCACTGTGCCTGGCC is placed at position 50 (0-based)
    with 50 deterministic bases before and after it on chromosome chr1.

    pyfaidx will write a .fai index file alongside the FASTA during read.

    Args:
        tmp_path: pytest fixture providing a temporary directory.

    Returns:
        Path to the FASTA file.
    """
    fasta_path = tmp_path / "reference.fasta"

    prefix = make_fasta_padding(50)
    target = "GTGAGCCACTGTGCCTGGCC"
    suffix = make_fasta_padding(50)

    fasta_content = f">chr1\n{prefix}{target}{suffix}\n"
    fasta_path.write_text(fasta_content)

    return fasta_path


@pytest.fixture
def target_sequence() -> str:
    """Return the Cas9 target sequence used in synthetic data.

    Returns:
        The 20bp Cas9 target sequence.
    """
    return "GTGAGCCACTGTGCCTGGCC"


@pytest.fixture
def synthetic_sam(tmp_path: Path, synthetic_reference_genome: Path) -> Path:
    """Create a SAM file with reads in umi-tools format.

    The SAM contains multiple reads at the on-target position (~50 on chr1) with
    umi-tools format read names (READ_XXXXXXXX), proper SAM header matching the
    reference genome, and appropriate flags and quality scores for filtering tests.

    Includes:
    - Reads with MAPQ >= 50 (quality control passed)
    - Reads with low MAPQ (30) for filtering
    - Reads without valid UMI suffix for edge cases
    - Both forward and reverse strand reads with appropriate TLEN values

    Args:
        tmp_path: pytest fixture providing a temporary directory.
        synthetic_reference_genome: Path to reference FASTA (provides chr1 length).

    Returns:
        Path to the SAM file.
    """
    sam_path = tmp_path / "reads.sam"

    # Get reference genome length by reading the FASTA
    fasta_content = synthetic_reference_genome.read_text()
    lines = fasta_content.strip().split("\n")
    seq = lines[1] if len(lines) > 1 else ""
    genome_length = len(seq)

    # SAM header lines
    sam_lines = [
        "@HD\tVN:1.0\tSO:coordinate",
        f"@SQ\tSN:chr1\tLN:{genome_length}",
    ]

    # Position for on-target reads (where the target sequence is)
    target_pos = 50 + 1  # 1-based SAM coordinate

    # Read sequence and quality (20bp as specified)
    read_seq = "ACGTACGTACGTACGTACGT"
    read_qual = "I" * 20

    # Create paired-end reads where both pairs map to the same ~10bp window.
    # The analyze() function considers:
    #   - For forward reads (TLEN > 0): read_position = position
    #   - For reverse reads (TLEN < 0): read_position = mate_position + |TLEN| - 1
    #
    # Strategy: Create a reverse read that maps back near the forward reads.
    # Forward read 1 at position target_pos (51), mate at target_pos+10 (61)
    # Reverse read at position target_pos+10 (61), mate at target_pos (51)
    #   read_position = 51 + |TLEN| - 1. If TLEN = -10, then read_position = 51 + 10 - 1 = 60
    # So both forward (51) and reverse (60) are within 10bp, same window!
    #
    # Forward strand read at position 51 (second in pair): flag 128
    sam_lines.append(
        "\t".join([
            "READ1_AACCGGTT",
            "128",  # Second in pair
            "chr1",
            str(target_pos),  # Position 51
            "60",
            "20M",
            "=",
            str(target_pos + 10),  # Mate at 61
            "10",  # TLEN = 10 (positive, so forward)
            read_seq,
            read_qual,
        ])
    )

    # Reverse strand read (second in pair): flag 144 (128 + 16)
    # Position 61, mate at 51. When TLEN < 0:
    # read_position = mate_position + |TLEN| - 1 = 51 + 10 - 1 = 60
    # This is within 10bp of position 51, so they'll be grouped in the same window!
    sam_lines.append(
        "\t".join([
            "READ2_TTGGCCAA",
            "144",  # Second in pair + reverse (128 + 16)
            "chr1",
            str(target_pos + 10),  # Position 61
            "60",
            "20M",
            "=",
            str(target_pos),  # Mate at 51
            "-10",  # TLEN = -10 (negative, so reverse)
            read_seq,
            read_qual,
        ])
    )

    # Another forward strand read at similar position (second in pair): flag 128
    sam_lines.append(
        "\t".join([
            "READ3_GGAACCTT",
            "128",  # Second in pair
            "chr1",
            str(target_pos + 3),  # Position 54
            "60",
            "20M",
            "=",
            str(target_pos + 13),  # Mate at 64
            "10",  # TLEN = 10
            read_seq,
            read_qual,
        ])
    )

    # Low MAPQ read (should be filtered out)
    sam_lines.append(
        "\t".join([
            "READ4_TTAACCGG",
            "128",
            "chr1",
            str(target_pos),
            "30",
            "20M",
            "=",
            str(target_pos + 200),
            "200",
            read_seq,
            read_qual,
        ])
    )

    # Read without valid UMI suffix (edge case)
    sam_lines.append(
        "\t".join([
            "READ5",
            "128",
            "chr1",
            str(target_pos - 10),
            "60",
            "20M",
            "=",
            str(target_pos + 190),
            "200",
            read_seq,
            read_qual,
        ])
    )

    sam_content = "\n".join(sam_lines) + "\n"
    sam_path.write_text(sam_content)

    return sam_path


@pytest.fixture
def test_data_dir() -> Path:
    """Return the path to tests/data/fastq/ directory relative to project root.

    Returns:
        Path object pointing to tests/data/fastq/.
    """
    return Path(__file__).resolve().parent.parent / "data" / "fastq"


@pytest.fixture
def annotations() -> Dict[str, str]:
    """Return annotation dict matching what analyze() expects.

    Returns:
        Dictionary with Description, Targetsite, and Sequence keys.
    """
    return {
        "Description": "test_cas9",
        "Targetsite": "cas9_site1",
        "Sequence": "GTGAGCCACTGTGCCTGGCC",
    }
