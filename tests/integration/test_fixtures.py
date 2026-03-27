"""Basic tests to verify integration test fixtures are functional."""

from pathlib import Path
from typing import Dict

import pytest

pytestmark = pytest.mark.integration


def test_synthetic_reference_genome(synthetic_reference_genome: Path) -> None:
    """Verify synthetic_reference_genome fixture creates a valid FASTA file.

    Args:
        synthetic_reference_genome: Fixture providing path to generated FASTA.
    """
    assert synthetic_reference_genome.exists(), "FASTA file not created"
    assert synthetic_reference_genome.suffix == ".fasta", "Should be a FASTA file"

    content = synthetic_reference_genome.read_text()
    lines = content.strip().split("\n")

    # Verify header
    assert lines[0] == ">chr1", "Header should be >chr1"

    # Verify sequence contains target
    sequence = lines[1]
    target = "GTGAGCCACTGTGCCTGGCC"
    assert target in sequence, f"Target sequence {target} not found in reference"

    # Verify position: target should start at position 50
    target_start = sequence.find(target)
    assert target_start == 50, f"Target should start at position 50, but found at {target_start}"

    # Verify padding: 50 bases before, target (20bp), 50 bases after = 120bp total
    assert len(sequence) == 120, f"Expected 120bp total, got {len(sequence)}"


def test_target_sequence(target_sequence: str) -> None:
    """Verify target_sequence fixture returns correct value.

    Args:
        target_sequence: Fixture providing the target sequence string.
    """
    assert target_sequence == "GTGAGCCACTGTGCCTGGCC", "Target sequence mismatch"
    assert len(target_sequence) == 20, "Target should be 20bp"


def test_synthetic_sam(synthetic_sam: Path, synthetic_reference_genome: Path) -> None:
    """Verify synthetic_sam fixture creates a valid SAM file.

    Args:
        synthetic_sam: Fixture providing path to generated SAM file.
        synthetic_reference_genome: Fixture providing reference genome path.
    """
    assert synthetic_sam.exists(), "SAM file not created"

    content = synthetic_sam.read_text()
    lines = content.strip().split("\n")

    # Verify SAM header
    assert lines[0].startswith("@HD"), "First line should be @HD header"
    assert lines[1].startswith("@SQ"), "Second line should be @SQ header"
    assert "SN:chr1" in lines[1], "SQ header should reference chr1"

    # Verify we have data lines (not just headers)
    data_lines = [l for l in lines if not l.startswith("@")]
    assert len(data_lines) >= 4, "Should have at least 4 data records"

    # Verify UMI format: check at least one read with proper UMI suffix
    umi_reads = [l for l in data_lines if "_" in l.split("\t")[0]]
    assert len(umi_reads) >= 3, "Should have reads with UMI format (READ_XXXXXXXX)"

    # Verify low MAPQ read exists
    mapq_values = [int(l.split("\t")[4]) for l in data_lines]
    assert 30 in mapq_values, "Should have a read with MAPQ=30 for filtering tests"
    assert any(m >= 60 for m in mapq_values), "Should have high-quality reads with MAPQ>=60"

    # Verify read without UMI suffix exists
    bare_reads = [l for l in data_lines if "_" not in l.split("\t")[0]]
    assert len(bare_reads) >= 1, "Should have at least one read without UMI suffix"


def test_test_data_dir(test_data_dir: Path) -> None:
    """Verify test_data_dir fixture points to correct location.

    Args:
        test_data_dir: Fixture providing path to tests/data/fastq/.
    """
    assert test_data_dir.exists(), "test_data_dir should exist"
    assert test_data_dir.name == "fastq", "Should point to fastq directory"
    assert "data" in str(test_data_dir), "Path should include 'data' directory"


def test_annotations(annotations: Dict[str, str]) -> None:
    """Verify annotations fixture returns correct structure.

    Args:
        annotations: Fixture providing annotation dictionary.
    """
    assert "Description" in annotations, "Missing Description key"
    assert "Targetsite" in annotations, "Missing Targetsite key"
    assert "Sequence" in annotations, "Missing Sequence key"

    assert annotations["Description"] == "test_cas9", "Description mismatch"
    assert annotations["Targetsite"] == "cas9_site1", "Targetsite mismatch"
    assert annotations["Sequence"] == "GTGAGCCACTGTGCCTGGCC", "Sequence mismatch"
