"""Unit tests for the fq() FASTQ reader against real downsampled gzipped data.

Exercises the fq() generator from the demultiplex module
against real sequencing data (Cas9 10K reads, FnCas12a 10K reads) to catch
issues that synthetic data cannot: malformed records, encoding quirks,
and off-by-one errors in position arithmetic on real sequencing output.
"""

import re
from itertools import islice
from pathlib import Path
from typing import List

import pytest


VALID_BASES_PATTERN = re.compile(r"^[ACGTN]+$")
EXPECTED_RECORD_COUNT = 10_000
CONSISTENCY_SAMPLE_SIZE = 100

ALL_FASTQ_FILES = [
    ("cas9", "SRR33421097_R1.fastq.gz"),
    ("cas9", "SRR33421097_R2.fastq.gz"),
    ("fncas12a", "GNPAM327_R1.fastq.gz"),
    ("fncas12a", "GNPAM327_R2.fastq.gz"),
]

READ_LENGTH_PARAMS = [
    ("cas9", 151),
    ("fncas12a", 150),
]

FILENAME_MAP = {
    "cas9": ("SRR33421097_R1.fastq.gz", "SRR33421097_R2.fastq.gz"),
    "fncas12a": ("GNPAM327_R1.fastq.gz", "GNPAM327_R2.fastq.gz"),
}


def format_fastq_file_id(parametrize_value: tuple) -> str:
    """Generate readable test ID from (subdir, filename) tuple.

    Args:
        parametrize_value: Tuple of (subdirectory, filename).

    Returns:
        Human-readable test identifier string.
    """
    return f"{parametrize_value[0]}/{parametrize_value[1]}"


@pytest.fixture(scope="module")
def test_data_dir() -> Path:
    """Return path to tests/data/fastq/ directory.

    Returns:
        Path object pointing to the FASTQ test data directory.
    """
    data_directory = Path(__file__).resolve().parent / "data" / "fastq"
    assert data_directory.is_dir(), f"Test data directory not found: {data_directory}"
    return data_directory


class TestFqReaderRealData:
    """Validate fq() reader behavior on real gzipped FASTQ files."""

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_yields_expected_record_count(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """Each file contains exactly 10,000 records."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        record_count = sum(1 for _ in fq(filepath))

        assert record_count == EXPECTED_RECORD_COUNT

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_each_record_is_four_element_list(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """Every record yielded by fq() has exactly 4 lines."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        for record_index, record in enumerate(fq(filepath)):
            assert isinstance(record, list), (
                f"Record {record_index} is not a list"
            )
            assert len(record) == 4, (
                f"Record {record_index} has {len(record)} elements, expected 4"
            )

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_header_line_starts_with_at_symbol(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """Every record's first line starts with '@'."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        for record_index, record in enumerate(fq(filepath)):
            assert record[0].startswith("@"), (
                f"Record {record_index} header does not start with '@': "
                f"{record[0][:50]}"
            )

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_plus_line_starts_with_plus(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """Every record's third line starts with '+'."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        for record_index, record in enumerate(fq(filepath)):
            assert record[2].startswith("+"), (
                f"Record {record_index} plus line does not start with '+': "
                f"{record[2][:50]}"
            )

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_all_lines_are_newline_terminated(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """All 4 lines in every record end with a newline character."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        for record_index, record in enumerate(fq(filepath)):
            for line_index, line in enumerate(record):
                assert line.endswith("\n"), (
                    f"Record {record_index}, line {line_index} missing "
                    f"trailing newline: {line[:50]!r}"
                )

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_sequence_and_quality_have_equal_length(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """Sequence (line 2) and quality (line 4) have equal length in every record."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        for record_index, record in enumerate(fq(filepath)):
            sequence_length = len(record[1].strip())
            quality_length = len(record[3].strip())
            assert sequence_length == quality_length, (
                f"Record {record_index}: sequence length {sequence_length} "
                f"!= quality length {quality_length}"
            )

    @pytest.mark.parametrize(
        "subdirectory,filename", ALL_FASTQ_FILES, ids=format_fastq_file_id
    )
    def test_sequence_contains_only_valid_bases(
        self, test_data_dir: Path, subdirectory: str, filename: str
    ) -> None:
        """All bases in every sequence are in {A, C, G, T, N}."""
        from umi.demultiplex import fq

        filepath = str(test_data_dir / subdirectory / filename)
        for record_index, record in enumerate(fq(filepath)):
            sequence = record[1].strip()
            assert VALID_BASES_PATTERN.match(sequence), (
                f"Record {record_index} contains invalid bases: "
                f"{sequence[:50]}"
            )


class TestFqReaderReadLengths:
    """Validate that read lengths match expected values for each dataset."""

    @pytest.mark.parametrize("dataset,expected_length", READ_LENGTH_PARAMS)
    def test_read_lengths_are_uniform(
        self, test_data_dir: Path, dataset: str, expected_length: int
    ) -> None:
        """All sequences in both R1 and R2 have the expected read length."""
        from umi.demultiplex import fq

        r1_filename, r2_filename = FILENAME_MAP[dataset]

        for read_label, filename in [("R1", r1_filename), ("R2", r2_filename)]:
            filepath = str(test_data_dir / dataset / filename)
            for record_index, record in enumerate(fq(filepath)):
                sequence_length = len(record[1].strip())
                assert sequence_length == expected_length, (
                    f"{dataset} {read_label} record {record_index}: "
                    f"length {sequence_length} != expected {expected_length}"
                )
