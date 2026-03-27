"""Tests for the umi demultiplex module."""

from pathlib import Path
from typing import Dict

import pytest

from .conftest import write_fastq_record, write_fastq_file, make_index_sequence


# -- Fixtures ----------------------------------------------------------------

BARCODE_SAMPLE_A = "ACGTACG"  # 7-base barcode half for sample A
BARCODE_SAMPLE_B = "TGCATGC"  # 7-base barcode half for sample B


@pytest.fixture
def sample_barcodes() -> Dict[str, str]:
    """Barcode-to-sample-name mapping for two samples."""
    combined = BARCODE_SAMPLE_A + BARCODE_SAMPLE_A
    return {combined: "sampleA"}


@pytest.fixture
def gzipped_fastq_quartet(tmp_path: Path) -> Dict[str, Path]:
    """Four gzipped FASTQ files (read1, read2, index1, index2) with 3 reads.

    All reads carry the sampleA barcode (ACGTACG + ACGTACG).

    Returns:
        Dict with keys 'read1', 'read2', 'index1', 'index2' mapping to file paths.
    """
    read_seq = "ATCGATCGATCGATCG"
    qual = "I" * len(read_seq)
    index_seq = make_index_sequence(BARCODE_SAMPLE_A)
    index_qual = "I" * len(index_seq)

    records_r1, records_r2, records_i1, records_i2 = [], [], [], []
    for i in range(3):
        records_r1.append(write_fastq_record(f"read{i}/1", read_seq, qual))
        records_r2.append(write_fastq_record(f"read{i}/2", read_seq, qual))
        records_i1.append(write_fastq_record(f"read{i}/i1", index_seq, index_qual))
        records_i2.append(write_fastq_record(f"read{i}/i2", index_seq, index_qual))

    return {
        "read1": write_fastq_file(tmp_path / "reads.r1.fastq.gz", records_r1, gzipped=True),
        "read2": write_fastq_file(tmp_path / "reads.r2.fastq.gz", records_r2, gzipped=True),
        "index1": write_fastq_file(tmp_path / "reads.i1.fastq.gz", records_i1, gzipped=True),
        "index2": write_fastq_file(tmp_path / "reads.i2.fastq.gz", records_i2, gzipped=True),
    }


# -- Tests -------------------------------------------------------------------


class TestFastqReader:
    """Tests for demultiplex.fq() gzipped FASTQ parsing."""

    def test_reads_correct_number_of_records(
        self, gzipped_fastq_quartet: Dict[str, Path]
    ) -> None:
        """fq() yields one 4-element list per FASTQ record."""
        from umi.demultiplex import fq

        records = list(fq(str(gzipped_fastq_quartet["read1"])))
        assert len(records) == 3

    def test_record_structure_is_four_lines(
        self, gzipped_fastq_quartet: Dict[str, Path]
    ) -> None:
        """Each yielded record is a list of 4 newline-terminated strings."""
        from umi.demultiplex import fq

        record = next(fq(str(gzipped_fastq_quartet["read1"])))
        assert len(record) == 4
        assert record[0].startswith("@")
        assert record[2].startswith("+")
        for line in record:
            assert line.endswith("\n")


class TestGetSampleId:
    """Tests for barcode-based sample identification."""

    def test_known_barcode_returns_sample_name(self) -> None:
        """When the barcode matches, return the sample name."""
        from umi.demultiplex import get_sample_id

        index_seq = make_index_sequence(BARCODE_SAMPLE_A) + "\n"
        i1 = ["@idx\n", index_seq, "+\n", "IIII\n"]
        i2 = ["@idx\n", index_seq, "+\n", "IIII\n"]
        sample_names = {BARCODE_SAMPLE_A + BARCODE_SAMPLE_A: "sampleA"}

        result = get_sample_id(i1, i2, sample_names)
        assert result == "sampleA"

    def test_unknown_barcode_returns_raw_barcode(self) -> None:
        """When the barcode is not in the lookup, return the raw barcode string."""
        from umi.demultiplex import get_sample_id

        index_seq = make_index_sequence(BARCODE_SAMPLE_B) + "\n"
        i1 = ["@idx\n", index_seq, "+\n", "IIII\n"]
        i2 = ["@idx\n", index_seq, "+\n", "IIII\n"]

        result = get_sample_id(i1, i2, {})
        assert result == BARCODE_SAMPLE_B + BARCODE_SAMPLE_B


class TestDemultiplex:
    """End-to-end tests for the demultiplex function."""

    def test_reads_below_min_reads_go_to_undetermined(
        self,
        gzipped_fastq_quartet: Dict[str, Path],
        sample_barcodes: Dict[str, str],
        tmp_path: Path,
    ) -> None:
        """With 3 reads and min_reads=5, all reads land in undetermined."""
        from umi.demultiplex import demultiplex

        out_dir = tmp_path / "demux_out"
        demultiplex(
            str(gzipped_fastq_quartet["read1"]),
            str(gzipped_fastq_quartet["read2"]),
            str(gzipped_fastq_quartet["index1"]),
            str(gzipped_fastq_quartet["index2"]),
            sample_barcodes,
            str(out_dir),
            min_reads=5,
        )

        undetermined = out_dir / "undetermined.r1.fastq"
        assert undetermined.exists()
        assert undetermined.read_text().count("@") == 3

        sample_file = out_dir / "sampleA.r1.fastq"
        assert not sample_file.exists()

    def test_reads_at_min_reads_produce_sample_fastq(
        self,
        gzipped_fastq_quartet: Dict[str, Path],
        sample_barcodes: Dict[str, str],
        tmp_path: Path,
    ) -> None:
        """With 3 reads and min_reads=3, sample FASTQ is created."""
        from umi.demultiplex import demultiplex

        out_dir = tmp_path / "demux_out"
        demultiplex(
            str(gzipped_fastq_quartet["read1"]),
            str(gzipped_fastq_quartet["read2"]),
            str(gzipped_fastq_quartet["index1"]),
            str(gzipped_fastq_quartet["index2"]),
            sample_barcodes,
            str(out_dir),
            min_reads=3,
        )

        sample_r1 = out_dir / "sampleA.r1.fastq"
        sample_r2 = out_dir / "sampleA.r2.fastq"
        assert sample_r1.exists()
        assert sample_r2.exists()
