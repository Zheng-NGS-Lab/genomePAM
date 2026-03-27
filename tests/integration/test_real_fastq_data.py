"""Integration tests exercising pipeline functions against real downsampled FASTQ data.

All classes are marked @pytest.mark.integration. Tests use real 10K-read
Cas9 and FnCas12a FASTQ files checked into tests/data/fastq/ to validate
barcode distributions and structural properties.
"""

from collections import Counter
from pathlib import Path

import pytest


@pytest.mark.integration
class TestRealDataStructuralProperties:
    """Validate structural properties of the real test data itself."""

    @pytest.mark.parametrize(
        "subdirectory,r1_filename,r2_filename",
        [
            ("cas9", "SRR33421097_R1.fastq.gz", "SRR33421097_R2.fastq.gz"),
            ("fncas12a", "GNPAM327_R1.fastq.gz", "GNPAM327_R2.fastq.gz"),
        ],
    )
    def test_r1_and_r2_have_equal_read_counts(
        self,
        test_data_dir: Path,
        subdirectory: str,
        r1_filename: str,
        r2_filename: str,
    ) -> None:
        """R1 and R2 files contain the same number of records."""
        from umi.demultiplex import fq

        r1_count = sum(1 for _ in fq(str(test_data_dir / subdirectory / r1_filename)))
        r2_count = sum(1 for _ in fq(str(test_data_dir / subdirectory / r2_filename)))

        assert r1_count == r2_count, (
            f"{subdirectory}: R1 has {r1_count} records, R2 has {r2_count}"
        )

    @pytest.mark.parametrize(
        "subdirectory,r1_filename,r2_filename",
        [
            ("cas9", "SRR33421097_R1.fastq.gz", "SRR33421097_R2.fastq.gz"),
            ("fncas12a", "GNPAM327_R1.fastq.gz", "GNPAM327_R2.fastq.gz"),
        ],
    )
    def test_r1_and_r2_headers_are_paired(
        self,
        test_data_dir: Path,
        subdirectory: str,
        r1_filename: str,
        r2_filename: str,
    ) -> None:
        """First 200 R1/R2 pairs share the same read ID (before the space)."""
        from itertools import islice
        from umi.demultiplex import fq

        r1_path = str(test_data_dir / subdirectory / r1_filename)
        r2_path = str(test_data_dir / subdirectory / r2_filename)

        pair_count = 200
        r1_records = list(islice(fq(r1_path), pair_count))
        r2_records = list(islice(fq(r2_path), pair_count))

        for pair_index in range(pair_count):
            r1_read_id = r1_records[pair_index][0].split()[0]
            r2_read_id = r2_records[pair_index][0].split()[0]
            assert r1_read_id == r2_read_id, (
                f"Pair {pair_index}: R1 ID {r1_read_id} != R2 ID {r2_read_id}"
            )

    def test_cas9_dominant_barcode_frequency(self, test_data_dir: Path) -> None:
        """AGTGACAC appears in >50% of Cas9 reads at positions [10:18]."""
        from umi.demultiplex import fq

        cas9_dominant_barcode = "AGTGACAC"
        filepath = str(test_data_dir / "cas9" / "SRR33421097_R1.fastq.gz")
        barcode_counter = Counter()
        total_records = 0

        for record in fq(filepath):
            sequence = record[1].strip()
            barcode = sequence[10:18]
            barcode_counter[barcode] += 1
            total_records += 1

        dominant_fraction = barcode_counter[cas9_dominant_barcode] / total_records

        assert dominant_fraction > 0.50, (
            f"Expected {cas9_dominant_barcode} in >50% of reads, "
            f"got {dominant_fraction:.1%} ({barcode_counter[cas9_dominant_barcode]}/{total_records})"
        )

    def test_fncas12a_barcode_distribution(self, test_data_dir: Path) -> None:
        """The 4 expected FnCas12a barcodes account for >95% of reads."""
        from umi.demultiplex import fq

        fncas12a_expected_barcodes = {"CCATTGCC", "CCAATGCC", "CCAGTGCC", "CCACTGCC"}
        filepath = str(test_data_dir / "fncas12a" / "GNPAM327_R1.fastq.gz")
        barcode_counter = Counter()
        total_records = 0

        for record in fq(filepath):
            sequence = record[1].strip()
            barcode = sequence[10:18]
            barcode_counter[barcode] += 1
            total_records += 1

        expected_count = sum(
            barcode_counter[barcode] for barcode in fncas12a_expected_barcodes
        )
        combined_fraction = expected_count / total_records

        assert combined_fraction > 0.95, (
            f"Expected 4 barcodes to cover >95% of reads, "
            f"got {combined_fraction:.1%} ({expected_count}/{total_records})"
        )
