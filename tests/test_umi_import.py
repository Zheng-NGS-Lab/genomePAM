"""Tests for umi package import chain and signature compatibility."""

import inspect

import pytest


class TestUmiImportChain:
    """Verify the umi package structure and importability."""

    def test_import_demultiplex_module(self) -> None:
        """The demultiplex module imports without error."""
        from umi import demultiplex

        assert hasattr(demultiplex, "demultiplex")
        assert hasattr(demultiplex, "fq")
        assert hasattr(demultiplex, "get_sample_id")

    def test_consolidate_import_skipped_without_htseq(self) -> None:
        """Consolidate requires HTSeq; verify the import fails gracefully outside Docker."""
        try:
            from umi import consolidate
            # If HTSeq is installed, just confirm the module loaded
            assert hasattr(consolidate, "consolidate")
        except ImportError:
            pytest.skip("HTSeq not available (expected outside Docker)")


class TestSignatureCompatibility:
    """Verify function signatures match the call sites in guideseq.py."""

    def test_demultiplex_accepts_seven_arguments(self) -> None:
        """guideseq.py calls demultiplex.demultiplex() with 7 positional args."""
        from umi.demultiplex import demultiplex

        sig = inspect.signature(demultiplex)
        params = list(sig.parameters.keys())
        assert params == [
            "read1", "read2", "index1", "index2",
            "sample_barcodes", "out_dir", "min_reads",
        ]

    def test_consolidate_accepts_four_arguments(self) -> None:
        """guideseq.py calls consolidate.consolidate() with 4 positional args."""
        try:
            from umi.consolidate import consolidate
        except ImportError:
            pytest.skip("HTSeq not available")

        sig = inspect.signature(consolidate)
        params = list(sig.parameters.keys())
        assert params == [
            "fastq_file", "consolidated_fastq_file", "min_qual", "min_freq",
        ]
