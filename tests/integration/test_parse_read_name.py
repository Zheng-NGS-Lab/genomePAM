"""Integration tests for the parseReadName function."""

import pytest


@pytest.mark.integration
class TestParseReadName:
    """Tests for parseReadName() function that extracts UMI from read names.

    The function parses umi-tools deduplicated read names in the format:
    READID_UMI, where UMI is a concatenated string of valid DNA bases [ACGTN].
    """

    def test_standard_umi_tools_format(self) -> None:
        """Parse standard umi-tools format with long UMI sequence.

        After umi-tools deduplication, read names have the format READID_UMI
        where the UMI is appended as a single concatenated token. This test
        validates that the function correctly extracts the UMI segment.
        """
        from identifyOfftargetSites import parseReadName

        read_name = "READ1_AACCGGTTACGTACGTTGCATGCA"
        umi, count = parseReadName(read_name, umi_len=24, idx_len=8)

        assert umi == "AACCGGTTACGTACGTTGCATGCA"
        assert count == 1

    def test_read_name_without_umi_returns_none(self) -> None:
        """Read name without underscore-prefixed bases returns None values.

        If a read name does not contain an underscore followed by valid DNA
        bases at the end, the function should return (None, None).
        """
        from identifyOfftargetSites import parseReadName

        read_name = "SRR12345.67890"
        umi, count = parseReadName(read_name, umi_len=0, idx_len=0)

        assert umi is None
        assert count is None

    def test_umi_with_only_acgtn_characters(self) -> None:
        """Only valid DNA bases [ACGTN] are matched as UMI.

        The regex pattern _([ACGTN]+)$ ensures that only sequences containing
        the valid DNA bases are recognized as UMI. This test confirms that
        a UMI with the complete set of valid bases is extracted.
        """
        from identifyOfftargetSites import parseReadName

        read_name = "MYREAD_ACGTNACGT"
        umi, count = parseReadName(read_name, umi_len=8, idx_len=4)

        assert umi == "ACGTNACGT"
        assert count == 1

    def test_numeric_suffix_not_mistaken_for_umi(self) -> None:
        """Numeric suffixes without underscore are not parsed as UMI.

        Read names with dot-separated numeric identifiers (e.g., typical
        SRA format) should not be parsed as containing UMI since the regex
        requires an underscore before the bases.
        """
        from identifyOfftargetSites import parseReadName

        read_name = "SRR123.456"
        umi, count = parseReadName(read_name, umi_len=0, idx_len=0)

        assert umi is None
        assert count is None

    def test_underscore_with_non_base_characters_returns_none(self) -> None:
        """Non-base characters after underscore prevent UMI extraction.

        The regex pattern _([ACGTN]+)$ requires that all characters after
        the underscore be valid DNA bases. If any invalid characters (digits,
        lowercase letters, special characters) are present, no match occurs.
        """
        from identifyOfftargetSites import parseReadName

        read_name = "READ1_ABC123"
        umi, count = parseReadName(read_name, umi_len=0, idx_len=0)

        assert umi is None
        assert count is None

    def test_multiple_underscores_extracts_last_segment(self) -> None:
        """When multiple underscores present, the last [ACGTN]+ segment is matched.

        The regex pattern with $ anchor matches from the rightmost position.
        This ensures that read names with multiple underscore-separated fields
        will match only the final valid DNA base sequence.
        """
        from identifyOfftargetSites import parseReadName

        read_name = "READ_001_ACGTACGT"
        umi, count = parseReadName(read_name, umi_len=8, idx_len=0)

        assert umi == "ACGTACGT"
        assert count == 1

    def test_count_always_one_after_deduplication(self) -> None:
        """The count field is always 1 for deduplicated reads.

        After umi-tools deduplication, each unique molecular ID appears
        exactly once. The count return value is always 1, independent of
        the umi_len and idx_len parameters (which are kept for compatibility).
        """
        from identifyOfftargetSites import parseReadName

        read_names = [
            "READ1_ACGTACGT",
            "READ2_AAAAAAAAAA",
            "READ3_TTTTGGGGCCCCAAAA",
        ]

        for read_name in read_names:
            umi, count = parseReadName(read_name, umi_len=10, idx_len=5)
            assert count == 1, f"Expected count=1 for {read_name}, got {count}"
