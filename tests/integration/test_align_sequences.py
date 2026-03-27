"""Integration tests for alignSequences() and regexFromSequence() functions."""

import pytest


@pytest.mark.integration
class TestAlignSequences:
    """Tests for alignSequences() function that performs fuzzy sequence alignment.

    This function aligns a targetsite sequence (guide RNA target) against a window
    sequence (genomic region), returning match information including position, strand,
    mismatches, and indel information. Results are a 15-element list structured as:
    [offtarget_sequence_no_bulge, mismatches, strand_m, start_no_bulge, end_no_bulge,
     bulged_offtarget_sequence, length, score, substitutions, insertions, deletions,
     strand_b, bulged_start, bulged_end, realigned_target]
    """

    def test_perfect_match_returns_zero_mismatches(self) -> None:
        """A perfect match in the window should return 0 mismatches.

        When the target sequence appears exactly in the window (with flanking
        padding), the function should identify it with 0 mismatches on the
        forward strand.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        # Embed target in window with 30bp of flanking A's on each side
        window = "A" * 30 + target + "A" * 30

        result = alignSequences(target, window, max_score=7)

        # result[0]: offtarget_sequence_no_bulge should match target
        assert result[0] == target, f"Expected {target}, got {result[0]}"
        # result[1]: mismatches should be 0
        assert result[1] == 0, f"Expected 0 mismatches, got {result[1]}"
        # result[2]: strand should be '+'
        assert result[2] == "+", f"Expected '+' strand, got {result[2]}"

    def test_single_mismatch_detected(self) -> None:
        """A single SNP should be detected and counted as 1 mismatch.

        When the target sequence has one base substitution in the window,
        the function should detect this with a mismatch count of 1.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        # Create window with target having first G -> A substitution
        modified_target = "ATGAGCCACTGTGCCTGGCC"
        window = "A" * 30 + modified_target + "A" * 30

        result = alignSequences(target, window, max_score=7)

        # result[0]: should match the modified sequence in the window
        assert result[0] == modified_target
        # result[1]: mismatches should be 1
        assert result[1] == 1, f"Expected 1 mismatch, got {result[1]}"

    def test_reverse_complement_match(self) -> None:
        """The reverse complement of the target should match on the minus strand.

        When the window contains the reverse complement of the target sequence,
        the alignment should be detected on the reverse strand ('-').
        """
        from identifyOfftargetSites import alignSequences, reverseComplement

        target = "GTGAGCCACTGTGCCTGGCC"
        # Reverse complement of target
        rc_target = reverseComplement(target)  # Should be "GGCCAGGCACAGTGGCTCAC"
        window = "A" * 30 + rc_target + "A" * 30

        result = alignSequences(target, window, max_score=7)

        # result[2]: strand should be '-' for reverse complement match
        assert result[2] == "-", f"Expected '-' strand for RC match, got {result[2]}"
        # result[0]: should contain the reverse complement
        assert rc_target in (result[0], reverseComplement(result[0]))

    def test_no_match_within_score_threshold(self) -> None:
        """Completely unrelated sequence should return empty result fields.

        When the window sequence is entirely unrelated to the target,
        with no matches within the max_score threshold, the function
        should return empty strings for non-match fields.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        # Window with unrelated sequence (all T's)
        window = "T" * 100

        result = alignSequences(target, window, max_score=7)

        # result[0]: offtarget_sequence_no_bulge should be empty
        assert result[0] == "", f"Expected empty string for no match, got {result[0]}"
        # result[1]: mismatches should be empty
        assert result[1] == "", f"Expected empty string for mismatches, got {result[1]}"
        # result[2]: strand should be empty
        assert result[2] == "", f"Expected empty strand, got {result[2]}"

    def test_indel_detection(self) -> None:
        """Insertions or deletions should be detected in bulged alignment.

        When the window contains an insertion relative to the target,
        the function should detect this via the bulged alignment fields.
        The target sequence should include 'N' as a PAM marker for proper
        indel handling.
        """
        from identifyOfftargetSites import alignSequences

        # Target with N as PAM marker (required for proper indel realignment)
        target = "GTGAGCCACTGTGCCTGGCCN"
        # Window has insertion: extra 'A' in the middle
        window_with_insertion = "A" * 30 + "GTGAGCCAACTGTGCCTGGCCN" + "A" * 30

        result = alignSequences(target, window_with_insertion, max_score=7)

        # Either bulged_offtarget_sequence is non-empty (result[5])
        # or insertions/deletions are detected (result[9] or result[10] are non-zero)
        has_bulge = result[5] != ""
        # Handle both int and empty string return types
        insertions = result[9] if isinstance(result[9], int) else 0
        deletions = result[10] if isinstance(result[10], int) else 0
        has_indels = insertions > 0 or deletions > 0

        assert has_bulge or has_indels, (
            f"Expected indel detection. bulged_seq='{result[5]}', "
            f"insertions={result[9]}, deletions={result[10]}"
        )

    def test_result_is_fifteen_elements(self) -> None:
        """alignSequences() must always return exactly 15 elements.

        The return value is always a 15-element list, regardless of whether
        a match is found or not.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        window_perfect = "A" * 30 + target + "A" * 30
        window_no_match = "T" * 100

        result_perfect = alignSequences(target, window_perfect, max_score=7)
        result_no_match = alignSequences(target, window_no_match, max_score=7)

        assert len(result_perfect) == 15, (
            f"Expected 15 elements for perfect match, got {len(result_perfect)}"
        )
        assert len(result_no_match) == 15, (
            f"Expected 15 elements for no match, got {len(result_no_match)}"
        )

    def test_multiple_mismatches_at_threshold(self) -> None:
        """Multiple mismatches up to max_score should be detected.

        With max_score=7, the function should find and report alignments
        with up to 7 substitutions.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        # Create sequence with 2 substitutions (G->A at pos 0, C->A at pos 8)
        window_2mm = "A" * 30 + "ATGAGCCAATGTGCCTGGCC" + "A" * 30

        result = alignSequences(target, window_2mm, max_score=7)

        # Should find match with 2 mismatches
        assert result[0] != "", "Should find match with 2 mismatches"
        assert result[1] == 2, f"Expected 2 mismatches, got {result[1]}"

    def test_custom_max_score_threshold(self) -> None:
        """max_score parameter controls the mismatch tolerance threshold.

        When max_score is set to a low value (e.g., 1), sequences with more
        than 1 mismatch should not be detected as matches.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        # Sequence with 2 mismatches (G->A at pos 0, C->A at pos 8)
        window_2mm = "A" * 30 + "ATGAGCCAATGTGCCTGGCC" + "A" * 30

        # With max_score=1, this should not match (has 2 mismatches)
        result = alignSequences(target, window_2mm, max_score=1)

        # Should return empty match
        assert result[0] == "", (
            f"With max_score=1, sequence with 2 mismatches should not match. "
            f"Got {result[0]} with {result[1]} mismatches"
        )

    def test_lowercase_input_is_handled(self) -> None:
        """Lowercase sequences should be converted to uppercase internally.

        The function calls window_sequence.upper(), so lowercase input
        should work correctly and return uppercase result.
        """
        from identifyOfftargetSites import alignSequences

        target = "GTGAGCCACTGTGCCTGGCC"
        window_lowercase = "a" * 30 + "gtgagccactgtgcctggcc" + "a" * 30

        result = alignSequences(target, window_lowercase, max_score=7)

        # Should find perfect match with uppercase result
        assert result[0] == target, (
            f"Expected uppercase result {target}, got {result[0]}"
        )
        assert result[1] == 0, f"Expected 0 mismatches for perfect match"


@pytest.mark.integration
class TestRegexFromSequence:
    """Tests for regexFromSequence() that generates fuzzy regex patterns.

    This function converts a DNA sequence (potentially with IUPAC ambiguous
    base codes) into two regex patterns for fuzzy matching: one for
    substitutions only, and one allowing indels.
    """

    def test_unambiguous_sequence(self) -> None:
        """Unambiguous bases should expand to single-base character classes.

        For sequence "ACGT", the pattern should contain [A], [C], [G], [T]
        for each position (or their IUPAC equivalents, which for unambiguous
        bases are just the base itself).
        """
        from identifyOfftargetSites import regexFromSequence

        pattern_std, pattern_gap = regexFromSequence("ACGT")

        # Check that pattern_std contains character classes for each base
        # The pattern will be (?b:[A][C][G][T]){s<=7}
        assert "[A]" in pattern_std, f"Expected [A] in pattern, got {pattern_std}"
        assert "[C]" in pattern_std, f"Expected [C] in pattern, got {pattern_std}"
        assert "[G]" in pattern_std, f"Expected [G] in pattern, got {pattern_std}"
        assert "[T]" in pattern_std, f"Expected [T] in pattern, got {pattern_std}"

    def test_ambiguous_bases_expanded(self) -> None:
        """Ambiguous IUPAC bases should be expanded to all possible bases.

        'N' should expand to [GATC], and 'Y' should expand to [CT]
        according to IUPAC codes.
        """
        from identifyOfftargetSites import regexFromSequence

        # Test N (expands to GATC)
        pattern_std_n, _ = regexFromSequence("N")
        # The pattern should contain [GATC] or similar (order may vary)
        assert "[" in pattern_std_n and "]" in pattern_std_n, (
            f"Expected character class in pattern, got {pattern_std_n}"
        )
        # Verify all bases are present (in any order)
        for base in ["G", "A", "T", "C"]:
            assert base in pattern_std_n, (
                f"Expected base {base} in N expansion, got {pattern_std_n}"
            )

        # Test Y (expands to CT)
        pattern_std_y, _ = regexFromSequence("Y")
        assert "C" in pattern_std_y and "T" in pattern_std_y, (
            f"Expected C and T in Y expansion, got {pattern_std_y}"
        )

    def test_returns_two_patterns(self) -> None:
        """regexFromSequence() must return a tuple of exactly 2 pattern strings.

        The function returns (pattern_standard, pattern_gap) for two different
        fuzzy matching modes.
        """
        from identifyOfftargetSites import regexFromSequence

        result = regexFromSequence("ACGT")

        assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
        assert len(result) == 2, f"Expected 2 elements, got {len(result)}"
        assert isinstance(result[0], str), f"Expected string for pattern_standard"
        assert isinstance(result[1], str), f"Expected string for pattern_gap"

    def test_error_threshold_in_pattern(self) -> None:
        """Custom error threshold should appear in both patterns.

        When errors=5 is passed, both patterns should contain {s<=5} and
        include the indel specifications in pattern_gap.
        """
        from identifyOfftargetSites import regexFromSequence

        pattern_std, pattern_gap = regexFromSequence("ACGT", errors=5)

        # Standard pattern should have {s<=5}
        assert "{s<=5}" in pattern_std, (
            f"Expected {{s<=5}} in standard pattern, got {pattern_std}"
        )

        # Gap pattern should have appropriate indel/error specs
        assert "{" in pattern_gap and "}" in pattern_gap, (
            f"Expected error specification in gap pattern, got {pattern_gap}"
        )
        # Should mention indel constraints and error threshold
        assert "i<=" in pattern_gap or "d<=" in pattern_gap, (
            f"Expected indel constraints in gap pattern, got {pattern_gap}"
        )

    def test_lookahead_atomic_grouping(self) -> None:
        """By default, lookahead=True should wrap pattern in atomic grouping.

        With lookahead=True (default), patterns should start with (?b:
        (possessive atomic grouping in regex module).
        """
        from identifyOfftargetSites import regexFromSequence

        pattern_std_lookahead, _ = regexFromSequence("ACGT", lookahead=True)
        pattern_std_no_lookahead, _ = regexFromSequence("ACGT", lookahead=False)

        assert pattern_std_lookahead.startswith("(?b:"), (
            f"Expected (?b: prefix with lookahead=True, got {pattern_std_lookahead}"
        )
        assert not pattern_std_no_lookahead.startswith("(?b:"), (
            f"Expected no (?b: prefix with lookahead=False, "
            f"got {pattern_std_no_lookahead}"
        )

    def test_indel_parameter_affects_gap_pattern(self) -> None:
        """indels parameter should control allowed indel count in gap pattern.

        With indels=2, the gap pattern should allow i<=2 and d<=2.
        With indels=0, no indels should be allowed.
        """
        from identifyOfftargetSites import regexFromSequence

        _, pattern_gap_2 = regexFromSequence("ACGT", indels=2, errors=7)
        _, pattern_gap_0 = regexFromSequence("ACGT", indels=0, errors=7)

        # Pattern with indels=2 should mention i<=2 and d<=2
        assert "i<=2" in pattern_gap_2, (
            f"Expected i<=2 in gap pattern with indels=2, got {pattern_gap_2}"
        )
        assert "d<=2" in pattern_gap_2, (
            f"Expected d<=2 in gap pattern with indels=2, got {pattern_gap_2}"
        )

        # Pattern with indels=0 should mention i<=0 and d<=0
        assert "i<=0" in pattern_gap_0, (
            f"Expected i<=0 in gap pattern with indels=0, got {pattern_gap_0}"
        )
        assert "d<=0" in pattern_gap_0, (
            f"Expected d<=0 in gap pattern with indels=0, got {pattern_gap_0}"
        )

    def test_mixed_ambiguous_and_unambiguous_bases(self) -> None:
        """A sequence with both ambiguous and unambiguous bases should work.

        Sequence like "ACGTN" should have [A], [C], [G], [T], and [GATC]
        (or similar) in the pattern.
        """
        from identifyOfftargetSites import regexFromSequence

        pattern_std, _ = regexFromSequence("ACGTN")

        # Should contain the four unambiguous bases
        assert "[A]" in pattern_std
        assert "[C]" in pattern_std
        assert "[G]" in pattern_std
        assert "[T]" in pattern_std

        # Should also have the ambiguous N base (all four bases)
        assert "[" in pattern_std and "]" in pattern_std, "Should have character classes"

    def test_pattern_strings_are_valid_format(self) -> None:
        """Returned patterns should be in valid regex format.

        Patterns should be syntactically valid and contain expected
        regex components like character classes and quantifiers.
        """
        from identifyOfftargetSites import regexFromSequence

        pattern_std, pattern_gap = regexFromSequence("ACGT", errors=7, indels=1)

        # Both patterns should contain character classes
        assert pattern_std.count("[") == pattern_std.count("]"), (
            f"Mismatched brackets in pattern_std: {pattern_std}"
        )
        assert pattern_gap.count("[") == pattern_gap.count("]"), (
            f"Mismatched brackets in pattern_gap: {pattern_gap}"
        )

        # Both should contain fuzzy match syntax
        assert "{" in pattern_std and "}" in pattern_std, (
            f"Missing fuzzy syntax in pattern_std: {pattern_std}"
        )
        assert "{" in pattern_gap and "}" in pattern_gap, (
            f"Missing fuzzy syntax in pattern_gap: {pattern_gap}"
        )

    def test_single_base_sequence(self) -> None:
        """A single-base sequence should produce valid patterns.

        Even a sequence of length 1 should produce valid patterns with
        appropriate error thresholds.
        """
        from identifyOfftargetSites import regexFromSequence

        pattern_std, pattern_gap = regexFromSequence("A", errors=3)

        # Should still have character class and error specification
        assert "[" in pattern_std and "]" in pattern_std
        assert "{s<=3}" in pattern_std
        assert isinstance(pattern_gap, str) and len(pattern_gap) > 0
