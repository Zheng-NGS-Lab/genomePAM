#!/usr/bin/env python3
"""
Compare current trim_tag_umi.nf approach vs gold-standard approach.

Simulates both UMI extraction and barcode filtering in pure Python:
- Approach A: Current pipeline (4x BBDuk + umi-tools extract)
  - Trims R2 adapters in two passes: first Pass 1 (Read2Tail + polyG),
    then Pass 2 (RC-FIXSEQ with restrictright constraint)
- Approach B: Gold standard (cutadapt + umi-tools extract)
  - Trims R2 adapters in single pass (all adapters at once)

No external dependencies beyond Python 3 stdlib.

Note: R2 sequences may differ slightly (typically <1%) due to adapter
trimming order. When multiple adapters have overlapping regions, the order
of trimming can affect final sequence length. This is expected behavior.
"""

import argparse
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import NamedTuple


class Quality(NamedTuple):
    """Phred quality score."""
    score: int

    @staticmethod
    def from_char(char: str) -> "Quality":
        """Convert ASCII character to phred score (Illumina 1.8+)."""
        return Quality(ord(char) - 33)

    def is_low(self, threshold: int = 10) -> bool:
        """Check if score is below threshold."""
        return self.score < threshold


@dataclass
class FastqRead:
    """Single FASTQ record."""
    header: str
    sequence: str
    plus: str
    qualities: str

    @property
    def read_id(self) -> str:
        """Extract read ID (before space)."""
        return self.header.split()[0]

    @property
    def read_length(self) -> int:
        """Get sequence length."""
        return len(self.sequence)

    def extract_barcode(self, pos_start: int, pos_end: int) -> str:
        """Extract barcode at positions [pos_start, pos_end)."""
        return self.sequence[pos_start - 1:pos_end]

    def extract_umi(self, pos_start: int, pos_end: int) -> str:
        """Extract UMI at positions [pos_start, pos_end)."""
        return self.sequence[pos_start - 1:pos_end]

    def has_valid_barcode(self, barcode: str) -> bool:
        """Check if barcode contains only standard bases and is not poly-homopolymer."""
        if not barcode or len(barcode) != 8:
            return False
        if not all(base in "ACGT" for base in barcode):
            return False
        if len(set(barcode)) == 1:  # Poly-homopolymer (CCCC, AAAA, etc)
            return False
        return True

    def trim_3prime_quality(self, threshold: int = 10) -> "FastqRead":
        """Trim low-quality bases from 3' end (right)."""
        idx = len(self.sequence)
        for i in range(len(self.qualities) - 1, -1, -1):
            if Quality.from_char(self.qualities[i]).is_low(threshold):
                idx = i
            else:
                break
        return FastqRead(
            self.header,
            self.sequence[:idx],
            self.plus,
            self.qualities[:idx],
        )

    def trim_left(self, num_bases: int) -> "FastqRead":
        """Trim from left (5' end)."""
        if num_bases >= len(self.sequence):
            return FastqRead(self.header, "", self.plus, "")
        return FastqRead(
            self.header,
            self.sequence[num_bases:],
            self.plus,
            self.qualities[num_bases:],
        )


@dataclass
class AdapterTrimResult:
    """Result of adapter trimming."""
    read: FastqRead
    was_trimmed: bool


def read_fastq_gz(path: Path) -> list[FastqRead]:
    """Read all records from gzipped FASTQ file."""
    reads = []
    with gzip.open(str(path), "rt") as f:
        while True:
            header = f.readline().rstrip("\n")
            if not header:
                break
            sequence = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")
            qualities = f.readline().rstrip("\n")
            reads.append(FastqRead(header, sequence, plus, qualities))
    return reads


def find_best_overlap(sequence: str, adapter: str, min_overlap: int = 8) -> tuple[int, int]:
    """
    Find best suffix overlap of adapter in sequence.

    Returns:
        (match_start_in_seq, trimmed_length) or (len(sequence), 0) if no match
    """
    best_pos = len(sequence)
    best_length = 0

    for overlap in range(len(adapter), min_overlap - 1, -1):
        adapter_suffix = adapter[-overlap:]
        seq_suffix = sequence[-overlap:]

        mismatches = sum(1 for a, b in zip(adapter_suffix, seq_suffix) if a != b)
        max_mismatches = max(1, overlap // 8)

        if mismatches <= max_mismatches:
            best_pos = len(sequence) - overlap
            best_length = overlap
            break

    return best_pos, best_length


def trim_3prime_adapter(sequence: str, qualities: str, adapters: list[str]) -> tuple[str, str]:
    """
    Trim 3' adapters (multiple passes).

    Parameters:
        sequence: DNA sequence
        qualities: Quality string
        adapters: List of adapter sequences to search for

    Returns:
        (trimmed_sequence, trimmed_qualities)
    """
    for adapter in adapters:
        match_pos, match_len = find_best_overlap(sequence, adapter, min_overlap=8)
        if match_len > 0:
            sequence = sequence[:match_pos]
            qualities = qualities[:match_pos]

    return sequence, qualities


def reverse_complement(sequence: str) -> str:
    """Get reverse complement of DNA sequence."""
    complement = str.maketrans("ATCG", "TAGC")
    return sequence.translate(complement)[::-1]


def levenshtein_distance(s1: str, s2: str) -> int:
    """Compute edit distance between two strings."""
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (c1 != c2)
            current.append(min(insertions, deletions, substitutions))
        previous = current

    return previous[-1]


def detect_fixseq(reads_r1: list[FastqRead], pos_start: int, pos_end: int) -> str:
    """
    Auto-detect FIXSEQ from barcode frequencies.

    Extract barcodes from positions [pos_start, pos_end), filter by frequency
    and quality, return most common.
    """
    barcode_counts = Counter()

    for read in reads_r1[:10000]:  # Use first 10k reads
        barcode = read.extract_barcode(pos_start, pos_end)
        if read.has_valid_barcode(barcode):
            barcode_counts[barcode] += 1

    most_common = barcode_counts.most_common(1)
    return most_common[0][0] if most_common else ""


def build_whitelist(detected_fixseq: str, barcode_counts: Counter) -> list[str]:
    """
    Build barcode whitelist using fuzzy matching.

    Include barcodes within edit distance <= 2 of detected FIXSEQ.
    """
    if not detected_fixseq:
        return []

    whitelist = []
    for barcode, count in barcode_counts.most_common(100):
        if count < 10:
            break
        if levenshtein_distance(barcode, detected_fixseq) <= 2:
            whitelist.append(barcode)

    return whitelist


def approach_a(
    reads_r1: list[FastqRead],
    reads_r2: list[FastqRead],
    read1_tail: str,
    read2_tail: str,
    pos1: int,
    pos2: int,
) -> tuple[list[FastqRead], list[FastqRead], dict]:
    """
    Current pipeline approach (simulate trim_tag_umi.nf).

    Steps:
    1. Trim R1 3' adapters (Read1Tail + polyG)
    2. Auto-detect FIXSEQ
    3. Build whitelist via fuzzy matching
    4. Trim R2 3' adapters (Read2Tail + polyG)
    5. Trim R2 3' RC-FIXSEQ (second pass)
    6. Extract UMI/barcode, filter by whitelist
    """
    stats = {"steps": {}}

    # Step 1: Trim R1 adapters
    adapters_r1 = [read1_tail, "G" * 20]
    trimmed_r1 = []
    for read in reads_r1:
        seq, qual = trim_3prime_adapter(read.sequence, read.qualities, adapters_r1)
        trimmed_r1.append(FastqRead(read.header, seq, read.plus, qual))

    stats["steps"]["r1_after_trim"] = len(trimmed_r1)
    stats["steps"]["r1_median_length_after_trim"] = (
        median([r.read_length for r in trimmed_r1]) if trimmed_r1 else 0
    )

    # Step 2: Auto-detect FIXSEQ
    detected_fixseq = detect_fixseq(trimmed_r1, pos1, pos2)
    stats["detected_fixseq"] = detected_fixseq

    # Step 3: Build barcode whitelist
    barcode_counts = Counter()
    for read in trimmed_r1[:10000]:
        barcode = read.extract_barcode(pos1, pos2)
        if read.has_valid_barcode(barcode):
            barcode_counts[barcode] += 1

    whitelist = build_whitelist(detected_fixseq, barcode_counts)
    stats["whitelist"] = whitelist

    # Step 4: Trim R2 first pass (Read2Tail + polyG)
    adapters_r2_pass1 = [read2_tail, "G" * 20]
    trimmed_r2_pass1 = []
    for read in reads_r2:
        seq, qual = trim_3prime_adapter(read.sequence, read.qualities, adapters_r2_pass1)
        trimmed_r2_pass1.append(FastqRead(read.header, seq, read.plus, qual))

    stats["steps"]["r2_after_trim_pass1"] = len(trimmed_r2_pass1)
    stats["steps"]["r2_median_length_after_trim_pass1"] = (
        median([r.read_length for r in trimmed_r2_pass1]) if trimmed_r2_pass1 else 0
    )

    # Step 5: Trim R2 second pass (RC-FIXSEQ with restrictright constraint)
    rc_fixseq = reverse_complement(detected_fixseq)
    trimmed_r2_pass2 = []
    for read in trimmed_r2_pass1:
        # BBDuk restrictright means only search in first pos2 bases
        search_region = read.sequence[:pos2] if len(read.sequence) > pos2 else read.sequence
        match_pos, match_len = find_best_overlap(search_region, rc_fixseq, min_overlap=7)

        if match_len > 0:
            trimmed_seq = read.sequence[:match_pos]
            trimmed_qual = read.qualities[:match_pos]
        else:
            trimmed_seq = read.sequence
            trimmed_qual = read.qualities

        trimmed_r2_pass2.append(FastqRead(read.header, trimmed_seq, read.plus, trimmed_qual))

    stats["steps"]["r2_after_trim_pass2"] = len(trimmed_r2_pass2)
    stats["steps"]["r2_median_length_after_trim_pass2"] = (
        median([r.read_length for r in trimmed_r2_pass2]) if trimmed_r2_pass2 else 0
    )

    # Step 6: Extract UMI/barcode and filter by whitelist
    output_r1 = []
    output_r2 = []
    for r1, r2 in zip(trimmed_r1, trimmed_r2_pass2):
        umi = r1.extract_umi(1, 11)  # First 10 bases
        barcode = r1.extract_barcode(pos1, pos2)

        if barcode in whitelist:
            # Format header as umi-tools output: @ID_BARCODE_UMI
            new_header = f"{r1.read_id}_{barcode}_{umi}"

            # Trim first 18 bases from R1
            trimmed_r1_final = r1.trim_left(pos2)
            trimmed_r1_final.header = new_header

            trimmed_r2_final = r2
            trimmed_r2_final.header = new_header

            output_r1.append(trimmed_r1_final)
            output_r2.append(trimmed_r2_final)

    stats["steps"]["after_barcode_filter"] = len(output_r1)

    return output_r1, output_r2, stats


def approach_b(
    reads_r1: list[FastqRead],
    reads_r2: list[FastqRead],
    read1_tail: str,
    read2_tail: str,
    pos1: int,
    pos2: int,
    detected_fixseq: str,
    whitelist: list[str],
) -> tuple[list[FastqRead], list[FastqRead], dict]:
    """
    Gold-standard approach (single-pass R2 trimming).

    Same logic as Approach A but trims R2 adapters in one pass
    (Read2Tail + polyG + RC-FIXSEQ together).
    """
    stats = {"steps": {}}

    # Step 1: Trim R1 adapters (same as A)
    adapters_r1 = [read1_tail, "G" * 20]
    trimmed_r1 = []
    for read in reads_r1:
        seq, qual = trim_3prime_adapter(read.sequence, read.qualities, adapters_r1)
        trimmed_r1.append(FastqRead(read.header, seq, read.plus, qual))

    stats["steps"]["r1_after_trim"] = len(trimmed_r1)
    stats["steps"]["r1_median_length_after_trim"] = (
        median([r.read_length for r in trimmed_r1]) if trimmed_r1 else 0
    )

    # Step 2: Trim R2 in single pass (all adapters at once)
    rc_fixseq = reverse_complement(detected_fixseq)
    adapters_r2_all = [read2_tail, "G" * 20, rc_fixseq]
    trimmed_r2 = []
    for read in reads_r2:
        seq, qual = trim_3prime_adapter(read.sequence, read.qualities, adapters_r2_all)
        trimmed_r2.append(FastqRead(read.header, seq, read.plus, qual))

    stats["steps"]["r2_after_trim"] = len(trimmed_r2)
    stats["steps"]["r2_median_length_after_trim"] = (
        median([r.read_length for r in trimmed_r2]) if trimmed_r2 else 0
    )

    # Step 3: Extract UMI/barcode and filter by whitelist
    output_r1 = []
    output_r2 = []
    for r1, r2 in zip(trimmed_r1, trimmed_r2):
        umi = r1.extract_umi(1, 11)  # First 10 bases
        barcode = r1.extract_barcode(pos1, pos2)

        if barcode in whitelist:
            # Format header as umi-tools output: @ID_BARCODE_UMI
            new_header = f"{r1.read_id}_{barcode}_{umi}"

            # Trim first 18 bases from R1
            trimmed_r1_final = r1.trim_left(pos2)
            trimmed_r1_final.header = new_header

            trimmed_r2_final = r2
            trimmed_r2_final.header = new_header

            output_r1.append(trimmed_r1_final)
            output_r2.append(trimmed_r2_final)

    stats["steps"]["after_barcode_filter"] = len(output_r1)

    return output_r1, output_r2, stats


def compare_outputs(
    a_r1: list[FastqRead],
    a_r2: list[FastqRead],
    b_r1: list[FastqRead],
    b_r2: list[FastqRead],
) -> dict:
    """Compare outputs from both approaches."""
    comparison = {}

    # Read count
    comparison["read_count_a"] = len(a_r1)
    comparison["read_count_b"] = len(b_r1)
    comparison["read_count_diff"] = abs(len(a_r1) - len(b_r1))

    if len(a_r1) == 0 or len(b_r1) == 0:
        return comparison

    # Sequence comparison (for common reads)
    min_count = min(len(a_r1), len(b_r1))
    r1_identical = sum(
        1 for i in range(min_count) if a_r1[i].sequence == b_r1[i].sequence
    )
    r2_identical = sum(
        1 for i in range(min_count) if a_r2[i].sequence == b_r2[i].sequence
    )

    comparison["r1_sequences_identical"] = r1_identical
    comparison["r1_sequences_identical_pct"] = (
        100 * r1_identical / min_count if min_count > 0 else 0
    )
    comparison["r2_sequences_identical"] = r2_identical
    comparison["r2_sequences_identical_pct"] = (
        100 * r2_identical / min_count if min_count > 0 else 0
    )

    # Header comparison
    headers_identical = sum(
        1 for i in range(min_count) if a_r1[i].header == b_r1[i].header
    )
    comparison["headers_identical"] = headers_identical
    comparison["headers_identical_pct"] = 100 * headers_identical / min_count if min_count > 0 else 0

    # Length differences
    r1_length_diffs = defaultdict(int)
    r2_length_diffs = defaultdict(int)
    for i in range(min_count):
        r1_diff = len(a_r1[i].sequence) - len(b_r1[i].sequence)
        r2_diff = len(a_r2[i].sequence) - len(b_r2[i].sequence)
        r1_length_diffs[r1_diff] += 1
        r2_length_diffs[r2_diff] += 1

    comparison["r1_length_diffs"] = dict(sorted(r1_length_diffs.items()))
    comparison["r2_length_diffs"] = dict(sorted(r2_length_diffs.items()))

    return comparison


def print_results(dataset_name: str, stats_a: dict, stats_b: dict, comparison: dict):
    """Print comparison results."""
    print(f"\n=== Dataset: {dataset_name} ===")
    print(f"Input read pairs: {stats_a['steps'].get('r1_after_trim', 0)}")
    print()

    print(f"FIXSEQ detected: {stats_a['detected_fixseq']}")
    whitelist = stats_a["whitelist"]
    print(f"Whitelist barcodes ({len(whitelist)}): {whitelist[:10]}")
    if len(whitelist) > 10:
        print(f"  ... and {len(whitelist) - 10} more")
    print()

    print("Approach A (current pipeline - 4x BBDuk + umi-tools):")
    print(f"  After R1 adapter trim: {stats_a['steps'].get('r1_after_trim', 0)} reads, "
          f"median len {stats_a['steps'].get('r1_median_length_after_trim', 0)}")
    print(f"  After R2 adapter trim (pass 1): {stats_a['steps'].get('r2_after_trim_pass1', 0)} reads, "
          f"median len {stats_a['steps'].get('r2_median_length_after_trim_pass1', 0)}")
    print(f"  After R2 RC-FIXSEQ trim (pass 2): {stats_a['steps'].get('r2_after_trim_pass2', 0)} reads, "
          f"median len {stats_a['steps'].get('r2_median_length_after_trim_pass2', 0)}")
    print(f"  After barcode filter: {stats_a['steps'].get('after_barcode_filter', 0)} reads")
    print()

    print("Approach B (gold standard - cutadapt + umi-tools):")
    print(f"  After R1 adapter trim: {stats_b['steps'].get('r1_after_trim', 0)} reads, "
          f"median len {stats_b['steps'].get('r1_median_length_after_trim', 0)}")
    print(f"  After R2 adapter trim (single pass): {stats_b['steps'].get('r2_after_trim', 0)} reads, "
          f"median len {stats_b['steps'].get('r2_median_length_after_trim', 0)}")
    print(f"  After barcode filter: {stats_b['steps'].get('after_barcode_filter', 0)} reads")
    print()

    print("Comparison:")
    read_count_diff = comparison["read_count_diff"]
    read_count_a = comparison["read_count_a"]
    pct = 100 * read_count_diff / read_count_a if read_count_a > 0 else 0
    print(f"  Read count difference: {read_count_diff} ({pct:.1f}%)")
    print(f"  R1 sequences identical: {comparison['r1_sequences_identical']}/{comparison['read_count_a']} "
          f"({comparison['r1_sequences_identical_pct']:.1f}%)")
    print(f"  R2 sequences identical: {comparison['r2_sequences_identical']}/{comparison['read_count_a']} "
          f"({comparison['r2_sequences_identical_pct']:.1f}%)")
    print(f"  Headers identical: {comparison['headers_identical']}/{comparison['read_count_a']} "
          f"({comparison['headers_identical_pct']:.1f}%)")

    r1_diffs = comparison["r1_length_diffs"]
    if r1_diffs:
        print(f"  R1 length difference distribution: {r1_diffs}")
    r2_diffs = comparison["r2_length_diffs"]
    if r2_diffs:
        print(f"  R2 length difference distribution: {r2_diffs}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compare UMI trimming approaches"
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tests/data/fastq"),
        help="Base directory for test data",
    )
    args = parser.parse_args()

    # Pipeline parameters
    read1_tail = "AGATCGGAAGAGCACACGTC"
    read2_tail = "AGATCGGAAGAGCGTCGTGT"
    pos1 = 11
    pos2 = 18

    datasets = [
        ("cas9", "cas9/SRR33421097_R1.fastq.gz", "cas9/SRR33421097_R2.fastq.gz"),
        ("fncas12a", "fncas12a/GNPAM327_R1.fastq.gz", "fncas12a/GNPAM327_R2.fastq.gz"),
    ]

    for dataset_name, r1_path, r2_path in datasets:
        r1_file = args.data_dir / r1_path
        r2_file = args.data_dir / r2_path

        if not r1_file.exists() or not r2_file.exists():
            print(f"Skipping {dataset_name}: files not found")
            continue

        print(f"Loading {dataset_name} data...")
        reads_r1 = read_fastq_gz(r1_file)
        reads_r2 = read_fastq_gz(r2_file)

        # Run Approach A
        print(f"Running Approach A ({dataset_name})...")
        a_r1, a_r2, stats_a = approach_a(reads_r1, reads_r2, read1_tail, read2_tail, pos1, pos2)

        # Run Approach B using same detected FIXSEQ and whitelist for fair comparison
        print(f"Running Approach B ({dataset_name})...")
        b_r1, b_r2, stats_b = approach_b(
            reads_r1,
            reads_r2,
            read1_tail,
            read2_tail,
            pos1,
            pos2,
            stats_a["detected_fixseq"],
            stats_a["whitelist"],
        )

        # Compare outputs
        comparison = compare_outputs(a_r1, a_r2, b_r1, b_r2)

        # Print results
        print_results(dataset_name, stats_a, stats_b, comparison)


if __name__ == "__main__":
    main()
