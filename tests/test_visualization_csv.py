"""Tests for CSV export in the visualization module."""

import csv
from pathlib import Path

import pytest

from modules.guideseq.visualization import exportOfftargetsCSV


class TestVisualizationCSV:
    """Tests for exportOfftargetsCSV function."""

    def test_export_csv_basic(self, tmp_path: Path) -> None:
        """Test basic CSV export with minimal data."""
        offtargets = [
            {
                'seq': 'ATGAGCCACTGTGCCTGGCC',
                'bulged_seq': '',
                'reads': 432,
                'target_seq': 'GTGAGCCACTGTGCCTGGCC',
                'realigned_target_seq': '',
            }
        ]
        target_seq = 'GTGAGCCACTGTGCCTGGCC'
        outfile = str(tmp_path / 'test_output')

        exportOfftargetsCSV(offtargets, target_seq, outfile)

        csv_path = Path(outfile + '.csv')
        assert csv_path.exists()

        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert rows[0]['sequence_type'] == 'no_bulge'
        assert rows[0]['offtarget_sequence'] == 'ATGAGCCACTGTGCCTGGCC'
        assert rows[0]['read_count'] == '432'
        assert rows[0]['target_sequence'] == 'GTGAGCCACTGTGCCTGGCC'
        assert rows[0]['realigned_target_sequence'] == ''
        assert rows[0]['genomic_range'] == ''
        assert rows[0]['region_type'] == ''
        assert rows[0]['gene'] == ''

    def test_export_csv_with_annotation(self, tmp_path: Path) -> None:
        """Test CSV includes annotation fields when provided."""
        offtargets = [
            {
                'seq': 'ATGAGCCACTGTGCCTGGCC',
                'bulged_seq': '',
                'reads': 432,
                'target_seq': 'GTGAGCCACTGTGCCTGGCC',
                'realigned_target_seq': '',
                'range': 'chr13:80238315-80238344',
                'type': 'Intergenic',
                'gene': 'LINC01080-SPRY2',
            }
        ]
        target_seq = 'GTGAGCCACTGTGCCTGGCC'
        outfile = str(tmp_path / 'test_output')

        exportOfftargetsCSV(offtargets, target_seq, outfile)

        csv_path = Path(outfile + '.csv')
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 1
        assert rows[0]['genomic_range'] == 'chr13:80238315-80238344'
        assert rows[0]['region_type'] == 'Intergenic'
        assert rows[0]['gene'] == 'LINC01080-SPRY2'

    def test_export_csv_both_bulge_and_nobulge(self, tmp_path: Path) -> None:
        """Test separate rows are created for bulged and non-bulged sequences."""
        offtargets = [
            {
                'seq': 'ATGAGCCACTGTGCCTGGCC',
                'bulged_seq': 'ATGAGCCAC-TGTGCCTGGCC',
                'reads': 432,
                'target_seq': 'GTGAGCCACTGTGCCTGGCC',
                'realigned_target_seq': 'GTGAGCCAC-TGTGCCTGGCC',
                'range': 'chr13:80238315-80238344',
                'type': 'Intergenic',
                'gene': 'LINC01080-SPRY2',
            }
        ]
        target_seq = 'GTGAGCCACTGTGCCTGGCC'
        outfile = str(tmp_path / 'test_output')

        exportOfftargetsCSV(offtargets, target_seq, outfile)

        csv_path = Path(outfile + '.csv')
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 2

        # First row should be no_bulge
        assert rows[0]['sequence_type'] == 'no_bulge'
        assert rows[0]['offtarget_sequence'] == 'ATGAGCCACTGTGCCTGGCC'
        assert rows[0]['read_count'] == '432'
        assert rows[0]['target_sequence'] == 'GTGAGCCACTGTGCCTGGCC'
        assert rows[0]['realigned_target_sequence'] == ''
        assert rows[0]['genomic_range'] == 'chr13:80238315-80238344'
        assert rows[0]['region_type'] == 'Intergenic'
        assert rows[0]['gene'] == 'LINC01080-SPRY2'

        # Second row should be bulge
        assert rows[1]['sequence_type'] == 'bulge'
        assert rows[1]['offtarget_sequence'] == 'ATGAGCCAC-TGTGCCTGGCC'
        assert rows[1]['read_count'] == '432'
        assert rows[1]['target_sequence'] == 'GTGAGCCACTGTGCCTGGCC'
        assert rows[1]['realigned_target_sequence'] == 'GTGAGCCAC-TGTGCCTGGCC'
        assert rows[1]['genomic_range'] == 'chr13:80238315-80238344'
        assert rows[1]['region_type'] == 'Intergenic'
        assert rows[1]['gene'] == 'LINC01080-SPRY2'
