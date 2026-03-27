"""Tests for FIXSEQ auto-detection logic in trim_tag_umi.nf.

The trim_tag_umi Nextflow process can auto-detect the FIXSEQ barcode
from the data or use an explicit override. This test verifies the
branching logic that resolves the effective FIXSEQ value.
"""

import subprocess
from pathlib import Path
from typing import Tuple

import pytest


class TestFixseqAutoDetection:
    """Verify FIXSEQ auto-detection logic resolves from data or uses explicit override."""

    FIXSEQ_SCRIPT = """
set -euo pipefail
# Write synthetic _barcodes file
cat > _barcodes <<'BARCODES'
{barcodes_content}
BARCODES

FIXSEQ="{fixseq_param}"
if [[ "$FIXSEQ" == "auto" ]]; then
    effective_fixseq=$(head -1 _barcodes)
    if [[ -z "$effective_fixseq" ]]; then
        echo "ERROR: FIXSEQ auto-detection failed -- no barcodes passed filtering at positions 1-2" >&2
        exit 1
    fi
    echo "AUTO:$effective_fixseq"
else
    effective_fixseq="$FIXSEQ"
    echo "EXPLICIT:$effective_fixseq"
fi
echo "$effective_fixseq" > _effective_fixseq
"""

    def _run_fixseq_logic(
        self, tmp_path: Path, barcodes_content: str, fixseq_param: str
    ) -> Tuple[int, str]:
        """Execute the FIXSEQ resolution logic via bash.

        Args:
            tmp_path: Pytest temp directory.
            barcodes_content: Content to write to _barcodes file.
            fixseq_param: Value of FIXSEQ parameter.

        Returns:
            Tuple of (return code, stdout).
        """
        script = self.FIXSEQ_SCRIPT.format(
            barcodes_content=barcodes_content,
            fixseq_param=fixseq_param,
        )
        result = subprocess.run(
            ["bash", "-c", script],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
        return result.returncode, result.stdout

    def test_auto_mode_uses_top_barcode(self, tmp_path: Path) -> None:
        """Auto mode detects FIXSEQ from the first barcode in the data."""
        barcodes = "CCAATGCC\nCCATTGCC\nCCAGTGCC\n"
        returncode, stdout = self._run_fixseq_logic(
            tmp_path, barcodes, "auto"
        )

        assert returncode == 0
        assert "AUTO:CCAATGCC" in stdout

        # Verify effective_fixseq was written to file
        effective_fixseq_file = tmp_path / "_effective_fixseq"
        assert effective_fixseq_file.exists()
        assert effective_fixseq_file.read_text().strip() == "CCAATGCC"

    def test_explicit_mode_uses_provided_value(self, tmp_path: Path) -> None:
        """Explicit mode uses the provided FIXSEQ value regardless of data."""
        barcodes = "CCAATGCC\nCCATTGCC\nCCAGTGCC\n"
        returncode, stdout = self._run_fixseq_logic(
            tmp_path, barcodes, "AGTGACAC"
        )

        assert returncode == 0
        assert "EXPLICIT:AGTGACAC" in stdout

        # Verify effective_fixseq was written to file
        effective_fixseq_file = tmp_path / "_effective_fixseq"
        assert effective_fixseq_file.exists()
        assert effective_fixseq_file.read_text().strip() == "AGTGACAC"

    def test_auto_mode_empty_barcodes_exits_with_error(self, tmp_path: Path) -> None:
        """Auto mode on empty _barcodes file exits with a diagnostic error."""
        script = self.FIXSEQ_SCRIPT.format(
            barcodes_content="",
            fixseq_param="auto",
        )
        result = subprocess.run(
            ["bash", "-c", script],
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "auto-detection failed" in result.stderr
