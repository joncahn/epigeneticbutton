"""
Pytest configuration and shared fixtures for EpigeneticButton tests.
"""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_chrom_sizes(temp_dir):
    """Create a minimal chrom.sizes file for testing."""
    chrom_sizes_path = temp_dir / "test.chrom.sizes"
    chrom_sizes_path.write_text(
        "Chr1\t30427671\n"
        "Chr2\t19698289\n"
        "Chr3\t23459830\n"
        "Chr4\t18585056\n"
        "Chr5\t26975502\n"
    )
    return str(chrom_sizes_path)


@pytest.fixture
def valid_bedmethyl_content():
    """Return valid bedMethyl format content (11 columns)."""
    return (
        "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\t75.5\n"
        "Chr1\t200\t201\tm\t800\t-\t200\t201\t255,0,0\t8\t62.3\n"
        "Chr2\t150\t151\tm\t950\t+\t150\t151\t255,0,0\t12\t85.0\n"
        "Chr3\t300\t301\tm\t700\t+\t300\t301\t255,0,0\t7\t50.0\n"
    )


@pytest.fixture
def valid_bedmethyl_10col_content():
    """Return valid bedMethyl format with 10 columns (no percent modified)."""
    return (
        "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\n"
        "Chr1\t200\t201\tm\t800\t-\t200\t201\t255,0,0\t8\n"
        "Chr2\t150\t151\tm\t950\t+\t150\t151\t255,0,0\t12\n"
    )


@pytest.fixture
def bedmethyl_with_header():
    """Return bedMethyl content with track header."""
    return (
        "track name=methylation\n"
        "# This is a comment\n"
        "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\t75.5\n"
        "Chr1\t200\t201\tm\t800\t-\t200\t201\t255,0,0\t8\t62.3\n"
    )


@pytest.fixture
def mock_bam_header():
    """Return a mock BAM header with chromosome information."""
    return (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:Chr1\tLN:30427671\n"
        "@SQ\tSN:Chr2\tLN:19698289\n"
        "@SQ\tSN:Chr3\tLN:23459830\n"
        "@PG\tID:minimap2\tPN:minimap2\tVN:2.24\n"
    )


@pytest.fixture
def mock_bam_read_with_mm_ml():
    """Return a mock BAM read line with MM and ML tags."""
    return (
        "read001\t0\tChr1\t1000\t60\t100M\t*\t0\t0\t"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
        "*" * 100 + "\t"
        "NM:i:0\tms:i:200\tAS:i:200\tnn:i:0\t"
        "MM:Z:C+m,10,5,8;\t"
        "ML:B:C,255,128,200,100,150\n"
    )


@pytest.fixture
def mock_bam_read_without_mm():
    """Return a mock BAM read line without MM tag."""
    return (
        "read001\t0\tChr1\t1000\t60\t100M\t*\t0\t0\t"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
        "*" * 100 + "\t"
        "NM:i:0\tms:i:200\tAS:i:200\tnn:i:0\n"
    )


@pytest.fixture
def dmc_sample_names():
    """Return sample name examples for direct methylation (dmC)."""
    return {
        "dmc_modbam": "mC__Col0__leaf__dmC__rep1__ColCEN",
        "dmc_bedmethyl": "mC__Col0__leaf__bedMethyl__rep1__ColCEN",
        "bismark_wgbs": "mC__Col0__leaf__WGBS__rep1__ColCEN",
        "bismark_pico": "mC__Col0__leaf__Pico__rep1__ColCEN",
        "bismark_emseq": "mC__Col0__leaf__EMseq__rep1__ColCEN",
        "bismark_default": "mC__Col0__leaf__sample__rep1__ColCEN",
    }
