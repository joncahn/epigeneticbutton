"""
Unit tests for workflow/scripts/validate_dmc_input.py

Tests the validation functions for dmC (direct methylation) inputs:
- validate_modbam: validates modBAM files for MM/ML tags
- validate_bedmethyl: validates bedMethyl file format
"""

import pytest
import gzip
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

# Import the validation functions
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'workflow', 'scripts'))
from validate_dmc_input import validate_modbam, validate_bedmethyl


class TestValidateBedMethyl:
    """Tests for validate_bedmethyl function."""

    def test_valid_bedmethyl_11_columns(self, temp_dir, valid_bedmethyl_content):
        """Test validation passes for valid 11-column bedMethyl."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(valid_bedmethyl_content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is True
        assert "Valid bedMethyl format" in message
        assert "4 lines validated" in message
        assert "3 chromosomes" in message

    def test_valid_bedmethyl_10_columns(self, temp_dir, valid_bedmethyl_10col_content):
        """Test validation passes for valid 10-column bedMethyl (no percent)."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(valid_bedmethyl_10col_content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is True
        assert "Valid bedMethyl format" in message

    def test_bedmethyl_with_header_and_comments(self, temp_dir, bedmethyl_with_header):
        """Test that headers and comments are properly skipped."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(bedmethyl_with_header)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is True
        assert "2 lines validated" in message  # Only data lines counted

    def test_bedmethyl_gzipped(self, temp_dir, valid_bedmethyl_content):
        """Test validation works with gzipped bedMethyl files."""
        bedmethyl_path = temp_dir / "test.bedmethyl.gz"
        with gzip.open(bedmethyl_path, 'wt') as f:
            f.write(valid_bedmethyl_content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is True
        assert "Valid bedMethyl format" in message

    def test_bedmethyl_file_not_found(self):
        """Test validation fails for non-existent file."""
        is_valid, message = validate_bedmethyl("/nonexistent/file.bedmethyl")

        assert is_valid is False
        assert "File not found" in message

    def test_bedmethyl_too_few_columns(self, temp_dir):
        """Test validation fails when fewer than 10 columns."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text("Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\n")

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Expected at least 10 columns" in message

    def test_bedmethyl_invalid_coordinates_non_numeric(self, temp_dir):
        """Test validation fails with non-numeric coordinates."""
        content = "Chr1\tabc\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Invalid coordinates" in message

    def test_bedmethyl_negative_start(self, temp_dir):
        """Test validation fails with negative start position."""
        content = "Chr1\t-10\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Start position cannot be negative" in message

    def test_bedmethyl_end_less_than_start(self, temp_dir):
        """Test validation fails when end <= start."""
        content = "Chr1\t200\t100\tm\t1000\t+\t100\t101\t255,0,0\t10\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "End (100) must be greater than start (200)" in message

    def test_bedmethyl_invalid_strand(self, temp_dir):
        """Test validation fails with invalid strand character."""
        content = "Chr1\t100\t101\tm\t1000\tX\t100\t101\t255,0,0\t10\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Invalid strand" in message

    def test_bedmethyl_negative_coverage(self, temp_dir):
        """Test validation fails with negative coverage."""
        content = "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t-5\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Coverage cannot be negative" in message

    def test_bedmethyl_invalid_coverage_value(self, temp_dir):
        """Test validation fails with non-numeric coverage."""
        content = "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\tabc\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Invalid coverage value" in message

    def test_bedmethyl_percent_out_of_range(self, temp_dir):
        """Test validation fails when percent modified > 100."""
        content = "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\t150.0\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Percent modified should be 0-100" in message

    def test_bedmethyl_invalid_percent_value(self, temp_dir):
        """Test validation fails with non-numeric percent modified."""
        content = "Chr1\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\tabc\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "Invalid percent modified value" in message

    def test_bedmethyl_empty_file(self, temp_dir):
        """Test validation fails for empty file."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text("")

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "No valid data lines found" in message

    def test_bedmethyl_only_headers(self, temp_dir):
        """Test validation fails when only headers/comments present."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text("track name=methylation\n# Comment line\n")

        is_valid, message = validate_bedmethyl(str(bedmethyl_path))

        assert is_valid is False
        assert "No valid data lines found" in message

    def test_bedmethyl_with_chrom_sizes(self, temp_dir, valid_bedmethyl_content, sample_chrom_sizes):
        """Test validation with chromosome size validation."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(valid_bedmethyl_content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path), sample_chrom_sizes)

        assert is_valid is True
        assert "Valid bedMethyl format" in message

    def test_bedmethyl_position_exceeds_chrom_length(self, temp_dir, sample_chrom_sizes):
        """Test validation fails when position exceeds chromosome length."""
        # Chr1 length is 30427671 in sample_chrom_sizes
        content = "Chr1\t100\t50000000\tm\t1000\t+\t100\t50000000\t255,0,0\t10\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path), sample_chrom_sizes)

        assert is_valid is False
        assert "exceeds chromosome Chr1 length" in message

    def test_bedmethyl_no_matching_chromosomes(self, temp_dir, sample_chrom_sizes):
        """Test validation fails when no chromosomes match reference."""
        content = "chrZ\t100\t101\tm\t1000\t+\t100\t101\t255,0,0\t10\t75.5\n"
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(content)

        is_valid, message = validate_bedmethyl(str(bedmethyl_path), sample_chrom_sizes)

        assert is_valid is False
        assert "No matching chromosomes" in message


class TestValidateModBAM:
    """Tests for validate_modbam function."""

    @patch('subprocess.run')
    def test_valid_modbam_with_mm_ml_tags(self, mock_run, temp_dir, mock_bam_header, mock_bam_read_with_mm_ml):
        """Test validation passes for modBAM with MM and ML tags."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        # Mock samtools view output
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout=mock_bam_read_with_mm_ml * 10,  # 10 reads
            stderr=""
        )

        is_valid, message = validate_modbam(str(bam_path))

        assert is_valid is True
        assert "Valid modBAM with MM/ML tags" in message
        assert mock_run.called

    @patch('subprocess.run')
    def test_modbam_without_mm_tag(self, mock_run, temp_dir, mock_bam_read_without_mm):
        """Test validation fails when MM tag is missing."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        mock_run.return_value = MagicMock(
            returncode=0,
            stdout=mock_bam_read_without_mm * 10,
            stderr=""
        )

        is_valid, message = validate_modbam(str(bam_path))

        assert is_valid is False
        assert "No MM (methylation) tags found" in message

    @patch('subprocess.run')
    def test_modbam_without_ml_tag(self, mock_run, temp_dir):
        """Test validation fails when ML tag is missing."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        # Read with MM but no ML tag
        read_with_mm_only = (
            "read001\t0\tChr1\t1000\t60\t100M\t*\t0\t0\t"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\t"
            "*" * 42 + "\t"
            "MM:Z:C+m,10,5,8;\n"
        )

        mock_run.return_value = MagicMock(
            returncode=0,
            stdout=read_with_mm_only * 10,
            stderr=""
        )

        is_valid, message = validate_modbam(str(bam_path))

        assert is_valid is False
        assert "No ML (methylation likelihood) tags found" in message

    def test_modbam_file_not_found(self):
        """Test validation fails for non-existent file."""
        is_valid, message = validate_modbam("/nonexistent/file.bam")

        assert is_valid is False
        assert "File not found" in message

    @patch('subprocess.run')
    def test_modbam_empty_file(self, mock_run, temp_dir):
        """Test validation fails for BAM with no reads."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="",
            stderr=""
        )

        is_valid, message = validate_modbam(str(bam_path))

        assert is_valid is False
        assert "BAM file contains no reads" in message

    @patch('subprocess.run')
    def test_modbam_samtools_error(self, mock_run, temp_dir):
        """Test validation fails when samtools returns error."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        mock_run.return_value = MagicMock(
            returncode=1,
            stdout="",
            stderr="Error: truncated file"
        )

        is_valid, message = validate_modbam(str(bam_path))

        assert is_valid is False
        assert "Failed to read BAM file" in message

    @patch('subprocess.run')
    def test_modbam_timeout(self, mock_run, temp_dir):
        """Test validation handles timeout gracefully."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        mock_run.side_effect = subprocess.TimeoutExpired(cmd="samtools", timeout=120)

        is_valid, message = validate_modbam(str(bam_path))

        assert is_valid is False
        assert "Timeout while reading BAM file" in message

    @patch('subprocess.run')
    def test_modbam_with_chrom_sizes_matching(self, mock_run, temp_dir,
                                               mock_bam_header, mock_bam_read_with_mm_ml,
                                               sample_chrom_sizes):
        """Test validation with chromosome size checking - matching chromosomes."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        # First call: samtools view (for reads)
        # Second call: samtools view -H (for header)
        mock_run.side_effect = [
            MagicMock(returncode=0, stdout=mock_bam_read_with_mm_ml * 10, stderr=""),
            MagicMock(returncode=0, stdout=mock_bam_header, stderr="")
        ]

        is_valid, message = validate_modbam(str(bam_path), sample_chrom_sizes)

        assert is_valid is True
        assert "Valid modBAM with MM/ML tags" in message

    @patch('subprocess.run')
    def test_modbam_no_matching_chromosomes(self, mock_run, temp_dir,
                                            mock_bam_read_with_mm_ml,
                                            sample_chrom_sizes):
        """Test validation fails when BAM chromosomes don't match reference."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        # Header with completely different chromosomes
        mismatched_header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chrZ\tLN:1000000\n"
            "@SQ\tSN:chrY\tLN:2000000\n"
        )

        mock_run.side_effect = [
            MagicMock(returncode=0, stdout=mock_bam_read_with_mm_ml * 10, stderr=""),
            MagicMock(returncode=0, stdout=mismatched_header, stderr="")
        ]

        is_valid, message = validate_modbam(str(bam_path), sample_chrom_sizes)

        assert is_valid is False
        assert "No matching chromosomes between BAM and reference" in message

    @patch('subprocess.run')
    def test_modbam_low_chromosome_coverage(self, mock_run, temp_dir,
                                            mock_bam_read_with_mm_ml,
                                            sample_chrom_sizes):
        """Test validation fails when < 50% of reference chromosomes are in BAM."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        # Header with only 1 out of 5 chromosomes matching (20% < 50%)
        low_coverage_header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:Chr1\tLN:30427671\n"
            "@SQ\tSN:chrX\tLN:1000000\n"
            "@SQ\tSN:chrY\tLN:2000000\n"
        )

        mock_run.side_effect = [
            MagicMock(returncode=0, stdout=mock_bam_read_with_mm_ml * 10, stderr=""),
            MagicMock(returncode=0, stdout=low_coverage_header, stderr="")
        ]

        is_valid, message = validate_modbam(str(bam_path), sample_chrom_sizes)

        assert is_valid is False
        assert "of reference chromosomes found in BAM" in message

    @patch('subprocess.run')
    def test_modbam_header_read_error(self, mock_run, temp_dir, mock_bam_read_with_mm_ml):
        """Test validation handles header read errors."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        mock_run.side_effect = [
            MagicMock(returncode=0, stdout=mock_bam_read_with_mm_ml * 10, stderr=""),
            MagicMock(returncode=1, stdout="", stderr="Error reading header")
        ]

        chrom_sizes = str(temp_dir / "chrom.sizes")
        Path(chrom_sizes).write_text("Chr1\t1000000\n")

        is_valid, message = validate_modbam(str(bam_path), chrom_sizes)

        assert is_valid is False
        assert "Failed to read BAM header" in message


class TestValidationMain:
    """Tests for the main() function and CLI behavior."""

    @patch('subprocess.run')
    def test_main_modbam_success(self, mock_run, temp_dir, mock_bam_read_with_mm_ml):
        """Test main function with modBAM input type."""
        bam_path = temp_dir / "test.bam"
        bam_path.touch()

        mock_run.return_value = MagicMock(
            returncode=0,
            stdout=mock_bam_read_with_mm_ml * 10,
            stderr=""
        )

        # Import and test main
        import validate_dmc_input
        sys.argv = ['validate_dmc_input.py', 'modBAM', str(bam_path)]

        with pytest.raises(SystemExit) as exc_info:
            validate_dmc_input.main()

        assert exc_info.value.code == 0

    def test_main_bedmethyl_success(self, temp_dir, valid_bedmethyl_content):
        """Test main function with bedMethyl input type."""
        bedmethyl_path = temp_dir / "test.bedmethyl"
        bedmethyl_path.write_text(valid_bedmethyl_content)

        import validate_dmc_input
        sys.argv = ['validate_dmc_input.py', 'bedMethyl', str(bedmethyl_path)]

        with pytest.raises(SystemExit) as exc_info:
            validate_dmc_input.main()

        assert exc_info.value.code == 0

    def test_main_unknown_input_type(self):
        """Test main function with unknown input type."""
        import validate_dmc_input
        sys.argv = ['validate_dmc_input.py', 'unknown', '/tmp/test.file']

        with pytest.raises(SystemExit) as exc_info:
            validate_dmc_input.main()

        assert exc_info.value.code == 1

    def test_main_insufficient_arguments(self):
        """Test main function with too few arguments."""
        import validate_dmc_input
        sys.argv = ['validate_dmc_input.py', 'modBAM']

        with pytest.raises(SystemExit) as exc_info:
            validate_dmc_input.main()

        assert exc_info.value.code == 1
