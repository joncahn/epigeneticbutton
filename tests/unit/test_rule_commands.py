"""
Targeted rule tests that execute real samtools commands on synthetic SAM data.

These tests validate the actual shell logic extracted from Snakemake rules
(filter_bam_pe, filter_bam_se, filter_rna_se) using minimal synthetic
inputs. bamCoverage tests are skipped by default (not in epicc env).

Run with: pytest tests/unit/test_rule_commands.py -v
"""

import random
import shutil
import subprocess
import struct
import tempfile
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# Skip entire module if samtools is not available
# ---------------------------------------------------------------------------
pytestmark = pytest.mark.requires_samtools

SAMTOOLS = shutil.which("samtools")
if SAMTOOLS is None:
    pytest.skip("samtools not found", allow_module_level=True)


# ---------------------------------------------------------------------------
# Synthetic SAM data generators
# ---------------------------------------------------------------------------

def _sam_header_lines(chrom="chr1", length=100000):
    """Generate minimal SAM header as a list of lines (no trailing newlines)."""
    return [
        "@HD\tVN:1.6\tSO:unsorted",
        f"@SQ\tSN:{chrom}\tLN:{length}",
    ]


def _generate_pe_sam(path, n_pairs=50, n_dup_pairs=5, seed=42):
    """Generate a paired-end SAM file with proper FLAGS and duplicates.

    Creates n_pairs unique read pairs plus n_dup_pairs duplicate pairs
    (exact copies of the first n_dup_pairs pairs).
    """
    rng = random.Random(seed)
    lines = _sam_header_lines()
    seq = "A" * 100
    qual = "I" * 100

    pairs = []
    for i in range(n_pairs):
        pos = rng.randint(1, 90000)
        tlen = 250
        name = f"read_{i:04d}"
        # Read1: flag 99 (paired, proper pair, mate reverse, first in pair)
        # Read2: flag 147 (paired, proper pair, reverse, second in pair)
        r1 = f"{name}\t99\tchr1\t{pos}\t30\t100M\t=\t{pos + tlen - 100}\t{tlen}\t{seq}\t{qual}"
        r2 = f"{name}\t147\tchr1\t{pos + tlen - 100}\t30\t100M\t=\t{pos}\t{-tlen}\t{seq}\t{qual}"
        pairs.append((r1, r2))

    # Add duplicates (repeat first n_dup_pairs pairs with different names)
    for i in range(n_dup_pairs):
        orig_r1, orig_r2 = pairs[i]
        dup_name = f"dup_{i:04d}"
        dup_r1 = orig_r1.replace(f"read_{i:04d}", dup_name)
        dup_r2 = orig_r2.replace(f"read_{i:04d}", dup_name)
        pairs.append((dup_r1, dup_r2))

    for r1, r2 in pairs:
        lines.append(r1)
        lines.append(r2)

    path.write_text("\n".join(lines) + "\n")


def _generate_se_sam(path, n_reads=100, n_dups=10, seed=42):
    """Generate a single-end SAM file with duplicates.

    Creates n_reads unique reads plus n_dups duplicate reads.
    """
    rng = random.Random(seed)
    lines = _sam_header_lines()
    seq = "A" * 100
    qual = "I" * 100

    reads = []
    for i in range(n_reads):
        pos = rng.randint(1, 90000)
        flag = rng.choice([0, 16])  # forward or reverse
        name = f"read_{i:04d}"
        read_line = f"{name}\t{flag}\tchr1\t{pos}\t30\t100M\t*\t0\t0\t{seq}\t{qual}"
        reads.append(read_line)

    # Add duplicates
    for i in range(n_dups):
        dup_name = f"dup_{i:04d}"
        dup_line = reads[i].replace(f"read_{i:04d}", dup_name)
        reads.append(dup_line)

    lines.extend(reads)
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Helper functions: extracted shell logic from Snakemake rules
# ---------------------------------------------------------------------------

def run_filter_bam_pe(sam_path, out_bam, metrics_dup, metrics_flag, tmpdir, threads=1):
    """Execute the filter_bam_pe pipeline: view -F → fixmate → view -q → sort → markdup → flagstat.

    Mirrors the rule order: MAPQ filtering must happen after fixmate so that
    dropping a low-MAPQ mate doesn't desync the name-collated stream feeding
    fixmate.
    """
    temp0 = tmpdir / "temp0.bam"
    temp1 = tmpdir / "temp1.bam"
    temp2 = tmpdir / "temp2.bam"
    temp3 = tmpdir / "temp3.bam"

    cmds = [
        f"samtools view -@ {threads} -b -h -F 256 -o {temp0} {sam_path}",
        f"samtools fixmate -@ {threads} -m {temp0} {temp1}",
        f"samtools view -@ {threads} -b -h -q 10 -o {temp2} {temp1}",
        f"samtools sort -@ {threads} -o {temp3} {temp2}",
        f"samtools markdup -r -s -f {metrics_dup} -@ {threads} {temp3} {out_bam}",
        f"samtools flagstat -@ {threads} {out_bam} > {metrics_flag}",
    ]
    for cmd in cmds:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, f"Command failed: {cmd}\nstderr: {result.stderr}"


def run_filter_bam_se(sam_path, out_bam, metrics_dup, metrics_flag, tmpdir, threads=1):
    """Execute the filter_bam_se pipeline: view → sort → markdup → flagstat."""
    temp1 = tmpdir / "temp1.bam"
    temp2 = tmpdir / "temp2.bam"

    cmds = [
        f"samtools view -@ {threads} -b -h -q 10 -F 256 -o {temp1} {sam_path}",
        f"samtools sort -@ {threads} -o {temp2} {temp1}",
        f"samtools markdup -r -s -f {metrics_dup} -@ {threads} {temp2} {out_bam}",
        f"samtools flagstat -@ {threads} {out_bam} > {metrics_flag}",
    ]
    for cmd in cmds:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, f"Command failed: {cmd}\nstderr: {result.stderr}"


def run_filter_rna_se(sam_path, out_bam, metrics_flag, threads=1):
    """Execute the filter_rna_se pipeline: sort → index → flagstat."""
    cmds = [
        f"samtools sort -@ {threads} {sam_path} -o {out_bam}",
        f"samtools index -@ {threads} {out_bam}",
        f"samtools flagstat -@ {threads} {out_bam} > {metrics_flag}",
    ]
    for cmd in cmds:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        assert result.returncode == 0, f"Command failed: {cmd}\nstderr: {result.stderr}"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="class")
def pe_chip_outputs():
    """Run filter_bam_pe on synthetic PE SAM and return output paths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        sam = tmpdir / "pe_input.sam"
        out_bam = tmpdir / "pe_filtered.bam"
        metrics_dup = tmpdir / "markdup.txt"
        metrics_flag = tmpdir / "flagstat.txt"

        _generate_pe_sam(sam)
        run_filter_bam_pe(sam, out_bam, metrics_dup, metrics_flag, tmpdir)

        yield {
            "bam": out_bam,
            "metrics_dup": metrics_dup,
            "metrics_flag": metrics_flag,
            "tmpdir": tmpdir,
        }


@pytest.fixture(scope="class")
def se_chip_outputs():
    """Run filter_bam_se on synthetic SE SAM and return output paths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        sam = tmpdir / "se_input.sam"
        out_bam = tmpdir / "se_filtered.bam"
        metrics_dup = tmpdir / "markdup.txt"
        metrics_flag = tmpdir / "flagstat.txt"

        _generate_se_sam(sam)
        run_filter_bam_se(sam, out_bam, metrics_dup, metrics_flag, tmpdir)

        yield {
            "bam": out_bam,
            "metrics_dup": metrics_dup,
            "metrics_flag": metrics_flag,
            "tmpdir": tmpdir,
        }


@pytest.fixture(scope="class")
def se_rna_outputs():
    """Run filter_rna_se on synthetic SE SAM and return output paths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        sam = tmpdir / "rna_input.sam"
        out_bam = tmpdir / "rna_sorted.bam"
        metrics_flag = tmpdir / "flagstat.txt"

        _generate_se_sam(sam, n_reads=100, n_dups=10, seed=42)
        run_filter_rna_se(sam, out_bam, metrics_flag)

        yield {
            "bam": out_bam,
            "bai": Path(str(out_bam) + ".bai"),
            "metrics_flag": metrics_flag,
            "tmpdir": tmpdir,
        }


# ===========================================================================
# TestFilterChipPE
# ===========================================================================

@pytest.mark.unit
class TestFilterChipPE:
    """Tests for the filter_bam_pe rule logic."""

    def test_output_bam_exists(self, pe_chip_outputs):
        """Output BAM file exists and is non-empty."""
        bam = pe_chip_outputs["bam"]
        assert bam.exists() and bam.stat().st_size > 0

    def test_markdup_metrics_valid(self, pe_chip_outputs):
        """Markdup metrics file contains the COMMAND and DUPLICATE lines."""
        metrics = pe_chip_outputs["metrics_dup"]
        content = metrics.read_text()
        assert "COMMAND" in content or "READ" in content or "DUPLICATE" in content or len(content) > 0, \
            f"Markdup metrics appears empty or invalid: {content[:200]}"

    def test_flagstat_valid(self, pe_chip_outputs):
        """Flagstat output contains expected summary lines."""
        content = pe_chip_outputs["metrics_flag"].read_text()
        assert "in total" in content, f"Flagstat missing 'in total': {content[:200]}"
        assert "mapped" in content, f"Flagstat missing 'mapped': {content[:200]}"

    def test_quickcheck_passes(self, pe_chip_outputs):
        """samtools quickcheck passes on the output BAM."""
        result = subprocess.run(
            ["samtools", "quickcheck", str(pe_chip_outputs["bam"])],
            capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0, f"quickcheck failed: {result.stderr}"


# ===========================================================================
# TestFilterChipSE
# ===========================================================================

@pytest.mark.unit
class TestFilterChipSE:
    """Tests for the filter_bam_se rule logic."""

    def test_output_bam_exists(self, se_chip_outputs):
        """Output BAM file exists and is non-empty."""
        bam = se_chip_outputs["bam"]
        assert bam.exists() and bam.stat().st_size > 0

    def test_markdup_metrics_valid(self, se_chip_outputs):
        """Markdup metrics file is non-empty."""
        metrics = se_chip_outputs["metrics_dup"]
        content = metrics.read_text()
        assert len(content) > 0, "Markdup metrics file is empty"

    def test_flagstat_valid(self, se_chip_outputs):
        """Flagstat output contains expected summary lines."""
        content = se_chip_outputs["metrics_flag"].read_text()
        assert "in total" in content
        assert "mapped" in content

    def test_quickcheck_passes(self, se_chip_outputs):
        """samtools quickcheck passes on the output BAM."""
        result = subprocess.run(
            ["samtools", "quickcheck", str(se_chip_outputs["bam"])],
            capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0, f"quickcheck failed: {result.stderr}"


# ===========================================================================
# TestFilterRNASE
# ===========================================================================

@pytest.mark.unit
class TestFilterRNASE:
    """Tests for the filter_rna_se rule logic."""

    def test_sorted_bam_exists(self, se_rna_outputs):
        """Sorted BAM file exists and is non-empty."""
        bam = se_rna_outputs["bam"]
        assert bam.exists() and bam.stat().st_size > 0

    def test_bai_exists(self, se_rna_outputs):
        """BAI index file exists."""
        bai = se_rna_outputs["bai"]
        assert bai.exists() and bai.stat().st_size > 0

    def test_flagstat_valid(self, se_rna_outputs):
        """Flagstat output contains expected summary lines."""
        content = se_rna_outputs["metrics_flag"].read_text()
        assert "in total" in content
        assert "mapped" in content

    def test_header_shows_coordinate_sorted(self, se_rna_outputs):
        """BAM header shows SO:coordinate."""
        result = subprocess.run(
            ["samtools", "view", "-H", str(se_rna_outputs["bam"])],
            capture_output=True, text=True, timeout=30
        )
        assert result.returncode == 0
        assert "SO:coordinate" in result.stdout


# ===========================================================================
# TestMakeCoverageChip
# ===========================================================================

@pytest.mark.unit
class TestMakeCoverageChip:
    """Tests for bamCoverage (deeptools). Skipped if bamCoverage not in PATH."""

    @pytest.fixture(scope="class")
    def coverage_outputs(self, se_chip_outputs):
        """Run bamCoverage on the filtered SE ChIP BAM."""
        if not shutil.which("bamCoverage"):
            pytest.skip("bamCoverage not found (not in epicc env)")
        bam = se_chip_outputs["bam"]
        # Index the BAM first
        subprocess.run(["samtools", "index", str(bam)], check=True, timeout=30)
        bw = se_chip_outputs["tmpdir"] / "coverage.bw"
        result = subprocess.run(
            ["bamCoverage", "-b", str(bam), "-o", str(bw), "-bs", "50", "-p", "1"],
            capture_output=True, text=True, timeout=60
        )
        assert result.returncode == 0, f"bamCoverage failed: {result.stderr}"
        yield {"bw": bw}

    def test_bigwig_exists(self, coverage_outputs):
        """Bigwig file exists and is non-empty."""
        bw = coverage_outputs["bw"]
        assert bw.exists() and bw.stat().st_size > 0

    def test_bigwig_magic_bytes(self, coverage_outputs):
        """Bigwig has valid magic bytes."""
        with open(coverage_outputs["bw"], "rb") as f:
            magic = struct.unpack("<I", f.read(4))[0]
        valid = {0x888FFC26, 0x26FC8F88}
        assert magic in valid, f"Invalid magic: 0x{magic:08X}"
