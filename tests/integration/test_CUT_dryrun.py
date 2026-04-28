"""
Integration tests for CUT&RUN / CUT&Tag workflows using Snakemake dry-run.

Verifies that the DAG resolves for all four CUT&x assay variants and
that the per-assay-family default peak callers (epic2 for *_broad,
SEACR for *_narrow) are reachable. Also exercises the CLI override path
that swaps the default for MACS2.

Test samples (from test_samples_CUT.tsv):
  CUT_RUN_broad  (H3K27me3 PE)  WT × 2 reps + IgG (shared Control)
  CUT_RUN_narrow (TF1 PE)       WT × 2 reps      (shared CR IgG)
  CUT_TAG_broad  (H3K27me3 SE)  WT × 2 reps + IgG (shared Control)
  CUT_TAG_narrow (TF1 SE)       WT × 2 reps      (shared CT IgG)
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun,
)


pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_CUT.yaml")


# ---------------------------------------------------------------------------
# Per-replicate peak targets (one per assay-variant × layout)
# ---------------------------------------------------------------------------
CR_BROAD_PE_PEAKS  = f"{_OUTPUT_DIR}/ChIP/peaks/peaks_pe__final__WT_CR_H3K27me3_rep1_peaks.broadPeak"
CR_NARROW_PE_PEAKS = f"{_OUTPUT_DIR}/ChIP/peaks/peaks_pe__final__WT_CR_TF1_rep1_peaks.narrowPeak"
CT_BROAD_SE_PEAKS  = f"{_OUTPUT_DIR}/ChIP/peaks/peaks_se__final__WT_CT_H3K27me3_rep1_peaks.broadPeak"
CT_NARROW_SE_PEAKS = f"{_OUTPUT_DIR}/ChIP/peaks/peaks_se__final__WT_CT_TF1_rep1_peaks.narrowPeak"

# Analysis-level (merged) peak targets — exercise replicate merging
CR_BROAD_ANALYSIS  = "CUT_RUN_broad__WT__H3K27me3__test_genome"
CR_NARROW_ANALYSIS = "CUT_RUN_narrow__WT__TF1__test_genome"
CT_BROAD_ANALYSIS  = "CUT_TAG_broad__WT__H3K27me3__test_genome"
CT_NARROW_ANALYSIS = "CUT_TAG_narrow__WT__TF1__test_genome"

CR_BROAD_SELECTED  = f"{_OUTPUT_DIR}/ChIP/peaks/selected_peaks__{CR_BROAD_ANALYSIS}.bedPeak"
CR_NARROW_SELECTED = f"{_OUTPUT_DIR}/ChIP/peaks/selected_peaks__{CR_NARROW_ANALYSIS}.bedPeak"
CT_BROAD_SELECTED  = f"{_OUTPUT_DIR}/ChIP/peaks/selected_peaks__{CT_BROAD_ANALYSIS}.bedPeak"
CT_NARROW_SELECTED = f"{_OUTPUT_DIR}/ChIP/peaks/selected_peaks__{CT_NARROW_ANALYSIS}.bedPeak"


@pytest.fixture(scope="module")
def test_options(repo_root):
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_CUT.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


# ===========================================================================
# Basic
# ===========================================================================

class TestBasic:

    def test_snakemake_installed(self, snakemake_available):
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_CUT.tsv"
        assert sample_file.exists()


# ===========================================================================
# CUT&RUN broad / narrow (PE)
# ===========================================================================

class TestCUTRUN:
    """CUT&RUN workflow dry-run (PE: paired-end fragments → bedgraph for SEACR;
    epic2 reads BAM directly for broad)."""

    def test_cr_broad_pe_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CR_BROAD_PE_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_cr_broad_uses_calling_peaks_pe(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, CR_BROAD_PE_PEAKS, ["--printshellcmds"]
        )
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "calling_peaks_pe" in output, \
            "CUT_RUN_broad PE should dispatch through calling_peaks_pe rule"

    def test_cr_broad_default_caller_is_epic2(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, CR_BROAD_PE_PEAKS, ["--printshellcmds"]
        )
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        # The dispatcher emits 'caller=epic2' in its banner via printf
        assert "caller=epic2" in output or "epic2" in output, \
            "CUT_RUN_broad default caller should be epic2"

    def test_cr_narrow_pe_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CR_NARROW_PE_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_cr_narrow_default_caller_is_seacr(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, CR_NARROW_PE_PEAKS, ["--printshellcmds"]
        )
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "SEACR" in output or "seacr" in output, \
            "CUT_RUN_narrow default caller should be SEACR"


# ===========================================================================
# CUT&Tag broad / narrow (SE)
# ===========================================================================

class TestCUTTag:
    """CUT&Tag SE: SEACR uses direct genomecov; epic2 reads BAM directly."""

    def test_ct_broad_se_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CT_BROAD_SE_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_ct_broad_uses_calling_peaks_se(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, CT_BROAD_SE_PEAKS, ["--printshellcmds"]
        )
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "calling_peaks_se" in output, \
            "CUT_TAG_broad SE should dispatch through calling_peaks_se rule"

    def test_ct_narrow_se_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CT_NARROW_SE_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# Replicate merging — CUT&x with 2 reps must reach pseudoreps + selected_peaks
# ===========================================================================

class TestReplicateMerging:

    def test_cr_broad_selected_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CR_BROAD_SELECTED)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_cr_narrow_selected_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CR_NARROW_SELECTED)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_ct_broad_selected_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CT_BROAD_SELECTED)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_ct_narrow_selected_peaks_dryrun(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, CT_NARROW_SELECTED)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# Caller override — switch to MACS2 via --config
# ===========================================================================

class TestCallerOverride:
    """A user can swap the default caller per family via cut_callpeaks.{broad,narrow}_caller.

    These tests pass --config overrides to verify the override is honored
    without rewriting the options file."""

    def test_cr_broad_can_use_macs2(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, CR_BROAD_PE_PEAKS,
            ["--config", "cut_callpeaks={'broad_caller':'macs2','narrow_caller':'seacr'}",
             "--printshellcmds"],
        )
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "caller=macs2" in output or "macs2 callpeak" in output, \
            "CUT_RUN_broad with broad_caller=macs2 override should use MACS2"

    def test_cr_narrow_can_use_macs2(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(
            repo_root, test_options, CR_NARROW_PE_PEAKS,
            ["--config", "cut_callpeaks={'broad_caller':'epic2','narrow_caller':'macs2'}",
             "--printshellcmds"],
        )
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "caller=macs2" in output or "macs2 callpeak" in output, \
            "CUT_RUN_narrow with narrow_caller=macs2 override should use MACS2"
