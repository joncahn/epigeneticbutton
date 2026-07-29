"""
Integration tests for ATAC-seq using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for ATAC-seq samples without actually executing the rules.

Test samples (from test_samples_ATAC.tsv):
  WT_ATAC_rep1        ATAC  PE  genotype:WT
  WT_ATAC_rep2        ATAC  PE  genotype:WT
  mutant_ATAC_rep1    ATAC  PE  genotype:mutant
  mutant_ATAC_rep2    ATAC  PE  genotype:mutant
  mutant2_ATAC_rep1   ATAC  PE  genotype:mutant2
  mutant2_ATAC_rep2   ATAC  PE  genotype:mutant2
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_ATAC.yaml")


# ---------------------------------------------------------------------------
# Target path constants
# ---------------------------------------------------------------------------

# Per-replicate coverage bigwig targets
WT_REP1_BIGWIG = f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__WT_ATAC_rep1.bw"
WT_REP2_BIGWIG = f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__WT_ATAC_rep2.bw"
MUTANT_REP1_BIGWIG = f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__mutant_ATAC_rep1.bw"

# Per-replicate peak targets
WT_REP1_PEAKS = f"{_OUTPUT_DIR}/ATAC/peaks/peaks_atac__final__WT_ATAC_rep1_peaks.narrowPeak"
MUTANT_REP1_PEAKS = f"{_OUTPUT_DIR}/ATAC/peaks/peaks_atac__final__mutant_ATAC_rep1_peaks.narrowPeak"

# Analysis-level names (Assay__levels_label__Genome)
ATAC_WT_ANALYSIS = "ATAC__WT__test_genome"
ATAC_MUTANT_ANALYSIS = "ATAC__mutant__test_genome"
ATAC_MUTANT2_ANALYSIS = "ATAC__mutant2__test_genome"

# Merged coverage bigwig targets (analysis-level, 2+ reps)
WT_MERGED_BIGWIG = f"{_OUTPUT_DIR}/ATAC/tracks/coverage__merged__{ATAC_WT_ANALYSIS}.bw"
MUTANT_MERGED_BIGWIG = f"{_OUTPUT_DIR}/ATAC/tracks/coverage__merged__{ATAC_MUTANT_ANALYSIS}.bw"

# Selected peaks (analysis-level, after pseudo-replicate filtering)
WT_SELECTED_PEAKS = f"{_OUTPUT_DIR}/ATAC/peaks/selected_peaks__{ATAC_WT_ANALYSIS}.bedPeak"

# Differential peaks (MAnorm pairwise comparison)
DIFFPEAKS_WT_VS_MUTANT = f"{_OUTPUT_DIR}/ATAC/peaks/{ATAC_WT_ANALYSIS}_vs_{ATAC_MUTANT_ANALYSIS}/{ATAC_WT_ANALYSIS}_vs_{ATAC_MUTANT_ANALYSIS}_all_MAvalues.xls"


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_ATAC.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


# ===========================================================================
# TestBasic
# ===========================================================================

class TestBasic:
    """Basic sanity checks for the ATAC test configuration."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_ATAC.tsv"
        assert sample_file.exists()


# ===========================================================================
# TestATACWorkflow
# ===========================================================================

class TestATACWorkflow:
    """Test ATAC-seq workflow dry-run for per-replicate targets."""

    def test_atac_pe_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for ATAC PE sample bigwig."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_BIGWIG)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_atac_includes_shift_bam(self, snakemake_available, repo_root, test_options):
        """Test that ATAC workflow includes atac_shift_bam rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "atac_shift_bam" in output, "ATAC workflow should use atac_shift_bam rule"

    def test_atac_includes_bam_to_bed(self, snakemake_available, repo_root, test_options):
        """Test that ATAC peak calling workflow includes atac_bam_to_bed rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_PEAKS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "atac_bam_to_bed" in output, "ATAC peak calling should use atac_bam_to_bed rule"

    def test_atac_includes_make_coverage(self, snakemake_available, repo_root, test_options):
        """Test that ATAC bigwig generation includes make_coverage_atac rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "make_coverage_atac" in output, "ATAC bigwig should use make_coverage_atac rule"


# ===========================================================================
# TestATACPeakCalling
# ===========================================================================

class TestATACPeakCalling:
    """Test ATAC-seq peak calling workflow."""

    def test_peak_target_resolves(self, snakemake_available, repo_root, test_options):
        """Test that ATAC peak target resolves in dry-run."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_peak_calling_uses_macs2(self, snakemake_available, repo_root, test_options):
        """Test that ATAC peak calling uses calling_peaks_atac rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_REP1_PEAKS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "calling_peaks_atac" in output, "ATAC peaks should use calling_peaks_atac rule"

    def test_selected_peaks_resolve(self, snakemake_available, repo_root, test_options):
        """Test that analysis-level selected peaks resolve."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_SELECTED_PEAKS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# TestReplicateMerging
# ===========================================================================

class TestReplicateMerging:
    """Test replicate merging rules for ATAC samples with multiple replicates."""

    def test_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """WT has 2 reps -- merged replicate bigwig should trigger merging."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_MERGED_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_bam_replicates" in output, \
            "ATAC with 2 reps should trigger merging_bam_replicates rule"

    def test_merged_bigwig_triggers_shift(self, snakemake_available, repo_root, test_options):
        """Merged bigwig should trigger atac_shift_bam on the merged BAM."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WT_MERGED_BIGWIG, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "atac_shift_bam" in output, \
            "Merged ATAC bigwig should trigger atac_shift_bam on merged BAM"

    def test_mutant_merged_bigwig_resolves(self, snakemake_available, repo_root, test_options):
        """Mutant has 2 reps -- merged bigwig should resolve."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, MUTANT_MERGED_BIGWIG)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# TestDifferentialPeaks
# ===========================================================================

class TestDifferentialPeaks:
    """Test MAnorm differential peak analysis between ATAC conditions."""

    def test_diffpeaks_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that MAnorm WT vs mutant dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DIFFPEAKS_WT_VS_MUTANT)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_diffpeaks_uses_manorm(self, snakemake_available, repo_root, test_options):
        """Test that differential peaks use perform_pairwise_diff_peaks rule."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DIFFPEAKS_WT_VS_MUTANT, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "perform_pairwise_diff_peaks" in output, \
            "Differential peaks should use perform_pairwise_diff_peaks rule"


# ===========================================================================
# TestWildcardResolution
# ===========================================================================

class TestWildcardResolution:
    """Test wildcard resolution for ATAC-seq targets."""

    def test_all_replicate_bigwigs_resolve(self, snakemake_available, repo_root, test_options):
        """Test that bigwigs resolve for all per-replicate samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__WT_ATAC_rep1.bw",
            f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__WT_ATAC_rep2.bw",
            f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__mutant_ATAC_rep1.bw",
            f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__mutant_ATAC_rep2.bw",
            f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__mutant2_ATAC_rep1.bw",
            f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__mutant2_ATAC_rep2.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_all_replicate_peaks_resolve(self, snakemake_available, repo_root, test_options):
        """Test that peaks resolve for all per-replicate samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/ATAC/peaks/peaks_atac__final__WT_ATAC_rep1_peaks.narrowPeak",
            f"{_OUTPUT_DIR}/ATAC/peaks/peaks_atac__final__mutant_ATAC_rep1_peaks.narrowPeak",
            f"{_OUTPUT_DIR}/ATAC/peaks/peaks_atac__final__mutant2_ATAC_rep1_peaks.narrowPeak",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_analysis_level_targets_resolve(self, snakemake_available, repo_root, test_options):
        """Test that analysis-level selected peaks resolve for all genotypes."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/ATAC/peaks/selected_peaks__{ATAC_WT_ANALYSIS}.bedPeak",
            f"{_OUTPUT_DIR}/ATAC/peaks/selected_peaks__{ATAC_MUTANT_ANALYSIS}.bedPeak",
            f"{_OUTPUT_DIR}/ATAC/peaks/selected_peaks__{ATAC_MUTANT2_ANALYSIS}.bedPeak",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"


# ===========================================================================
# TestDAGStructure
# ===========================================================================

class TestDAGStructure:
    """Test DAG generation and structure for ATAC-seq rules."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_atac_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG contains ATAC-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("atac_shift_bam", "make_coverage_atac"):
            assert rule in dag_output, f"Rule '{rule}' not found in ATAC DAG"

    def test_dag_contains_shared_chip_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG contains shared ChIP rules used by ATAC (filter_bam_pe, dispatch_final_bam)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("filter_bam_pe", "dispatch_final_bam"):
            assert rule in dag_output, f"Shared rule '{rule}' not found in ATAC DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WT_REP1_BIGWIG)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert "->" in result.stdout, "DAG should contain rule dependencies"


# ===========================================================================
# TestErrorHandling
# ===========================================================================

class TestErrorHandling:
    """Test error handling for invalid configurations."""

    def test_missing_sample_target(self, snakemake_available, repo_root, test_options):
        """Test that requesting a target for a non-existent sample fails gracefully."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/ATAC/tracks/coverage__final__nonexistent_sample.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for non-existent sample"

    def test_invalid_env_target(self, snakemake_available, repo_root, test_options):
        """Test that requesting a ChIP target fails (only ATAC samples in config)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/ChIP/tracks/coverage__final__WT_ATAC_rep1.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for wrong environment target"
