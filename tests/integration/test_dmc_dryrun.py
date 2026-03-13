"""
Integration tests for dmC (direct methylation) workflow using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for dmC samples (from ONT or other direct methylation platforms)
without actually executing the rules.

Test samples (from test_samples_dmc.tsv, new format):
  Sample_ID                 Assay  Genome       Levels
  WT_leaf_dmC_rep1          dmC    test_genome  genotype:WT,tissue:leaf
  WT_leaf_dmC_rep2          dmC    test_genome  genotype:WT,tissue:leaf
  WT_root_bedMethyl_rep1    dmC    test_genome  genotype:WT,tissue:root
  mutant_leaf_dmC_rep1      dmC    test_genome  genotype:mutant,tissue:leaf
  mutant_leaf_dmC_rep2      dmC    test_genome  genotype:mutant,tissue:leaf

Output paths use Sample_ID directly:
  results/mC/tracks/{Sample_ID}__{context}.bw

Analysis-level names (for DMRs) use build_analysis_name():
  {Assay}__{levels_label}__{IP_target}__{Genome} (empty parts omitted)
  e.g. dmC__WT_leaf__test_genome
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_dmc.yaml")


# ---------------------------------------------------------------------------
# Target path constants (derived from test_samples_dmc.tsv)
# ---------------------------------------------------------------------------

# Per-replicate bigwig targets
DMC_MODBAM_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__CG.bw"
DMC_MODBAM_REP2_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep2__CG.bw"
BEDMETHYL_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_root_bedMethyl_rep1__CG.bw"
MUTANT_TARGET = f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_dmC_rep1__CG.bw"

# Analysis-level names (empty parts omitted)
WT_LEAF_ANALYSIS = "dmC__WT_leaf__test_genome"
MUTANT_LEAF_ANALYSIS = "dmC__mutant_leaf__test_genome"

# DMR target
DMR_TARGET = f"{_OUTPUT_DIR}/mC/DMRs/summary__{WT_LEAF_ANALYSIS}__vs__{MUTANT_LEAF_ANALYSIS}__DMRs.txt"


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_dmc.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


class TestDmcDryRunBasic:
    """Basic dry-run tests for dmC (direct methylation) workflow."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_dmc.tsv"
        assert sample_file.exists()


class TestDmcModBAMWorkflow:
    """Test dmC modBAM workflow dry-run."""

    @pytest.fixture
    def dmc_modbam_target(self):
        """Return target for dmC modBAM bigwig output."""
        return DMC_MODBAM_TARGET

    def test_dmc_modbam_dryrun_succeeds(self, snakemake_available, repo_root, test_options, dmc_modbam_target):
        """Test that dry-run succeeds for dmC modBAM sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, dmc_modbam_target)

        # Check that dry-run completed successfully
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmc_modbam_includes_expected_rules(self, snakemake_available, repo_root, test_options, dmc_modbam_target):
        """Test that dmC modBAM workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, dmc_modbam_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        # Combine stdout and stderr (Snakemake outputs to stderr)
        output = result.stdout + result.stderr

        # Check for dmC-specific rules
        expected_rules = [
            "get_dmc_input",
            "dmc_input_checkpoint",
            "convert_bedmethyl_to_cx_report",
            "make_mc_bigwig_files"
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected dmC rule '{rule}' not found in dry-run output"

    def test_dmc_modbam_excludes_bismark_rules(self, snakemake_available, repo_root, test_options, dmc_modbam_target):
        """Test that dmC workflow does not trigger Bismark rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, dmc_modbam_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check that Bismark rules are NOT present
        bismark_rules = [
            "bismark_map_pe",
            "bismark_map_se",
            "make_bismark_indices"
        ]

        for rule in bismark_rules:
            assert rule not in output, f"Bismark rule '{rule}' should not be in dmC workflow"

    def test_dmc_modbam_all_contexts(self, snakemake_available, repo_root, test_options):
        """Test that all three methylation contexts are generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        contexts = ["CG", "CHG", "CHH"]

        for context in contexts:
            target = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__{context}.bw"
            result = run_snakemake_dryrun(repo_root, test_options, target)

            assert result.returncode == 0, f"Dry-run failed for {context} context: {result.stderr}"


class TestBedMethylWorkflow:
    """Test bedMethyl input workflow dry-run."""

    @pytest.fixture
    def bedmethyl_target(self):
        """Return target for bedMethyl bigwig output."""
        return BEDMETHYL_TARGET

    def test_bedmethyl_dryrun_succeeds(self, snakemake_available, repo_root, test_options, bedmethyl_target):
        """Test that dry-run succeeds for bedMethyl sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, bedmethyl_target)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_bedmethyl_includes_expected_rules(self, snakemake_available, repo_root, test_options, bedmethyl_target):
        """Test that bedMethyl workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, bedmethyl_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check for bedMethyl-specific rules
        # bedMethyl uses get_dmc_input (auto-detects input type),
        # then copy_bedmethyl_input, merge_pileup_sources, convert to CX_report
        expected_rules = [
            "get_dmc_input",
            "dmc_input_checkpoint",
            "convert_bedmethyl_to_cx_report",
            "make_mc_bigwig_files"
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected bedMethyl rule '{rule}' not found in dry-run output"

    def test_bedmethyl_skips_alignment_rules(self, snakemake_available, repo_root, test_options, bedmethyl_target):
        """Test that bedMethyl workflow skips alignment steps."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, bedmethyl_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # bedMethyl should not trigger alignment or pileup from BAM
        skipped_rules = [
            "prepare_modbam_for_pileup",
            "modkit_pileup_dmc"
        ]

        for rule in skipped_rules:
            assert rule not in output, f"Rule '{rule}' should be skipped for bedMethyl input"


class TestDmcDMRWorkflow:
    """Test dmC DMR calling workflow dry-run.

    Note: DMR workflow needs refactoring to work with the new modkit_pileup --motif
    approach. The rule currently expects combined bedMethyl files but we now generate
    context-specific files directly.
    """

    @pytest.fixture
    def dmr_target(self):
        """Return target for DMR analysis output.

        Note: DMR targets use analysis-level names (without replicate).
        """
        return DMR_TARGET

    def test_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_options, dmr_target):
        """Test that dry-run succeeds for DMR analysis."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, dmr_target)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmr_uses_dmrcaller_by_default(self, snakemake_available, repo_root, test_options, dmr_target):
        """Test that DMR workflow uses DMRcaller (same as Bismark) for dmC samples by default."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, dmr_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Default uses DMRcaller (call_DMRs_pairwise) with bedMethyl-to-CX_report conversion
        assert "call_DMRs_pairwise" in output, "Expected DMRcaller rule for dmC samples by default"
        assert "convert_bedmethyl_to_cx_report" in output, "Expected bedMethyl conversion for DMRcaller"
        assert "R_call_DMRs.R" in output, "Expected DMRcaller R script"


class TestDAGStructure:
    """Test DAG generation and structure."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_options, DMC_MODBAM_TARGET)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_dmc_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG contains dmC-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_options, DMC_MODBAM_TARGET)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        # DAG output is in DOT format
        dag_output = result.stdout

        # Check for dmC rule nodes in DAG
        expected_rules = [
            "get_dmc_input",
            "dmc_input_checkpoint",
            "convert_bedmethyl_to_cx_report",
            "make_mc_bigwig_files"
        ]

        for rule in expected_rules:
            assert rule in dag_output, f"Rule '{rule}' not found in DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_options, DMC_MODBAM_TARGET)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Check that rule dependencies are present
        # format: "rule1" -> "rule2"
        assert "->" in dag_output, "DAG should contain rule dependencies"


class TestWildcardResolution:
    """Test wildcard resolution for dmC samples."""

    def test_dmc_sample_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards are correctly resolved for dmC samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test various wildcard combinations
        targets = [
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__CG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep2__CHG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_dmC_rep1__CHH.bw",
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_bedmethyl_sample_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards are correctly resolved for bedMethyl samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, BEDMETHYL_TARGET)

        assert result.returncode == 0, f"Wildcard resolution failed for bedMethyl: {result.stderr}"


class TestErrorHandling:
    """Test error handling for invalid configurations."""

    def test_missing_sample_target(self, snakemake_available, repo_root, test_options):
        """Test that requesting a target for a non-existent sample fails gracefully."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request target with a sample name not in the sample sheet
        target = f"{_OUTPUT_DIR}/mC/tracks/nonexistent_sample__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)

        # Should fail because sample is not in the sample sheet
        assert result.returncode != 0, "Should fail for non-existent sample"

    def test_invalid_context(self, snakemake_available, repo_root, test_options):
        """Test that invalid methylation context fails."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request target with invalid context
        target = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__INVALID.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)

        # Should fail because INVALID is not a valid context
        assert result.returncode != 0, "Should fail for invalid methylation context"


class TestAllMCTarget:
    """Test the all_mc rule with dmC samples.

    Note: all_mc rule needs updates to handle merged replicates properly.
    Currently skipped until merged replicate handling is implemented.
    """

    @pytest.mark.skip(reason="all_mc rule needs merged replicate handling for dmC")
    def test_all_mc_rule_with_dmc(self, snakemake_available, repo_root, test_options):
        """Test that all_mc rule works with dmC samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test the all_mc checkpoint rule
        target = f"{_OUTPUT_DIR}/mC/chkpts/mC_analysis__test_dmc__test_genome.done"
        result = run_snakemake_dryrun(repo_root, test_options, target)

        assert result.returncode == 0, f"all_mc rule failed with dmC samples: {result.stderr}"

    @pytest.mark.skip(reason="all_mc rule needs merged replicate handling for dmC")
    def test_all_mc_includes_dmc_outputs(self, snakemake_available, repo_root, test_options):
        """Test that all_mc includes dmC-specific outputs."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = f"{_OUTPUT_DIR}/mC/chkpts/mC_analysis__test_dmc__test_genome.done"
        result = run_snakemake_dryrun(repo_root, test_options, target, ["--printshellcmds"])

        assert result.returncode == 0, f"all_mc rule failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check for dmC summary outputs
        assert "modkit_summary" in output or "WT_leaf_dmC_rep1" in output, \
            "all_mc should include dmC outputs"


class TestCxReportConversion:
    """Test CX_report conversion for dmC samples."""

    def test_cx_report_generated_for_dmc(self, snakemake_available, repo_root, test_options):
        """Test that CX_report conversion is in the DAG for dmC bigwig generation."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, DMC_MODBAM_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr
        assert "convert_bedmethyl_to_cx_report" in output, \
            "CX_report conversion should be in the DAG for dmC samples"

    def test_cx_report_generated_for_bedmethyl(self, snakemake_available, repo_root, test_options):
        """Test that CX_report conversion is in the DAG for bedMethyl bigwig generation."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_options, BEDMETHYL_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr
        assert "convert_bedmethyl_to_cx_report" in output, \
            "CX_report conversion should be in the DAG for bedMethyl samples"


class TestMultipleReplicates:
    """Test handling of multiple replicates."""

    def test_multiple_dmc_replicates(self, snakemake_available, repo_root, test_options):
        """Test that multiple dmC replicates are processed."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request both replicates
        targets = [
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__CG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep2__CG.bw"
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Replicate processing failed for {target}: {result.stderr}"

    def test_merged_replicates_dmr(self, snakemake_available, repo_root, test_options):
        """Test that DMR calling works with multiple replicates (merged automatically)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # DMR analysis uses analysis-level names (replicates merged automatically in the rule)
        result = run_snakemake_dryrun(repo_root, test_options, DMR_TARGET)

        assert result.returncode == 0, f"Merged replicate DMR failed: {result.stderr}"
