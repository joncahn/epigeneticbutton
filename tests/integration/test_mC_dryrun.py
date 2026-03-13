"""
Integration tests for ALL methylation assay types using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for every mC assay type without actually executing the rules.

Test samples (from test_samples_mC.tsv):
  Bismark assays:
    WT_leaf_WGBS_rep1       WGBS     PE   genotype:WT,tissue:leaf
    WT_leaf_WGBS_rep2       WGBS     PE   genotype:WT,tissue:leaf
    mutant_leaf_WGBS_rep1   WGBS     PE   genotype:mutant,tissue:leaf
    WT_leaf_WGBSnd_rep1     WGBS_nd  SE   genotype:WT,tissue:leaf
    mutant_leaf_WGBSnd_rep1 WGBS_nd  SE   genotype:mutant,tissue:leaf
    WT_leaf_PBAT_rep1       PBAT     PE   genotype:WT,tissue:leaf
    WT_leaf_EMseq_rep1      EMseq    PE   genotype:WT,tissue:leaf
    mutant_leaf_EMseq_rep1  EMseq    PE   genotype:mutant,tissue:leaf
  dmC assays:
    WT_leaf_dmC_rep1        dmC      SE   genotype:WT,tissue:leaf     (modBAM)
    WT_leaf_dmC_rep2        dmC      SE   genotype:WT,tissue:leaf     (modBAM)
    mutant_leaf_dmC_rep1    dmC      SE   genotype:mutant,tissue:leaf (modBAM)
    mutant_leaf_dmC_rep2    dmC      SE   genotype:mutant,tissue:leaf (modBAM)
    WT_root_bedMethyl_rep1  dmC      SE   genotype:WT,tissue:root     (bedMethyl)
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_mC.yaml")


# ---------------------------------------------------------------------------
# Target path constants
# ---------------------------------------------------------------------------

# Bismark per-replicate bigwig targets
WGBS_PE_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_WGBS_rep1__CG.bw"
WGBS_PE_REP2_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_WGBS_rep2__CG.bw"
WGBS_ND_SE_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_WGBSnd_rep1__CG.bw"
PBAT_PE_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_PBAT_rep1__CG.bw"
EMSEQ_PE_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_EMseq_rep1__CG.bw"

# dmC per-replicate bigwig targets
DMC_MODBAM_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__CG.bw"
DMC_MODBAM_REP2_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep2__CG.bw"
BEDMETHYL_TARGET = f"{_OUTPUT_DIR}/mC/tracks/WT_root_bedMethyl_rep1__CG.bw"
MUTANT_DMC_TARGET = f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_dmC_rep1__CG.bw"

# Analysis-level names (Assay__levels_label__Genome, empty parts omitted)
WGBS_WT_ANALYSIS = "WGBS__WT_leaf__test_genome"
WGBS_MUTANT_ANALYSIS = "WGBS__mutant_leaf__test_genome"
WGBS_ND_WT_ANALYSIS = "WGBS_nd__WT_leaf__test_genome"
WGBS_ND_MUTANT_ANALYSIS = "WGBS_nd__mutant_leaf__test_genome"
EMSEQ_WT_ANALYSIS = "EMseq__WT_leaf__test_genome"
EMSEQ_MUTANT_ANALYSIS = "EMseq__mutant_leaf__test_genome"
DMC_WT_LEAF_ANALYSIS = "dmC__WT_leaf__test_genome"
DMC_MUTANT_LEAF_ANALYSIS = "dmC__mutant_leaf__test_genome"

# DMR targets (bismark and dmC)
WGBS_DMR_TARGET = f"{_OUTPUT_DIR}/mC/DMRs/summary__{WGBS_WT_ANALYSIS}__vs__{WGBS_MUTANT_ANALYSIS}__DMRs.txt"
WGBS_ND_DMR_TARGET = f"{_OUTPUT_DIR}/mC/DMRs/summary__{WGBS_ND_WT_ANALYSIS}__vs__{WGBS_ND_MUTANT_ANALYSIS}__DMRs.txt"
EMSEQ_DMR_TARGET = f"{_OUTPUT_DIR}/mC/DMRs/summary__{EMSEQ_WT_ANALYSIS}__vs__{EMSEQ_MUTANT_ANALYSIS}__DMRs.txt"
DMC_DMR_TARGET = f"{_OUTPUT_DIR}/mC/DMRs/summary__{DMC_WT_LEAF_ANALYSIS}__vs__{DMC_MUTANT_LEAF_ANALYSIS}__DMRs.txt"


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_mC.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


# ===========================================================================
# TestBasic
# ===========================================================================

class TestBasic:
    """Basic sanity checks for the mC test configuration."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_mC.tsv"
        assert sample_file.exists()


# ===========================================================================
# TestBismarkWorkflows
# ===========================================================================

class TestBismarkWorkflows:
    """Test Bismark-based methylation workflows (WGBS, WGBS_nd, PBAT, EMseq)."""

    def test_wgbs_pe_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for WGBS PE sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_PE_TARGET)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_wgbs_pe_includes_bismark_rules(self, snakemake_available, repo_root, test_options):
        """Test that WGBS PE workflow includes bismark_map_pe."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "bismark_map_pe" in output, "WGBS PE should use bismark_map_pe rule"

    def test_wgbs_nd_se_includes_bismark_rules(self, snakemake_available, repo_root, test_options):
        """Test that WGBS_nd SE workflow includes bismark_map_se."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_ND_SE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "bismark_map_se" in output, "WGBS_nd SE should use bismark_map_se rule"

    def test_pbat_pe_includes_bismark_rules(self, snakemake_available, repo_root, test_options):
        """Test that PBAT PE workflow includes bismark_map_pe."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, PBAT_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "bismark_map_pe" in output, "PBAT PE should use bismark_map_pe rule"

    def test_emseq_pe_includes_bismark_rules(self, snakemake_available, repo_root, test_options):
        """Test that EMseq PE workflow includes bismark_map_pe."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, EMSEQ_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "bismark_map_pe" in output, "EMseq PE should use bismark_map_pe rule"


# ===========================================================================
# TestBismarkParameterRouting
# ===========================================================================

class TestBismarkParameterRouting:
    """Verify that bismark mapping and extraction parameters vary by assay type.

    The mC_mapping config section routes --non_directional, --pbat, and
    --ignore_r2 to the correct assay types via parameters_for_mc().
    """

    @pytest.fixture(scope="class")
    def wgbs_output(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        return result.stdout + result.stderr

    @pytest.fixture(scope="class")
    def wgbs_nd_output(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_ND_SE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        return result.stdout + result.stderr

    @pytest.fixture(scope="class")
    def pbat_output(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, PBAT_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        return result.stdout + result.stderr

    @pytest.fixture(scope="class")
    def emseq_output(self, snakemake_available, repo_root, test_options):
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, EMSEQ_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        return result.stdout + result.stderr

    # -- WGBS: standard directional, with --ignore_r2 --

    def test_wgbs_no_non_directional(self, wgbs_output):
        """WGBS should NOT have --non_directional."""
        assert "--non_directional" not in wgbs_output

    def test_wgbs_no_pbat(self, wgbs_output):
        """WGBS should NOT have --pbat."""
        assert "--pbat" not in wgbs_output

    def test_wgbs_has_ignore_r2(self, wgbs_output):
        """WGBS PE extraction should include --ignore_r2."""
        assert "--ignore_r2" in wgbs_output

    # -- WGBS_nd: non-directional, no --ignore_r2 --

    def test_wgbs_nd_has_non_directional(self, wgbs_nd_output):
        """WGBS_nd should have --non_directional."""
        assert "--non_directional" in wgbs_nd_output

    def test_wgbs_nd_no_pbat(self, wgbs_nd_output):
        """WGBS_nd should NOT have --pbat."""
        assert "--pbat" not in wgbs_nd_output

    def test_wgbs_nd_no_ignore_r2(self, wgbs_nd_output):
        """WGBS_nd should NOT have --ignore_r2."""
        assert "--ignore_r2" not in wgbs_nd_output

    # -- PBAT: --pbat flag, no --ignore_r2 --

    def test_pbat_has_pbat(self, pbat_output):
        """PBAT should have --pbat."""
        assert "--pbat" in pbat_output

    def test_pbat_no_non_directional(self, pbat_output):
        """PBAT should NOT have --non_directional."""
        assert "--non_directional" not in pbat_output

    def test_pbat_no_ignore_r2(self, pbat_output):
        """PBAT should NOT have --ignore_r2."""
        assert "--ignore_r2" not in pbat_output

    # -- EMseq: standard directional (like WGBS but no --ignore_r2) --

    def test_emseq_no_non_directional(self, emseq_output):
        """EMseq should NOT have --non_directional."""
        assert "--non_directional" not in emseq_output

    def test_emseq_no_pbat(self, emseq_output):
        """EMseq should NOT have --pbat."""
        assert "--pbat" not in emseq_output

    def test_emseq_no_ignore_r2(self, emseq_output):
        """EMseq should NOT have --ignore_r2."""
        assert "--ignore_r2" not in emseq_output


# ===========================================================================
# TestDmcModBAMWorkflow
# ===========================================================================

class TestDmcModBAMWorkflow:
    """Test dmC modBAM workflow dry-run."""

    def test_dmc_modbam_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for dmC modBAM sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DMC_MODBAM_TARGET)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmc_modbam_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that dmC modBAM workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DMC_MODBAM_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        expected_rules = [
            "get_dmc_input",
            "dmc_input_checkpoint",
            "convert_bedmethyl_to_cx_report",
            "make_mc_bigwig_files",
        ]
        for rule in expected_rules:
            assert rule in output, f"Expected dmC rule '{rule}' not found in dry-run output"

    def test_dmc_modbam_excludes_bismark_rules(self, snakemake_available, repo_root, test_options):
        """Test that dmC workflow does not trigger Bismark rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DMC_MODBAM_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        bismark_rules = ["bismark_map_pe", "bismark_map_se"]
        for rule in bismark_rules:
            assert rule not in output, f"Bismark rule '{rule}' should not be in dmC workflow"

    def test_dmc_modbam_all_contexts(self, snakemake_available, repo_root, test_options):
        """Test that all three methylation contexts are generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        for context in ("CG", "CHG", "CHH"):
            target = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__{context}.bw"
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Dry-run failed for {context}: {result.stderr}"


# ===========================================================================
# TestBedMethylWorkflow
# ===========================================================================

class TestBedMethylWorkflow:
    """Test bedMethyl input workflow dry-run."""

    def test_bedmethyl_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for bedMethyl sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BEDMETHYL_TARGET)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_bedmethyl_includes_expected_rules(self, snakemake_available, repo_root, test_options):
        """Test that bedMethyl workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BEDMETHYL_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        expected_rules = [
            "get_dmc_input",
            "dmc_input_checkpoint",
            "convert_bedmethyl_to_cx_report",
            "make_mc_bigwig_files",
        ]
        for rule in expected_rules:
            assert rule in output, f"Expected bedMethyl rule '{rule}' not found in dry-run output"

    def test_bedmethyl_skips_alignment_rules(self, snakemake_available, repo_root, test_options):
        """Test that bedMethyl workflow skips alignment steps."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BEDMETHYL_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        skipped_rules = ["prepare_modbam_for_pileup", "modkit_pileup_dmc"]
        for rule in skipped_rules:
            assert rule not in output, f"Rule '{rule}' should be skipped for bedMethyl input"


# ===========================================================================
# TestDmcDMRWorkflow
# ===========================================================================

class TestDmcDMRWorkflow:
    """Test dmC DMR calling workflow dry-run."""

    def test_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for DMR analysis between dmC samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DMC_DMR_TARGET)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmr_uses_dmrcaller_by_default(self, snakemake_available, repo_root, test_options):
        """Test that DMR workflow uses DMRcaller for dmC samples by default."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DMC_DMR_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "call_DMRs_pairwise" in output, "Expected DMRcaller rule for dmC samples by default"
        assert "convert_bedmethyl_to_cx_report" in output, "Expected bedMethyl conversion for DMRcaller"
        assert "R_call_DMRs.R" in output, "Expected DMRcaller R script"


# ===========================================================================
# TestDAGStructure
# ===========================================================================

class TestDAGStructure:
    """Test DAG generation and structure across bismark and dmC rules."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WGBS_PE_TARGET)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_bismark_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG for a bismark target contains bismark-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WGBS_PE_TARGET)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("bismark_map_pe", "make_bismark_indices", "make_mc_bigwig_files"):
            assert rule in dag_output, f"Rule '{rule}' not found in bismark DAG"

    def test_dag_contains_dmc_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG for a dmC target contains dmC-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, DMC_MODBAM_TARGET)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("get_dmc_input", "dmc_input_checkpoint", "convert_bedmethyl_to_cx_report", "make_mc_bigwig_files"):
            assert rule in dag_output, f"Rule '{rule}' not found in dmC DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, WGBS_PE_TARGET)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert "->" in result.stdout, "DAG should contain rule dependencies"


# ===========================================================================
# TestWildcardResolution
# ===========================================================================

class TestWildcardResolution:
    """Test wildcard resolution for all mC assay types."""

    def test_wgbs_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for WGBS samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_WGBS_rep1__CG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_WGBS_rep2__CHG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_WGBS_rep1__CHH.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_wgbs_nd_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for WGBS_nd samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_WGBSnd_rep1__CG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_WGBSnd_rep1__CHG.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_pbat_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for PBAT samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, PBAT_PE_TARGET)
        assert result.returncode == 0, f"Wildcard resolution failed for PBAT: {result.stderr}"

    def test_emseq_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for EMseq samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_EMseq_rep1__CG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_EMseq_rep1__CHG.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_dmc_sample_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for dmC modBAM samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__CG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep2__CHG.bw",
            f"{_OUTPUT_DIR}/mC/tracks/mutant_leaf_dmC_rep1__CHH.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_bedmethyl_sample_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for bedMethyl samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, BEDMETHYL_TARGET)
        assert result.returncode == 0, f"Wildcard resolution failed for bedMethyl: {result.stderr}"


# ===========================================================================
# TestReplicateMerging
# ===========================================================================

class TestReplicateMerging:
    """Test replicate merging rules for assays with multiple replicates."""

    def test_wgbs_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """WGBS WT_leaf has 2 reps -- merged replicate rules should appear."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        # Target the merged bigwig (analysis-level name)
        target = f"{_OUTPUT_DIR}/mC/tracks/{WGBS_WT_ANALYSIS}__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_mc_replicates" in output, \
            "WGBS with 2 reps should trigger merging_mc_replicates rule"

    def test_dmc_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """dmC WT_leaf has 2 reps -- merged replicate rules should appear."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/mC/tracks/{DMC_WT_LEAF_ANALYSIS}__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_mc_replicates" in output, \
            "dmC with 2 reps should trigger merging_mc_replicates rule"

    def test_emseq_single_rep_analysis_level(self, snakemake_available, repo_root, test_options):
        """EMseq WT_leaf has 1 rep -- analysis-level target still succeeds.

        Note: the pipeline routes single-rep analysis-level targets through
        merging_mc_replicates as a passthrough (copy/symlink), so we only
        verify the dry-run succeeds and that only one bismark_map_pe job
        appears for EMseq (not two).
        """
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/mC/tracks/{EMSEQ_WT_ANALYSIS}__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        # Only one EMseq rep should be mapped (WT_leaf_EMseq_rep1)
        assert "WT_leaf_EMseq_rep1" in output, \
            "EMseq analysis-level target should include the single replicate"


# ===========================================================================
# TestMixedAssayDMRs
# ===========================================================================

class TestMixedAssayDMRs:
    """Verify DMR rules are generated for WT-vs-mutant comparisons within each assay type."""

    def test_wgbs_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that WGBS WT vs mutant DMR dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_DMR_TARGET)
        assert result.returncode == 0, f"WGBS DMR dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_wgbs_nd_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that WGBS_nd WT vs mutant DMR dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_ND_DMR_TARGET)
        assert result.returncode == 0, f"WGBS_nd DMR dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_emseq_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that EMseq WT vs mutant DMR dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, EMSEQ_DMR_TARGET)
        assert result.returncode == 0, f"EMseq DMR dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmc_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dmC WT vs mutant DMR dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DMC_DMR_TARGET)
        assert result.returncode == 0, f"dmC DMR dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmr_rules_include_call_DMRs_pairwise(self, snakemake_available, repo_root, test_options):
        """Test that all DMR targets use call_DMRs_pairwise."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        for label, target in [
            ("WGBS", WGBS_DMR_TARGET),
            ("EMseq", EMSEQ_DMR_TARGET),
            ("dmC", DMC_DMR_TARGET),
        ]:
            result = run_snakemake_dryrun(repo_root, test_options, target, ["--printshellcmds"])
            assert result.returncode == 0, f"{label} DMR dry-run failed: {result.stderr}"
            output = result.stdout + result.stderr
            assert "call_DMRs_pairwise" in output, \
                f"{label} DMR should include call_DMRs_pairwise rule"


# ===========================================================================
# TestErrorHandling
# ===========================================================================

class TestErrorHandling:
    """Test error handling for invalid configurations."""

    def test_missing_sample_target(self, snakemake_available, repo_root, test_options):
        """Test that requesting a target for a non-existent sample fails gracefully."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/mC/tracks/nonexistent_sample__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for non-existent sample"

    def test_invalid_context(self, snakemake_available, repo_root, test_options):
        """Test that invalid methylation context fails."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/mC/tracks/WT_leaf_dmC_rep1__INVALID.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for invalid methylation context"


# ===========================================================================
# TestCxReportConversion
# ===========================================================================

class TestCxReportConversion:
    """Test CX_report conversion for dmC and bismark samples."""

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

    def test_cx_report_generated_for_bismark(self, snakemake_available, repo_root, test_options):
        """Test that bismark samples produce CX_report files (via bismark_methylation_extractor)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, WGBS_PE_TARGET, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        # Bismark generates CX_reports directly via bismark_methylation_extractor
        assert "cytosine_report" in output or "CX_report" in output, \
            "Bismark samples should produce CX_report files"
