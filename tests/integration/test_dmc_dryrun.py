"""
Integration tests for dmC (direct methylation) workflow using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for dmC samples (from ONT or other direct methylation platforms)
without actually executing the rules.
"""

import pytest
import subprocess
import re
import os
from pathlib import Path


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def repo_root():
    """Get the repository root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="module")
def test_config(repo_root):
    """Get the path to the test config file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_config_dmc.yaml"
    assert config_path.exists(), f"Test config not found at {config_path}"
    return str(config_path)


@pytest.fixture(scope="module")
def snakemake_available():
    """Check if snakemake is available."""
    try:
        result = subprocess.run(
            ["snakemake", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def run_snakemake_dryrun(repo_root, config_file, target=None, extra_args=None):
    """
    Run snakemake in dry-run mode with the given config.

    Args:
        repo_root: Path to repository root
        config_file: Path to config file
        target: Optional target file/rule to request
        extra_args: Optional list of additional arguments

    Returns:
        subprocess.CompletedProcess object
    """
    cmd = [
        "snakemake",
        "--dry-run",
        "--configfile", config_file,
        "--cores", "1",
    ]

    if target:
        cmd.append(target)

    # Add --quiet at end to avoid it consuming the target (Snakemake 9 behavior)
    cmd.append("--quiet")
    cmd.append("progress")

    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(
        cmd,
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=60
    )

    return result


def run_snakemake_dag(repo_root, config_file, target=None):
    """
    Generate the Snakemake DAG.

    Args:
        repo_root: Path to repository root
        config_file: Path to config file
        target: Optional target file/rule to request

    Returns:
        subprocess.CompletedProcess object
    """
    cmd = [
        "snakemake",
        "--dag",
        "--configfile", config_file,
        "--cores", "1"
    ]

    if target:
        cmd.append(target)

    result = subprocess.run(
        cmd,
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=60
    )

    return result


class TestDmcDryRunBasic:
    """Basic dry-run tests for dmC (direct methylation) workflow."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_config_file_exists(self, test_config):
        """Test that the test config file exists."""
        assert Path(test_config).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_dmc.tsv"
        assert sample_file.exists()


class TestDmcModBAMWorkflow:
    """Test dmC modBAM workflow dry-run."""

    @pytest.fixture
    def dmc_modbam_target(self):
        """Return target for dmC modBAM bigwig output."""
        return "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw"

    def test_dmc_modbam_dryrun_succeeds(self, snakemake_available, repo_root, test_config, dmc_modbam_target):
        """Test that dry-run succeeds for dmC modBAM sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmc_modbam_target)

        # Check that dry-run completed successfully
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmc_modbam_includes_expected_rules(self, snakemake_available, repo_root, test_config, dmc_modbam_target):
        """Test that dmC modBAM workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmc_modbam_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        # Combine stdout and stderr (Snakemake outputs to stderr)
        output = result.stdout + result.stderr

        # Check for dmC-specific rules
        # Note: modkit_pileup now does context filtering directly with --motif,
        # so split_bedmethyl_by_context and make_modkit_context_beds are not needed for dmC
        expected_rules = [
            "get_modbam",
            "align_modbam",
            "modkit_pileup",
            "make_dmc_bigwig_files"
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected dmC rule '{rule}' not found in dry-run output"

    def test_dmc_modbam_excludes_bismark_rules(self, snakemake_available, repo_root, test_config, dmc_modbam_target):
        """Test that dmC workflow does not trigger Bismark rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmc_modbam_target, ["--printshellcmds"])

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

    def test_dmc_modbam_all_contexts(self, snakemake_available, repo_root, test_config):
        """Test that all three methylation contexts are generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        contexts = ["CG", "CHG", "CHH"]

        for context in contexts:
            target = f"results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__{context}.bw"
            result = run_snakemake_dryrun(repo_root, test_config, target)

            assert result.returncode == 0, f"Dry-run failed for {context} context: {result.stderr}"


class TestBedMethylWorkflow:
    """Test bedMethyl input workflow dry-run."""

    @pytest.fixture
    def bedmethyl_target(self):
        """Return target for bedMethyl bigwig output."""
        return "results/mC/tracks/mC__WT__root__bedMethyl__rep1__test_genome__CG.bw"

    def test_bedmethyl_dryrun_succeeds(self, snakemake_available, repo_root, test_config, bedmethyl_target):
        """Test that dry-run succeeds for bedMethyl sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, bedmethyl_target)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_bedmethyl_includes_expected_rules(self, snakemake_available, repo_root, test_config, bedmethyl_target):
        """Test that bedMethyl workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, bedmethyl_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check for bedMethyl-specific rules
        expected_rules = [
            "get_bedmethyl",
            "copy_bedmethyl_for_pileup",
            "split_bedmethyl_by_context",
            "make_dmc_bigwig_files",
            "make_modkit_context_beds"
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected bedMethyl rule '{rule}' not found in dry-run output"

    def test_bedmethyl_skips_alignment_rules(self, snakemake_available, repo_root, test_config, bedmethyl_target):
        """Test that bedMethyl workflow skips alignment steps."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, bedmethyl_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # bedMethyl should not trigger alignment or pileup from BAM
        skipped_rules = [
            "align_modbam",
            "modkit_pileup"
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
        return "results/mC/DMRs/summary__mC__WT__leaf__dmC__test_genome__vs__mC__mutant__leaf__dmC__test_genome__DMRs.txt"

    def test_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_config, dmr_target):
        """Test that dry-run succeeds for DMR analysis."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmr_target)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmr_uses_dmrcaller_by_default(self, snakemake_available, repo_root, test_config, dmr_target):
        """Test that DMR workflow uses DMRcaller (same as Bismark) for dmC samples by default."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmr_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Default uses DMRcaller (call_DMRs_pairwise) with bedMethyl-to-CX_report conversion
        assert "call_DMRs_pairwise" in output, "Expected DMRcaller rule for dmC samples by default"
        assert "convert_bedmethyl_to_cx_report" in output, "Expected bedMethyl conversion for DMRcaller"
        assert "R_call_DMRs.R" in output, "Expected DMRcaller R script"


class TestDAGStructure:
    """Test DAG generation and structure."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_dmc_rules(self, snakemake_available, repo_root, test_config):
        """Test that DAG contains dmC-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        # DAG output is in DOT format
        dag_output = result.stdout

        # Check for dmC rule nodes in DAG
        # Note: modkit_pileup now does context filtering directly with --motif,
        # so split_bedmethyl_by_context is not needed for dmC samples
        expected_rules = [
            "get_modbam",
            "align_modbam",
            "modkit_pileup",
            "make_dmc_bigwig_files"
        ]

        for rule in expected_rules:
            assert rule in dag_output, f"Rule '{rule}' not found in DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_config):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Check that rule dependencies are present
        # format: "rule1" -> "rule2"
        assert "->" in dag_output, "DAG should contain rule dependencies"


class TestWildcardResolution:
    """Test wildcard resolution for dmC samples."""

    def test_dmc_sample_wildcards_resolve(self, snakemake_available, repo_root, test_config):
        """Test that wildcards are correctly resolved for dmC samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test various wildcard combinations
        targets = [
            "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw",
            "results/mC/tracks/mC__WT__leaf__dmC__rep2__test_genome__CHG.bw",
            "results/mC/tracks/mC__mutant__leaf__dmC__rep1__test_genome__CHH.bw",
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_config, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_bedmethyl_sample_wildcards_resolve(self, snakemake_available, repo_root, test_config):
        """Test that wildcards are correctly resolved for bedMethyl samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__root__bedMethyl__rep1__test_genome__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        assert result.returncode == 0, f"Wildcard resolution failed for bedMethyl: {result.stderr}"


class TestErrorHandling:
    """Test error handling for invalid configurations."""

    def test_missing_reference_genome(self, snakemake_available, repo_root, test_config):
        """Test that requesting non-existent reference genome fails gracefully."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request target with non-existent reference genome
        target = "results/mC/tracks/mC__WT__leaf__dmC__rep1__nonexistent_genome__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        # Should fail because reference genome is not in config
        assert result.returncode != 0, "Should fail for non-existent reference genome"

    def test_invalid_context(self, snakemake_available, repo_root, test_config):
        """Test that invalid methylation context fails."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request target with invalid context
        target = "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__INVALID.bw"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        # Should fail because INVALID is not a valid context
        assert result.returncode != 0, "Should fail for invalid methylation context"


class TestAllMCTarget:
    """Test the all_mc rule with dmC samples.

    Note: all_mc rule needs updates to handle merged replicates properly.
    Currently skipped until merged replicate handling is implemented.
    """

    @pytest.mark.skip(reason="all_mc rule needs merged replicate handling for dmC")
    def test_all_mc_rule_with_dmc(self, snakemake_available, repo_root, test_config):
        """Test that all_mc rule works with dmC samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test the all_mc checkpoint rule
        target = "results/mC/chkpts/mC_analysis__test_dmc__test_genome.done"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        assert result.returncode == 0, f"all_mc rule failed with dmC samples: {result.stderr}"

    @pytest.mark.skip(reason="all_mc rule needs merged replicate handling for dmC")
    def test_all_mc_includes_dmc_outputs(self, snakemake_available, repo_root, test_config):
        """Test that all_mc includes dmC-specific outputs."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/chkpts/mC_analysis__test_dmc__test_genome.done"
        result = run_snakemake_dryrun(repo_root, test_config, target, ["--printshellcmds"])

        assert result.returncode == 0, f"all_mc rule failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check for dmC summary outputs
        assert "summary__mC__WT__leaf__dmC__rep1__test_genome.txt" in output or "modkit_summary" in output, \
            "all_mc should include dmC summary outputs"


class TestContextBedGeneration:
    """Test methylation context BED file generation."""

    def test_context_beds_generation(self, snakemake_available, repo_root, test_config):
        """Test that context BED files are generated for reference genome."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        contexts = ["CG", "CHG", "CHH"]

        for context in contexts:
            target = f"genomes/test_genome/modkit_{context}.bed.gz"
            result = run_snakemake_dryrun(repo_root, test_config, target)

            assert result.returncode == 0, f"Context BED generation failed for {context}: {result.stderr}"

    def test_context_beds_are_dependencies_for_bedmethyl(self, snakemake_available, repo_root, test_config):
        """Test that context BED files are dependencies for bedMethyl splitting.

        Note: For dmC samples, modkit_pileup uses --motif filtering directly,
        so context BEDs are only needed for bedMethyl sample_type inputs.
        """
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test bedMethyl sample (not dmC) - context BEDs should be dependencies
        target = "results/mC/dmc/context__mC__WT__root__bedMethyl__rep1__test_genome__CG.bed.gz"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Check that modkit context bed generation is in the DAG for bedMethyl samples
        assert "make_modkit_context_beds" in dag_output, \
            "Context BED generation should be a dependency for bedMethyl samples"


class TestMultipleReplicates:
    """Test handling of multiple replicates."""

    def test_multiple_dmc_replicates(self, snakemake_available, repo_root, test_config):
        """Test that multiple dmC replicates are processed."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request both replicates
        targets = [
            "results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw",
            "results/mC/tracks/mC__WT__leaf__dmC__rep2__test_genome__CG.bw"
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_config, target)
            assert result.returncode == 0, f"Replicate processing failed for {target}: {result.stderr}"

    def test_merged_replicates_dmr(self, snakemake_available, repo_root, test_config):
        """Test that DMR calling works with multiple replicates (merged automatically)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # DMR analysis uses analysis-level names (replicates merged automatically in the rule)
        target = "results/mC/DMRs/summary__mC__WT__leaf__dmC__test_genome__vs__mC__mutant__leaf__dmC__test_genome__DMRs.txt"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        assert result.returncode == 0, f"Merged replicate DMR failed: {result.stderr}"
