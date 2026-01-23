"""
Integration tests for ONT methylation workflow using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for ONT samples without actually executing the rules.
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
    config_path = repo_root / "tests" / "integration" / "data" / "test_config_ont.yaml"
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


class TestONTDryRunBasic:
    """Basic dry-run tests for ONT methylation workflow."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_config_file_exists(self, test_config):
        """Test that the test config file exists."""
        assert Path(test_config).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_ont.tsv"
        assert sample_file.exists()


class TestONTModBAMWorkflow:
    """Test ONT modBAM workflow dry-run."""

    @pytest.fixture
    def ont_modbam_target(self):
        """Return target for ONT modBAM bigwig output."""
        return "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__CG.bw"

    def test_ont_modbam_dryrun_succeeds(self, snakemake_available, repo_root, test_config, ont_modbam_target):
        """Test that dry-run succeeds for ONT modBAM sample."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, ont_modbam_target)

        # Check that dry-run completed successfully
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_ont_modbam_includes_expected_rules(self, snakemake_available, repo_root, test_config, ont_modbam_target):
        """Test that ONT modBAM workflow includes expected rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, ont_modbam_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        # Combine stdout and stderr (Snakemake outputs to stderr)
        output = result.stdout + result.stderr

        # Check for ONT-specific rules
        expected_rules = [
            "get_modbam",
            "align_modbam",
            "modkit_pileup",
            "split_bedmethyl_by_context",
            "make_ont_bigwig_files",
            "make_modkit_context_beds"
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected ONT rule '{rule}' not found in dry-run output"

    def test_ont_modbam_excludes_bismark_rules(self, snakemake_available, repo_root, test_config, ont_modbam_target):
        """Test that ONT workflow does not trigger Bismark rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, ont_modbam_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check that Bismark rules are NOT present
        bismark_rules = [
            "bismark_map_pe",
            "bismark_map_se",
            "make_bismark_indices"
        ]

        for rule in bismark_rules:
            assert rule not in output, f"Bismark rule '{rule}' should not be in ONT workflow"

    def test_ont_modbam_all_contexts(self, snakemake_available, repo_root, test_config):
        """Test that all three methylation contexts are generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        contexts = ["CG", "CHG", "CHH"]

        for context in contexts:
            target = f"results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__{context}.bw"
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
            "make_ont_bigwig_files",
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


class TestONTDMRWorkflow:
    """Test ONT DMR calling workflow dry-run."""

    @pytest.fixture
    def dmr_target(self):
        """Return target for DMR analysis output."""
        return "results/mC/DMRs/summary__mC__WT__leaf__ONT__merged__test_genome__vs__mC__mutant__leaf__ONT__merged__test_genome__DMRs.txt"

    def test_dmr_dryrun_succeeds(self, snakemake_available, repo_root, test_config, dmr_target):
        """Test that dry-run succeeds for DMR analysis."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmr_target)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_dmr_includes_modkit_dmr_rule(self, snakemake_available, repo_root, test_config, dmr_target):
        """Test that DMR workflow uses modkit for ONT samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, dmr_target, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check for modkit DMR rule
        assert "call_DMRs_modkit" in output, "Expected modkit DMR rule for ONT samples"


class TestDAGStructure:
    """Test DAG generation and structure."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__CG.bw"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_ont_rules(self, snakemake_available, repo_root, test_config):
        """Test that DAG contains ONT-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__CG.bw"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        # DAG output is in DOT format
        dag_output = result.stdout

        # Check for ONT rule nodes in DAG
        expected_rules = [
            "get_modbam",
            "align_modbam",
            "modkit_pileup",
            "split_bedmethyl_by_context",
            "make_ont_bigwig_files"
        ]

        for rule in expected_rules:
            assert rule in dag_output, f"Rule '{rule}' not found in DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_config):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__CG.bw"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Check that rule dependencies are present
        # format: "rule1" -> "rule2"
        assert "->" in dag_output, "DAG should contain rule dependencies"


class TestWildcardResolution:
    """Test wildcard resolution for ONT samples."""

    def test_ont_sample_wildcards_resolve(self, snakemake_available, repo_root, test_config):
        """Test that wildcards are correctly resolved for ONT samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test various wildcard combinations
        targets = [
            "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__CG.bw",
            "results/mC/tracks/mC__WT__leaf__ONT__rep2__test_genome__CHG.bw",
            "results/mC/tracks/mC__mutant__leaf__ONT__rep1__test_genome__CHH.bw",
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
        target = "results/mC/tracks/mC__WT__leaf__ONT__rep1__nonexistent_genome__CG.bw"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        # Should fail because reference genome is not in config
        assert result.returncode != 0, "Should fail for non-existent reference genome"

    def test_invalid_context(self, snakemake_available, repo_root, test_config):
        """Test that invalid methylation context fails."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request target with invalid context
        target = "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__INVALID.bw"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        # Should fail because INVALID is not a valid context
        assert result.returncode != 0, "Should fail for invalid methylation context"


class TestAllMCTarget:
    """Test the all_mc rule with ONT samples."""

    def test_all_mc_rule_with_ont(self, snakemake_available, repo_root, test_config):
        """Test that all_mc rule works with ONT samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Test the all_mc checkpoint rule
        target = "results/mC/chkpts/mC_analysis__test_ont__test_genome.done"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        assert result.returncode == 0, f"all_mc rule failed with ONT samples: {result.stderr}"

    def test_all_mc_includes_ont_outputs(self, snakemake_available, repo_root, test_config):
        """Test that all_mc includes ONT-specific outputs."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/chkpts/mC_analysis__test_ont__test_genome.done"
        result = run_snakemake_dryrun(repo_root, test_config, target, ["--printshellcmds"])

        assert result.returncode == 0, f"all_mc rule failed: {result.stderr}"

        output = result.stdout + result.stderr

        # Check for ONT summary outputs
        assert "summary__mC__WT__leaf__ONT__rep1__test_genome.txt" in output or "modkit_summary" in output, \
            "all_mc should include ONT summary outputs"


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

    def test_context_beds_are_dependencies(self, snakemake_available, repo_root, test_config):
        """Test that context BED files are dependencies for splitting."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        target = "results/mC/ont/context__mC__WT__leaf__ONT__rep1__test_genome__CG.bed.gz"
        result = run_snakemake_dag(repo_root, test_config, target)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Check that modkit context bed generation is in the DAG
        assert "make_modkit_context_beds" in dag_output, \
            "Context BED generation should be a dependency"


class TestMultipleReplicates:
    """Test handling of multiple replicates."""

    def test_multiple_ont_replicates(self, snakemake_available, repo_root, test_config):
        """Test that multiple ONT replicates are processed."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # Request both replicates
        targets = [
            "results/mC/tracks/mC__WT__leaf__ONT__rep1__test_genome__CG.bw",
            "results/mC/tracks/mC__WT__leaf__ONT__rep2__test_genome__CG.bw"
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_config, target)
            assert result.returncode == 0, f"Replicate processing failed for {target}: {result.stderr}"

    def test_merged_replicates_dmr(self, snakemake_available, repo_root, test_config):
        """Test that merged replicates work for DMR calling."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        # DMR analysis uses merged replicates
        target = "results/mC/DMRs/summary__mC__WT__leaf__ONT__merged__test_genome__vs__mC__mutant__leaf__ONT__merged__test_genome__DMRs.txt"
        result = run_snakemake_dryrun(repo_root, test_config, target)

        assert result.returncode == 0, f"Merged replicate DMR failed: {result.stderr}"
