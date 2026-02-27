"""
Integration tests for S. pombe multi-assay workflow using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for all four assay types without actually executing the rules.

Test samples (from test_samples_pombe.tsv):
  ChIP_broad (8): WT/dcr1 x H3K9me2/me3 + 2 Inputs, PE
  ChIP_narrow (3): veg H3K4me3 x 2 reps + Input, SE
  RNAseq (4): WT/dcr1 x 2 reps, SE
  sRNA (3): WT x 2 reps + dcr1 x 1 rep, SE

Output paths use Sample_ID directly for per-replicate files.
Analysis-level names use build_analysis_name():
  {Assay}__{levels_label}__{IP_target}__{Genome}
"""

import pytest
import subprocess
from pathlib import Path


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration


# ---------------------------------------------------------------------------
# Target path constants (derived from test_samples_pombe.tsv)
# ---------------------------------------------------------------------------

# Per-replicate targets
CHIP_BROAD_TARGET = "results/ChIP/tracks/coverage__WT_cell_H3K9me2_rep1.bw"
CHIP_BROAD_TARGET_REP2 = "results/ChIP/tracks/coverage__WT_cell_H3K9me2_rep2.bw"
CHIP_NARROW_TARGET = "results/ChIP/tracks/coverage__veg_cell_H3K4me3_rep1.bw"
RNA_TARGET = "results/RNA/mapped/final__WT_cell_RNA_rep1.bam"
SRNA_TARGET = "results/sRNA/mapped/WT_cell_sRNA_rep1/Results.txt"

# Control sample targets
INPUT_BROAD_TARGET = "results/ChIP/tracks/coverage__WT_cell_Input_rep1.bw"
INPUT_NARROW_TARGET = "results/ChIP/tracks/coverage__veg_cell_Input_rep1.bw"

# Analysis-level names (Assay__levels_label__IP_target__Genome)
CHIP_BROAD_ANALYSIS = "ChIP_broad__WT_cell__H3K9me2__Spombe"
CHIP_NARROW_ANALYSIS = "ChIP_narrow__veg_cell__H3K4me3__Spombe"
RNA_ANALYSIS_WT = "RNAseq__WT_cell____Spombe"
SRNA_ANALYSIS_WT = "sRNA__WT_cell____Spombe"

# Env checkpoint targets
CHIP_CHECKPOINT = "results/ChIP/chkpts/ChIP_analysis__test_pombe__Spombe.done"
RNA_CHECKPOINT = "results/RNA/chkpts/RNA_analysis__test_pombe__Spombe.done"
SRNA_CHECKPOINT = "results/sRNA/chkpts/sRNA_analysis__test_pombe__Spombe.done"
FINAL_CHECKPOINT = "results/combined/chkpts/final_analysis__test_pombe.done"


@pytest.fixture(scope="module")
def repo_root():
    """Get the repository root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture(scope="module")
def test_config(repo_root):
    """Get the path to the test config file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_config_pombe.yaml"
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
        "--quiet", "progress",
    ]

    if extra_args:
        cmd.extend(extra_args)

    # Use -- separator to prevent Snakemake 9 from interpreting targets as options
    if target:
        cmd.append("--")
        cmd.append(target)

    result = subprocess.run(
        cmd,
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=120
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
        "--cores", "1",
    ]

    # Use -- separator to prevent Snakemake 9 from interpreting targets as options
    if target:
        cmd.append("--")
        cmd.append(target)

    result = subprocess.run(
        cmd,
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=120
    )

    return result


class TestPombeDryRunBasic:
    """Basic setup verification tests."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_config_file_exists(self, test_config):
        """Test that the test config file exists."""
        assert Path(test_config).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_pombe.tsv"
        assert sample_file.exists()


class TestChIPBroadWorkflow:
    """Test ChIP-seq broad-peak workflow dry-run (PE samples)."""

    def test_chip_broad_dryrun_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that dry-run succeeds for ChIP broad coverage bigwig."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_BROAD_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_chip_broad_includes_expected_rules(self, snakemake_available, repo_root, test_config):
        """Test that ChIP broad workflow includes PE mapping and broad peak rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_BROAD_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "bowtie2_map_pe",
            "filter_chip_pe",
            "make_coverage_chip",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_chip_broad_excludes_se_mapping(self, snakemake_available, repo_root, test_config):
        """Test that PE ChIP broad samples do not trigger SE mapping rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_BROAD_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "bowtie2_map_se" not in output, \
            "SE mapping rule should not appear for PE ChIP broad sample"

    def test_chip_broad_peak_calling(self, snakemake_available, repo_root, test_config):
        """Test that ChIP broad uses PE peak calling with broad mode."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "calling_peaks_macs2_pe" in output, \
            "PE peak calling rule should be in ChIP broad workflow"
        assert "--broad" in output, \
            "ChIP broad should use --broad flag in MACS2"


class TestChIPNarrowWorkflow:
    """Test ChIP-seq narrow-peak workflow dry-run (SE samples)."""

    def test_chip_narrow_dryrun_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that dry-run succeeds for ChIP narrow coverage bigwig."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_NARROW_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_chip_narrow_includes_expected_rules(self, snakemake_available, repo_root, test_config):
        """Test that ChIP narrow workflow includes SE mapping rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_NARROW_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "bowtie2_map_se",
            "filter_chip_se",
            "make_coverage_chip",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_chip_narrow_excludes_pe_mapping(self, snakemake_available, repo_root, test_config):
        """Test that SE ChIP narrow samples do not trigger PE mapping rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_NARROW_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "bowtie2_map_pe" not in output, \
            "PE mapping rule should not appear for SE ChIP narrow sample"


class TestRNAseqWorkflow:
    """Test RNA-seq workflow dry-run (SE samples)."""

    def test_rnaseq_dryrun_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that dry-run succeeds for RNA-seq final BAM."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, RNA_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rnaseq_includes_expected_rules(self, snakemake_available, repo_root, test_config):
        """Test that RNA-seq workflow includes SE STAR mapping and filtering."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, RNA_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "STAR_map_se",
            "filter_rna_se",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_rnaseq_deg_target(self, snakemake_available, repo_root, test_config):
        """Test that DEG calling dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, RNA_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "call_all_DEGs" in output, \
            "RNA checkpoint should include DEG calling rule"


class TestSmallRNAWorkflow:
    """Test small RNA-seq workflow dry-run (SE samples)."""

    def test_srna_dryrun_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that dry-run succeeds for ShortStack results."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, SRNA_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_srna_includes_expected_rules(self, snakemake_available, repo_root, test_config):
        """Test that sRNA workflow includes ShortStack and cluster rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, SRNA_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        expected_rules = [
            "shortstack_map",
            "merging_srna_replicates",
        ]

        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_srna_structural_rna_depletion(self, snakemake_available, repo_root, test_config):
        """Test that structural RNA depletion rules are present."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, SRNA_TARGET, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "filter_structural_rna" in output, \
            "Structural RNA depletion should be present (structural_rna_depletion: true)"


class TestControlLinking:
    """Test that ChIP samples correctly resolve their Input controls."""

    def test_chip_broad_resolves_control(self, snakemake_available, repo_root, test_config):
        """Test that ChIP broad workflow references WT Input control in peak calling."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "WT_cell_Input_rep1" in output, \
            "ChIP broad workflow should reference WT_cell_Input_rep1 as control"

    def test_chip_narrow_resolves_control(self, snakemake_available, repo_root, test_config):
        """Test that ChIP narrow workflow references veg Input control in peak calling."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "veg_cell_Input_rep1" in output, \
            "ChIP workflow should reference veg_cell_Input_rep1 as control"

    def test_input_control_processed(self, snakemake_available, repo_root, test_config):
        """Test that Input control samples can be processed independently."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, INPUT_BROAD_TARGET)

        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


class TestDAGStructure:
    """Test DAG generation and structure."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_config):
        """Test that DAG can be generated for final checkpoint."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_config, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_all_assay_rules(self, snakemake_available, repo_root, test_config):
        """Test that DAG contains rules from all four assay types."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_config, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        dag_output = result.stdout

        # Rules from each assay type
        assay_rules = {
            "ChIP": "bowtie2_map",
            "RNA": "STAR_map",
            "sRNA": "shortstack_map",
        }

        for assay, rule in assay_rules.items():
            assert rule in dag_output, f"Rule '{rule}' for {assay} not found in DAG"

    def test_dag_has_dependencies(self, snakemake_available, repo_root, test_config):
        """Test that DAG has rule dependencies (edges)."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dag(repo_root, test_config, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"

        assert "->" in result.stdout, "DAG should contain rule dependencies"


class TestReplicateHandling:
    """Test handling of multiple replicates."""

    def test_chip_broad_replicates_resolve(self, snakemake_available, repo_root, test_config):
        """Test that both ChIP broad H3K9me2 replicates can be resolved."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        for target in [CHIP_BROAD_TARGET, CHIP_BROAD_TARGET_REP2]:
            result = run_snakemake_dryrun(repo_root, test_config, target)
            assert result.returncode == 0, f"Replicate failed for {target}: {result.stderr}"

    def test_idr_analysis_present(self, snakemake_available, repo_root, test_config):
        """Test that IDR analysis rule is present for replicated ChIP samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_CHECKPOINT, ["--printshellcmds"])

        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"

        output = result.stdout + result.stderr

        assert "idr_analysis_replicates" in output, \
            "IDR analysis should be present for replicated ChIP samples"

    def test_rna_replicates_resolve(self, snakemake_available, repo_root, test_config):
        """Test that multiple RNA-seq replicates can be resolved."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        targets = [
            "results/RNA/mapped/final__WT_cell_RNA_rep1.bam",
            "results/RNA/mapped/final__WT_cell_RNA_rep2.bam",
        ]

        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_config, target)
            assert result.returncode == 0, f"Replicate failed for {target}: {result.stderr}"


class TestEnvCheckpoints:
    """Test environment-level checkpoint targets (validates full DAG resolution)."""

    def test_chip_checkpoint(self, snakemake_available, repo_root, test_config):
        """Test that ChIP analysis checkpoint dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, CHIP_CHECKPOINT)

        assert result.returncode == 0, f"ChIP checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rna_checkpoint(self, snakemake_available, repo_root, test_config):
        """Test that RNA analysis checkpoint dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, RNA_CHECKPOINT)

        assert result.returncode == 0, f"RNA checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_srna_checkpoint(self, snakemake_available, repo_root, test_config):
        """Test that sRNA analysis checkpoint dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, SRNA_CHECKPOINT)

        assert result.returncode == 0, f"sRNA checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_final_checkpoint(self, snakemake_available, repo_root, test_config):
        """Test that final combined checkpoint dry-run succeeds.

        This is the most comprehensive test — it forces Snakemake to resolve
        the entire multi-assay DAG (~257 steps) in dry-run mode.
        """
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, FINAL_CHECKPOINT)

        assert result.returncode == 0, f"Final checkpoint failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


class TestMapOnly:
    """Test intermediate target rules."""

    def test_map_only_target(self, snakemake_available, repo_root, test_config):
        """Test that map_only target dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, "map_only")

        assert result.returncode == 0, f"map_only failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_coverage_chip_target(self, snakemake_available, repo_root, test_config):
        """Test that coverage_chip target dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, "coverage_chip")

        assert result.returncode == 0, f"coverage_chip failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
