"""
Integration tests for RNA-seq and RAMPAGE workflows using Snakemake dry-run.

These tests verify that the Snakemake DAG can be correctly built and wildcards
resolved for RNAseq and RAMPAGE assay types without actually executing the rules.

Test samples (from test_samples_RNA.tsv):
  RNAseq:
    WT_RNA_rep1       RNAseq   PE   genotype:WT
    WT_RNA_rep2       RNAseq   PE   genotype:WT
    mutant_RNA_rep1   RNAseq   PE   genotype:mutant
    mutant_RNA_rep2   RNAseq   PE   genotype:mutant
  RAMPAGE:
    WT_RAMPAGE_rep1   RAMPAGE  PE   genotype:WT   (Control: WT_RNA_rep1)
    WT_RAMPAGE_rep2   RAMPAGE  PE   genotype:WT   (Control: WT_RNA_rep2)
"""

import pytest
from pathlib import Path

from tests.integration.conftest import (
    load_output_dir, run_snakemake_dryrun, run_snakemake_dag,
)


# Mark all tests in this module as integration tests
pytestmark = pytest.mark.integration

_OUTPUT_DIR = load_output_dir("test_options_RNA.yaml")


# ---------------------------------------------------------------------------
# Target path constants
# ---------------------------------------------------------------------------

# Per-replicate stranded bigwig targets (RNAseq = reverse stranded)
RNASEQ_WT_REP1_PLUS = f"{_OUTPUT_DIR}/RNA/tracks/WT_RNA_rep1__plus.bw"
RNASEQ_WT_REP1_MINUS = f"{_OUTPUT_DIR}/RNA/tracks/WT_RNA_rep1__minus.bw"
RNASEQ_WT_REP2_PLUS = f"{_OUTPUT_DIR}/RNA/tracks/WT_RNA_rep2__plus.bw"
RNASEQ_MUTANT_REP1_PLUS = f"{_OUTPUT_DIR}/RNA/tracks/mutant_RNA_rep1__plus.bw"

# Per-replicate stranded bigwig targets (RAMPAGE = forward stranded)
RAMPAGE_WT_REP1_PLUS = f"{_OUTPUT_DIR}/RNA/tracks/WT_RAMPAGE_rep1__plus.bw"
RAMPAGE_WT_REP1_MINUS = f"{_OUTPUT_DIR}/RNA/tracks/WT_RAMPAGE_rep1__minus.bw"

# Analysis-level names (Assay__levels_label__Genome)
RNASEQ_WT_ANALYSIS = "RNAseq__WT__test_genome"
RNASEQ_MUTANT_ANALYSIS = "RNAseq__mutant__test_genome"
RAMPAGE_WT_ANALYSIS = "RAMPAGE__WT__test_genome"

# Analysis-level merged bigwig targets (mutant, because WT RNAseq samples are
# treated as controls by RAMPAGE and excluded from analysis_samples)
RNASEQ_MUTANT_MERGED_PLUS = f"{_OUTPUT_DIR}/RNA/tracks/{RNASEQ_MUTANT_ANALYSIS}__plus.bw"
RNASEQ_MUTANT_MERGED_MINUS = f"{_OUTPUT_DIR}/RNA/tracks/{RNASEQ_MUTANT_ANALYSIS}__minus.bw"

# DEG targets
DEG_CHECKPOINT = f"{_OUTPUT_DIR}/RNA/chkpts/calling_DEGs__test_RNA__test_genome.done"
DEG_RPKM = f"{_OUTPUT_DIR}/RNA/DEG/genes_rpkm__test_RNA__test_genome.txt"
DEG_EXPRESSION_PLOT = f"{_OUTPUT_DIR}/RNA/plots/plot_expression__test_RNA__test_genome__unique_DEGs.pdf"

# RAMPAGE TSS targets
TSS_FINAL_REP1 = f"{_OUTPUT_DIR}/RNA/TSS/TSS__final__WT_RAMPAGE_rep1_peaks.narrowPeak"
TSS_FINAL_REP2 = f"{_OUTPUT_DIR}/RNA/TSS/TSS__final__WT_RAMPAGE_rep2_peaks.narrowPeak"
# Note: merged RAMPAGE TSS is not testable here because the WT RNAseq samples
# are treated as controls and excluded from analysis_samples, so assign_rna_input
# cannot resolve the merged control BAM for RAMPAGE__WT__test_genome.


@pytest.fixture(scope="module")
def test_options(repo_root):
    """Get the path to the test options file."""
    config_path = repo_root / "tests" / "integration" / "data" / "test_options_RNA.yaml"
    assert config_path.exists(), f"Test options file not found at {config_path}"
    return str(config_path)


# ===========================================================================
# TestBasic
# ===========================================================================

class TestBasic:
    """Basic sanity checks for the RNA test configuration."""

    def test_snakemake_installed(self, snakemake_available):
        """Test that snakemake is available."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed or not in PATH")

    def test_options_file_exists(self, test_options):
        """Test that the test options file exists."""
        assert Path(test_options).exists()

    def test_sample_file_exists(self, repo_root):
        """Test that the test sample file exists."""
        sample_file = repo_root / "tests" / "integration" / "data" / "test_samples_RNA.tsv"
        assert sample_file.exists()


# ===========================================================================
# TestRNAseqWorkflow
# ===========================================================================

class TestRNAseqWorkflow:
    """Test RNAseq workflow dry-run."""

    def test_rnaseq_pe_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that dry-run succeeds for RNAseq PE stranded bigwig."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, RNASEQ_WT_REP1_PLUS)
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rnaseq_pe_includes_star_rules(self, snakemake_available, repo_root, test_options):
        """Test that RNAseq PE workflow includes STAR_map_pe and filter_rna_pe."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, RNASEQ_WT_REP1_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "STAR_map_pe" in output, "RNAseq PE should use STAR_map_pe rule"
        assert "filter_rna_pe" in output, "RNAseq PE should use filter_rna_pe rule"

    def test_rnaseq_pe_includes_stranded_bigwig_rule(self, snakemake_available, repo_root, test_options):
        """Test that RNAseq PE workflow includes make_rna_stranded_bigwigs."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, RNASEQ_WT_REP1_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "make_rna_stranded_bigwigs" in output, "RNAseq should use make_rna_stranded_bigwigs rule"


# ===========================================================================
# TestStrandedBigwigs
# ===========================================================================

class TestStrandedBigwigs:
    """Test that stranded bigwig targets resolve for both plus and minus strands."""

    def test_plus_strand_resolves(self, snakemake_available, repo_root, test_options):
        """Test that plus-strand bigwig resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, RNASEQ_WT_REP1_PLUS)
        assert result.returncode == 0, f"Plus strand failed: {result.stderr}"

    def test_minus_strand_resolves(self, snakemake_available, repo_root, test_options):
        """Test that minus-strand bigwig resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, RNASEQ_WT_REP1_MINUS)
        assert result.returncode == 0, f"Minus strand failed: {result.stderr}"

    def test_rampage_stranded_bigwigs_resolve(self, snakemake_available, repo_root, test_options):
        """Test that RAMPAGE stranded bigwigs resolve."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        for target in (RAMPAGE_WT_REP1_PLUS, RAMPAGE_WT_REP1_MINUS):
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"RAMPAGE strand target failed for {target}: {result.stderr}"


# ===========================================================================
# TestReplicateMerging
# ===========================================================================

class TestReplicateMerging:
    """Test replicate merging rules for RNAseq with multiple replicates."""

    def test_rnaseq_two_reps_trigger_merging(self, snakemake_available, repo_root, test_options):
        """RNAseq mutant has 2 reps -- merged replicate rules should appear.

        Note: WT RNAseq samples are identified as controls (referenced by
        RAMPAGE samples) and excluded from analysis_samples, so we test
        merging with the mutant samples instead.
        """
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, RNASEQ_MUTANT_MERGED_PLUS, ["--printshellcmds"])
        assert result.returncode == 0, f"Dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        output = result.stdout + result.stderr
        assert "merging_rna_replicates" in output, \
            "RNAseq with 2 reps should trigger merging_rna_replicates rule"

    def test_merged_bigwig_resolves(self, snakemake_available, repo_root, test_options):
        """Test that merged analysis-level bigwig target resolves."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        for target in (RNASEQ_MUTANT_MERGED_PLUS, RNASEQ_MUTANT_MERGED_MINUS):
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Merged bigwig failed for {target}: {result.stderr}"


# ===========================================================================
# TestDEGAnalysis
# ===========================================================================

class TestDEGAnalysis:
    """Test DEG (differential expression) analysis dry-run."""

    def test_deg_checkpoint_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DEG calling dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DEG_CHECKPOINT)
        assert result.returncode == 0, f"DEG dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_deg_includes_r_script(self, snakemake_available, repo_root, test_options):
        """Test that DEG workflow includes the R_call_DEGs.R script."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DEG_CHECKPOINT, ["--printshellcmds"])
        assert result.returncode == 0, f"DEG dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "R_call_DEGs.R" in output, "DEG workflow should include R_call_DEGs.R"

    def test_rpkm_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that RPKM gathering dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DEG_RPKM)
        assert result.returncode == 0, f"RPKM dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_expression_plot_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that expression plot dry-run succeeds."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, DEG_EXPRESSION_PLOT)
        assert result.returncode == 0, f"Expression plot dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"


# ===========================================================================
# TestRAMPAGEWorkflow
# ===========================================================================

class TestRAMPAGEWorkflow:
    """Test RAMPAGE TSS calling workflow dry-run."""

    def test_rampage_tss_dryrun_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that RAMPAGE TSS dry-run succeeds for a per-replicate target."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, TSS_FINAL_REP1)
        assert result.returncode == 0, f"TSS dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rampage_uses_rnaseq_control(self, snakemake_available, repo_root, test_options):
        """Test that RAMPAGE TSS calling uses RNAseq as control input."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, TSS_FINAL_REP1, ["--printshellcmds"])
        assert result.returncode == 0, f"TSS dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "call_rampage_TSS" in output, "RAMPAGE workflow should include call_rampage_TSS rule"
        # The control BAM should reference the RNAseq sample
        assert "WT_RNA_rep1" in output, "RAMPAGE TSS should use WT_RNA_rep1 as control"

    def test_rampage_tss_rep2_resolves(self, snakemake_available, repo_root, test_options):
        """Test that RAMPAGE TSS resolves for rep2."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, TSS_FINAL_REP2)
        assert result.returncode == 0, f"TSS rep2 dry-run failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"

    def test_rampage_includes_macs2(self, snakemake_available, repo_root, test_options):
        """Test that RAMPAGE TSS calling uses macs2."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dryrun(repo_root, test_options, TSS_FINAL_REP1, ["--printshellcmds"])
        assert result.returncode == 0, f"TSS dry-run failed: {result.stderr}"
        output = result.stdout + result.stderr
        assert "macs2" in output, "RAMPAGE TSS calling should use macs2"


# ===========================================================================
# TestWildcardResolution
# ===========================================================================

class TestWildcardResolution:
    """Test wildcard resolution for all RNA sample types."""

    def test_rnaseq_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for all RNAseq samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/RNA/tracks/WT_RNA_rep1__plus.bw",
            f"{_OUTPUT_DIR}/RNA/tracks/WT_RNA_rep2__minus.bw",
            f"{_OUTPUT_DIR}/RNA/tracks/mutant_RNA_rep1__plus.bw",
            f"{_OUTPUT_DIR}/RNA/tracks/mutant_RNA_rep2__minus.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"

    def test_rampage_wildcards_resolve(self, snakemake_available, repo_root, test_options):
        """Test that wildcards resolve for RAMPAGE samples."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        targets = [
            f"{_OUTPUT_DIR}/RNA/tracks/WT_RAMPAGE_rep1__plus.bw",
            f"{_OUTPUT_DIR}/RNA/tracks/WT_RAMPAGE_rep2__minus.bw",
        ]
        for target in targets:
            result = run_snakemake_dryrun(repo_root, test_options, target)
            assert result.returncode == 0, f"Wildcard resolution failed for {target}: {result.stderr}"


# ===========================================================================
# TestDAGStructure
# ===========================================================================

class TestDAGStructure:
    """Test DAG generation and structure for RNA workflows."""

    def test_dag_generation_succeeds(self, snakemake_available, repo_root, test_options):
        """Test that DAG can be generated."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, RNASEQ_WT_REP1_PLUS)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        assert len(result.stdout) > 0, "DAG output is empty"

    def test_dag_contains_star_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG for an RNAseq target contains STAR-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, RNASEQ_WT_REP1_PLUS)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        for rule in ("make_STAR_indices", "STAR_map_pe", "make_rna_stranded_bigwigs"):
            assert rule in dag_output, f"Rule '{rule}' not found in RNA DAG"

    def test_dag_contains_tss_rules(self, snakemake_available, repo_root, test_options):
        """Test that DAG for a RAMPAGE TSS target contains TSS-specific rules."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, TSS_FINAL_REP1)
        assert result.returncode == 0, f"DAG generation failed: {result.stderr}"
        dag_output = result.stdout
        assert "call_rampage_TSS" in dag_output, "Rule 'call_rampage_TSS' not found in RAMPAGE DAG"

    def test_dag_rule_dependencies(self, snakemake_available, repo_root, test_options):
        """Test that DAG shows correct rule dependencies."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        result = run_snakemake_dag(repo_root, test_options, RNASEQ_WT_REP1_PLUS)
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
        target = f"{_OUTPUT_DIR}/RNA/tracks/nonexistent_sample__plus.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for non-existent sample"

    def test_invalid_strand_suffix(self, snakemake_available, repo_root, test_options):
        """Test that invalid strand suffix fails."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")
        target = f"{_OUTPUT_DIR}/RNA/tracks/WT_RNA_rep1__INVALID.bw"
        result = run_snakemake_dryrun(repo_root, test_options, target)
        assert result.returncode != 0, "Should fail for invalid strand suffix"
