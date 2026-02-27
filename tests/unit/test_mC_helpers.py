"""
Unit tests for mC-related helper functions.

Tests the helper functions used by workflow/rules/mC.smk:
- is_dmc_sample: checks if sample uses dmC (direct methylation) workflow
- parameters_for_mc: returns the correct parameter set for methylation calling

These functions look up sample metadata via the samples DataFrame (new format),
so tests construct a mock DataFrame with the current column schema.
"""

import pytest
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers under test – extracted from mC.smk logic
# These mirror the actual implementations which use get_sample_info_from_name()
# ---------------------------------------------------------------------------

def _make_samples_df(rows):
    """Build a mock samples DataFrame from a list of dicts."""
    df = pd.DataFrame(rows)
    df["sample_name"] = df["Sample_ID"]
    return df


def _get_sample_info(sample_name, df, field):
    """Look up a field value by sample_name (mirrors get_sample_info_from_name)."""
    match = df.loc[df["sample_name"] == sample_name]
    if match.empty:
        return None
    return match[field].iloc[0]


def is_dmc_sample(sample_name, df):
    """Check if a sample uses direct methylation (dmC) workflow."""
    return _get_sample_info(sample_name, df, "Assay") == "dmC"


def parameters_for_mc(sample_name, df):
    """Determine methylation calling parameters based on Assay."""
    assay = _get_sample_info(sample_name, df, "Assay")
    options = {"WGBS", "Pico", "EMseq", "dmC"}
    return assay if assay in options else "default"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mc_samples():
    """Return a samples DataFrame with various mC sample types."""
    rows = [
        {"Sample_ID": "WT_leaf_dmC_rep1", "Assay": "dmC", "Genome": "ColCEN",
         "Levels": "genotype:Col0,tissue:leaf", "Replicate_ID": "rep1",
         "Read_files": "mock.bam", "Read_layout": "SE", "IP_target": "", "Control": ""},
        {"Sample_ID": "WT_leaf_dmC_rep2", "Assay": "dmC", "Genome": "ColCEN",
         "Levels": "genotype:Col0,tissue:leaf", "Replicate_ID": "rep2",
         "Read_files": "mock2.bam", "Read_layout": "SE", "IP_target": "", "Control": ""},
        {"Sample_ID": "WT_root_WGBS_rep1", "Assay": "WGBS", "Genome": "ColCEN",
         "Levels": "genotype:Col0,tissue:root", "Replicate_ID": "rep1",
         "Read_files": "SRR12345", "Read_layout": "SE", "IP_target": "", "Control": ""},
        {"Sample_ID": "WT_leaf_Pico_rep1", "Assay": "Pico", "Genome": "ColCEN",
         "Levels": "genotype:Col0,tissue:leaf", "Replicate_ID": "rep1",
         "Read_files": "SRR12346", "Read_layout": "SE", "IP_target": "", "Control": ""},
        {"Sample_ID": "WT_leaf_EMseq_rep1", "Assay": "EMseq", "Genome": "ColCEN",
         "Levels": "genotype:Col0,tissue:leaf", "Replicate_ID": "rep1",
         "Read_files": "SRR12347", "Read_layout": "PE", "IP_target": "", "Control": ""},
        {"Sample_ID": "WT_leaf_mC_rep1", "Assay": "mC", "Genome": "ColCEN",
         "Levels": "genotype:Col0,tissue:leaf", "Replicate_ID": "rep1",
         "Read_files": "SRR12348", "Read_layout": "SE", "IP_target": "", "Control": ""},
    ]
    return _make_samples_df(rows)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestIsDmcSample:
    """Tests for the is_dmc_sample helper function."""

    def test_dmc_sample_is_dmc(self, mc_samples):
        assert is_dmc_sample("WT_leaf_dmC_rep1", mc_samples) is True

    def test_dmc_rep2_is_dmc(self, mc_samples):
        assert is_dmc_sample("WT_leaf_dmC_rep2", mc_samples) is True

    def test_wgbs_sample_is_not_dmc(self, mc_samples):
        assert is_dmc_sample("WT_root_WGBS_rep1", mc_samples) is False

    def test_pico_sample_is_not_dmc(self, mc_samples):
        assert is_dmc_sample("WT_leaf_Pico_rep1", mc_samples) is False

    def test_emseq_sample_is_not_dmc(self, mc_samples):
        assert is_dmc_sample("WT_leaf_EMseq_rep1", mc_samples) is False

    def test_default_mc_sample_is_not_dmc(self, mc_samples):
        assert is_dmc_sample("WT_leaf_mC_rep1", mc_samples) is False


class TestParametersForMc:
    """Tests for the parameters_for_mc helper function."""

    def test_dmc_returns_dmc(self, mc_samples):
        assert parameters_for_mc("WT_leaf_dmC_rep1", mc_samples) == "dmC"

    def test_wgbs_returns_wgbs(self, mc_samples):
        assert parameters_for_mc("WT_root_WGBS_rep1", mc_samples) == "WGBS"

    def test_pico_returns_pico(self, mc_samples):
        assert parameters_for_mc("WT_leaf_Pico_rep1", mc_samples) == "Pico"

    def test_emseq_returns_emseq(self, mc_samples):
        assert parameters_for_mc("WT_leaf_EMseq_rep1", mc_samples) == "EMseq"

    def test_unknown_assay_returns_default(self, mc_samples):
        assert parameters_for_mc("WT_leaf_mC_rep1", mc_samples) == "default"


class TestIntegrationScenarios:
    """Integration tests for realistic usage scenarios."""

    def test_complete_workflow_dmc(self, mc_samples):
        """Test dmC sample identification and parameter selection."""
        sample = "WT_leaf_dmC_rep1"
        assert is_dmc_sample(sample, mc_samples) is True
        assert parameters_for_mc(sample, mc_samples) == "dmC"

    def test_complete_workflow_bismark(self, mc_samples):
        """Test bisulfite sample identification and parameter selection."""
        sample = "WT_root_WGBS_rep1"
        assert is_dmc_sample(sample, mc_samples) is False
        assert parameters_for_mc(sample, mc_samples) == "WGBS"

    def test_multiple_replicates(self, mc_samples):
        """Test handling multiple replicates of dmC samples."""
        for sample in ["WT_leaf_dmC_rep1", "WT_leaf_dmC_rep2"]:
            assert is_dmc_sample(sample, mc_samples) is True
            assert parameters_for_mc(sample, mc_samples) == "dmC"

    def test_mixed_sample_types(self, mc_samples):
        """Test handling a mix of dmC and Bismark samples."""
        expected = {
            "WT_leaf_dmC_rep1": ("dmC", True),
            "WT_root_WGBS_rep1": ("WGBS", False),
            "WT_leaf_Pico_rep1": ("Pico", False),
            "WT_leaf_EMseq_rep1": ("EMseq", False),
            "WT_leaf_mC_rep1": ("default", False),
        }
        for sample, (exp_param, exp_dmc) in expected.items():
            assert is_dmc_sample(sample, mc_samples) == exp_dmc
            assert parameters_for_mc(sample, mc_samples) == exp_param
