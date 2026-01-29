"""
Unit tests for helper functions in workflow/rules/mC.smk

Tests the dmC-related helper functions:
- is_dmc_sample: checks if sample uses dmC (direct methylation) workflow
- get_dmc_input_type: determines input type (modBAM or bedMethyl)
- parameters_for_mc: returns the correct parameter set for methylation calling
"""

import pytest


def parse_sample_name(sample_name):
    """
    Mock implementation of parse_sample_name from Snakefile.

    Parses sample names in format:
    {data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}
    """
    parts = sample_name.split("__")
    if len(parts) != 6:
        raise ValueError(f"Invalid sample name format: {sample_name}")

    data_type, line, tissue, sample_type, rep, ref_genome = parts
    parsed = {
        "data_type": data_type,
        "line": line,
        "tissue": tissue,
        "sample_type": sample_type,
        "replicate": rep,
        "ref_genome": ref_genome
    }

    # Extract ChIP Input group from data_type
    if data_type.startswith("ChIP_"):
        chip_parts = data_type.split("_", 1)
        if len(chip_parts) == 2:
            parsed["group_label"] = chip_parts[1]

    # Extract TF name from data_type
    if data_type.startswith("TF_"):
        tf_parts = data_type.split("_", 1)
        if len(tf_parts) == 2:
            parsed["tf_name"] = tf_parts[1]

    return parsed


def parameters_for_mc(sample_name):
    """
    Determine methylation calling parameters based on sample type.

    Extracted from mC.smk for testing.
    """
    temp = parse_sample_name(sample_name)['sample_type']
    options = {"WGBS", "Pico", "EMseq", "dmC", "bedMethyl"}
    return temp if temp in options else "default"


def is_dmc_sample(sample_name):
    """
    Check if a sample uses dmC (direct methylation) workflow.

    Extracted from mC.smk for testing.
    """
    return parse_sample_name(sample_name)['sample_type'] in ["dmC", "bedMethyl"]


def get_dmc_input_type(sample_name):
    """
    Return the input type for dmC samples: 'bedMethyl' or 'modBAM'.

    Extracted from mC.smk for testing.
    """
    sample_type = parse_sample_name(sample_name)['sample_type']
    return "bedMethyl" if sample_type == "bedMethyl" else "modBAM"


class TestParseSampleName:
    """Tests for the parse_sample_name helper function."""

    def test_parse_standard_sample_name(self, dmc_sample_names):
        """Test parsing a standard sample name."""
        result = parse_sample_name(dmc_sample_names["dmc_modbam"])

        assert result["data_type"] == "mC"
        assert result["line"] == "Col0"
        assert result["tissue"] == "leaf"
        assert result["sample_type"] == "dmC"
        assert result["replicate"] == "rep1"
        assert result["ref_genome"] == "ColCEN"

    def test_parse_bedmethyl_sample_name(self, dmc_sample_names):
        """Test parsing a bedMethyl sample name."""
        result = parse_sample_name(dmc_sample_names["dmc_bedmethyl"])

        assert result["sample_type"] == "bedMethyl"

    def test_parse_chip_sample_with_group(self):
        """Test parsing ChIP sample with group label."""
        sample_name = "ChIP_H3K4me3__Col0__leaf__Input__rep1__ColCEN"
        result = parse_sample_name(sample_name)

        assert result["data_type"] == "ChIP_H3K4me3"
        assert result["group_label"] == "H3K4me3"

    def test_parse_tf_sample(self):
        """Test parsing TF sample with TF name."""
        sample_name = "TF_TB1__Col0__leaf__Input__rep1__ColCEN"
        result = parse_sample_name(sample_name)

        assert result["data_type"] == "TF_TB1"
        assert result["tf_name"] == "TB1"

    def test_parse_invalid_sample_name(self):
        """Test parsing fails with invalid format."""
        with pytest.raises(ValueError):
            parse_sample_name("invalid__sample__name")


class TestIsDmcSample:
    """Tests for the is_dmc_sample helper function."""

    def test_dmc_modbam_sample_is_dmc(self, dmc_sample_names):
        """Test that dmC modBAM sample is identified as dmC."""
        assert is_dmc_sample(dmc_sample_names["dmc_modbam"]) is True

    def test_bedmethyl_sample_is_dmc(self, dmc_sample_names):
        """Test that bedMethyl sample is identified as dmC."""
        assert is_dmc_sample(dmc_sample_names["dmc_bedmethyl"]) is True

    def test_wgbs_sample_is_not_dmc(self, dmc_sample_names):
        """Test that WGBS sample is not identified as dmC."""
        assert is_dmc_sample(dmc_sample_names["bismark_wgbs"]) is False

    def test_pico_sample_is_not_dmc(self, dmc_sample_names):
        """Test that Pico sample is not identified as dmC."""
        assert is_dmc_sample(dmc_sample_names["bismark_pico"]) is False

    def test_emseq_sample_is_not_dmc(self, dmc_sample_names):
        """Test that EMseq sample is not identified as dmC."""
        assert is_dmc_sample(dmc_sample_names["bismark_emseq"]) is False

    def test_default_sample_is_not_dmc(self, dmc_sample_names):
        """Test that default sample is not identified as dmC."""
        assert is_dmc_sample(dmc_sample_names["bismark_default"]) is False


class TestGetDmcInputType:
    """Tests for the get_dmc_input_type helper function."""

    def test_bedmethyl_sample_returns_bedmethyl(self, dmc_sample_names):
        """Test that bedMethyl sample type returns 'bedMethyl'."""
        result = get_dmc_input_type(dmc_sample_names["dmc_bedmethyl"])
        assert result == "bedMethyl"

    def test_dmc_sample_returns_modbam(self, dmc_sample_names):
        """Test that dmC sample type returns 'modBAM'."""
        result = get_dmc_input_type(dmc_sample_names["dmc_modbam"])
        assert result == "modBAM"

    def test_wgbs_sample_returns_modbam(self, dmc_sample_names):
        """Test that non-dmC sample types default to 'modBAM'."""
        # Note: This function should only be called on dmC samples,
        # but testing edge case behavior
        result = get_dmc_input_type(dmc_sample_names["bismark_wgbs"])
        assert result == "modBAM"


class TestParametersForMc:
    """Tests for the parameters_for_mc helper function."""

    def test_dmc_sample_returns_dmc(self, dmc_sample_names):
        """Test that dmC sample returns 'dmC' parameter set."""
        result = parameters_for_mc(dmc_sample_names["dmc_modbam"])
        assert result == "dmC"

    def test_bedmethyl_sample_returns_bedmethyl(self, dmc_sample_names):
        """Test that bedMethyl sample returns 'bedMethyl' parameter set."""
        result = parameters_for_mc(dmc_sample_names["dmc_bedmethyl"])
        assert result == "bedMethyl"

    def test_wgbs_sample_returns_wgbs(self, dmc_sample_names):
        """Test that WGBS sample returns 'WGBS' parameter set."""
        result = parameters_for_mc(dmc_sample_names["bismark_wgbs"])
        assert result == "WGBS"

    def test_pico_sample_returns_pico(self, dmc_sample_names):
        """Test that Pico sample returns 'Pico' parameter set."""
        result = parameters_for_mc(dmc_sample_names["bismark_pico"])
        assert result == "Pico"

    def test_emseq_sample_returns_emseq(self, dmc_sample_names):
        """Test that EMseq sample returns 'EMseq' parameter set."""
        result = parameters_for_mc(dmc_sample_names["bismark_emseq"])
        assert result == "EMseq"

    def test_unknown_sample_type_returns_default(self, dmc_sample_names):
        """Test that unknown sample type returns 'default' parameter set."""
        result = parameters_for_mc(dmc_sample_names["bismark_default"])
        assert result == "default"

    def test_other_sample_types_return_default(self):
        """Test that other sample types return 'default'."""
        sample_name = "mC__Col0__leaf__custom__rep1__ColCEN"
        result = parameters_for_mc(sample_name)
        assert result == "default"


class TestEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_malformed_sample_name_raises_error(self):
        """Test that malformed sample names raise appropriate errors."""
        with pytest.raises(ValueError):
            is_dmc_sample("malformed_sample_name")

    def test_empty_sample_name_raises_error(self):
        """Test that empty sample name raises error."""
        with pytest.raises(ValueError):
            is_dmc_sample("")

    def test_sample_name_with_underscores_in_fields(self):
        """Test sample names with underscores in individual fields."""
        # This should work - underscores are allowed within fields,
        # double underscores are field separators
        sample_name = "mC__Col_0__leaf_tissue__dmC__rep_1__Col_CEN"
        result = parse_sample_name(sample_name)

        assert result["line"] == "Col_0"
        assert result["tissue"] == "leaf_tissue"
        assert result["replicate"] == "rep_1"
        assert result["ref_genome"] == "Col_CEN"

    def test_case_sensitivity(self):
        """Test that sample type matching is case-sensitive."""
        # Should not be recognized as dmC (uppercase)
        sample_name = "mC__Col0__leaf__DMC__rep1__ColCEN"
        assert is_dmc_sample(sample_name) is False

        # Lowercase should work
        sample_name = "mC__Col0__leaf__dmC__rep1__ColCEN"
        assert is_dmc_sample(sample_name) is True


class TestIntegrationScenarios:
    """Integration tests for realistic usage scenarios."""

    def test_complete_workflow_dmc_modbam(self):
        """Test complete workflow for dmC modBAM sample."""
        sample_name = "mC__Col0__leaf__dmC__rep1__ColCEN"

        # Check if it's a dmC sample
        assert is_dmc_sample(sample_name) is True

        # Get the input type
        input_type = get_dmc_input_type(sample_name)
        assert input_type == "modBAM"

        # Get the parameter set
        params = parameters_for_mc(sample_name)
        assert params == "dmC"

    def test_complete_workflow_dmc_bedmethyl(self):
        """Test complete workflow for dmC bedMethyl sample."""
        sample_name = "mC__B73__root__bedMethyl__rep2__B73_v5"

        # Check if it's a dmC sample
        assert is_dmc_sample(sample_name) is True

        # Get the input type
        input_type = get_dmc_input_type(sample_name)
        assert input_type == "bedMethyl"

        # Get the parameter set
        params = parameters_for_mc(sample_name)
        assert params == "bedMethyl"

    def test_complete_workflow_bismark(self):
        """Test complete workflow for Bismark sample."""
        sample_name = "mC__Col0__leaf__WGBS__rep1__ColCEN"

        # Check if it's a dmC sample
        assert is_dmc_sample(sample_name) is False

        # Get the parameter set (should not call get_dmc_input_type for Bismark)
        params = parameters_for_mc(sample_name)
        assert params == "WGBS"

    def test_multiple_replicates(self):
        """Test handling multiple replicates of dmC samples."""
        samples = [
            "mC__Col0__leaf__dmC__rep1__ColCEN",
            "mC__Col0__leaf__dmC__rep2__ColCEN",
            "mC__Col0__leaf__dmC__rep3__ColCEN",
        ]

        for sample in samples:
            assert is_dmc_sample(sample) is True
            assert get_dmc_input_type(sample) == "modBAM"
            assert parameters_for_mc(sample) == "dmC"

    def test_mixed_sample_types(self):
        """Test handling a mix of dmC and Bismark samples."""
        samples = {
            "mC__Col0__leaf__dmC__rep1__ColCEN": ("dmC", True),
            "mC__Col0__leaf__WGBS__rep1__ColCEN": ("WGBS", False),
            "mC__Col0__leaf__bedMethyl__rep1__ColCEN": ("bedMethyl", True),
            "mC__Col0__leaf__Pico__rep1__ColCEN": ("Pico", False),
        }

        for sample_name, (expected_param, expected_is_dmc) in samples.items():
            assert is_dmc_sample(sample_name) == expected_is_dmc
            assert parameters_for_mc(sample_name) == expected_param
