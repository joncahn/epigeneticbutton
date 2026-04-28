"""Unit tests for workflow/scripts/sample_sheet.py"""

import pytest
import pandas as pd
from collections import OrderedDict
import sys
import os

# Add the workflow/scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "scripts"))

from sample_sheet import (
    VALID_ASSAYS,
    ASSAY_TO_ENV,
    ASSAY_TO_PEAKTYPE,
    NEW_COLNAMES,
    parse_levels,
    levels_to_label,
    levels_to_factors,
    parse_read_files,
    get_seq_id_and_path,
    build_analysis_key,
    build_analysis_name,
    build_analysis_to_replicates,
    identify_control_samples,
    get_control_sample_id,
    get_analysis_samples,
    get_sample_field,
    get_env,
    get_peaktype,
    add_compat_columns,
    read_sample_sheet,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def pombe_df():
    """A minimal pombe-like sample sheet in new format."""
    data = {
        "Sample_ID": [
            "WT_H3K9me2_rep1", "WT_H3K9me2_rep2",
            "WT_WCE_rep1",
            "dcr1_H3K9me2_rep1",
            "WT_RNA_rep1", "WT_RNA_rep2",
            "dcr1_RNA_rep1",
            "WT_sRNA_rep1",
        ],
        "Assay": [
            "ChIP_broad", "ChIP_broad",
            "ChIP_broad",
            "ChIP_broad",
            "RNAseq", "RNAseq",
            "RNAseq",
            "sRNA",
        ],
        "Genome": ["Spombe"] * 8,
        "Levels": [
            "genotype:WT", "genotype:WT",
            "genotype:WT",
            "genotype:dcr1",
            "genotype:WT", "genotype:WT",
            "genotype:dcr1",
            "genotype:WT",
        ],
        "Replicate_ID": [
            "rep1", "rep2", "rep1", "rep1",
            "rep1", "rep2", "rep1", "rep1",
        ],
        "Read_files": [
            "SRR20678305", "SRR20678333",
            "SRR5445712",
            "SRR20678308",
            "SRR30889044", "SRR30889043",
            "SRR30889039",
            "SRR20678362",
        ],
        "Read_layout": [
            "PE", "PE", "PE", "PE",
            "SE", "SE", "SE", "SE",
        ],
        "IP_target": [
            "H3K9me2", "H3K9me2", "WCE",
            "H3K9me2",
            "", "", "", "",
        ],
        "Control": [
            "WT_WCE_rep1", "WT_WCE_rep1", "",
            "WT_WCE_rep1",
            "", "", "", "",
        ],
    }
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# parse_levels
# ---------------------------------------------------------------------------

class TestParseLevels:
    def test_basic(self):
        result = parse_levels("genotype:WT,tissue:root")
        assert result == OrderedDict([("genotype", "WT"), ("tissue", "root")])

    def test_single_factor(self):
        result = parse_levels("genotype:WT")
        assert result == OrderedDict([("genotype", "WT")])

    def test_three_factors(self):
        result = parse_levels("genotype:WT,tissue:root,time:T0")
        assert len(result) == 3
        assert result["time"] == "T0"

    def test_whitespace_tolerant(self):
        result = parse_levels("genotype: WT , tissue: root")
        assert result == OrderedDict([("genotype", "WT"), ("tissue", "root")])

    def test_invalid_no_colon(self):
        with pytest.raises(ValueError, match="expected 'factor:level'"):
            parse_levels("genotype_WT")


class TestLevelsToLabel:
    def test_basic(self):
        assert levels_to_label("genotype:WT,tissue:root") == "WT_root"

    def test_single(self):
        assert levels_to_label("genotype:WT") == "WT"

    def test_three(self):
        assert levels_to_label("genotype:WT,tissue:root,time:T0") == "WT_root_T0"


class TestLevelsToFactors:
    def test_basic(self):
        assert levels_to_factors("genotype:WT,tissue:root") == ["genotype", "tissue"]


# ---------------------------------------------------------------------------
# parse_read_files
# ---------------------------------------------------------------------------

class TestParseReadFiles:
    def test_single_sra(self):
        parts, is_sra = parse_read_files("SRR20678305", "PE")
        assert parts == ["SRR20678305"]
        assert is_sra is True

    def test_merged_sra(self):
        parts, is_sra = parse_read_files("SRR111+SRR222+SRR333", "SE")
        assert parts == ["SRR111", "SRR222", "SRR333"]
        assert is_sra is True

    def test_ddbj_accession(self):
        parts, is_sra = parse_read_files("DRR400324", "PE")
        assert parts == ["DRR400324"]
        assert is_sra is True

    def test_ena_accession(self):
        parts, is_sra = parse_read_files("ERR123456", "SE")
        assert parts == ["ERR123456"]
        assert is_sra is True

    def test_local_se(self):
        parts, is_sra = parse_read_files("/path/to/reads.fq.gz", "SE")
        assert parts == ["/path/to/reads.fq.gz"]
        assert is_sra is False

    def test_local_pe(self):
        parts, is_sra = parse_read_files("/path/r1.fq.gz,/path/r2.fq.gz", "PE")
        assert parts == ["/path/r1.fq.gz,/path/r2.fq.gz"]
        assert is_sra is False


class TestGetSeqIdAndPath:
    def test_sra(self):
        seq_id, fq_path = get_seq_id_and_path("SRR20678305", "PE")
        assert seq_id == "SRR20678305"
        assert fq_path == "SRA"

    def test_sra_merged(self):
        seq_id, fq_path = get_seq_id_and_path("SRR111+SRR222", "SE")
        assert seq_id == "SRR111,SRR222"
        assert fq_path == "SRA"

    def test_ddbj(self):
        seq_id, fq_path = get_seq_id_and_path("DRR400324", "PE")
        assert seq_id == "DRR400324"
        assert fq_path == "SRA"

    def test_bam(self):
        seq_id, fq_path = get_seq_id_and_path("/path/to/sample.bam", "SE")
        assert seq_id == "sample"
        assert fq_path == "/path/to/sample.bam"

    def test_bam_url(self):
        url = "https://lemna.org/data/shimada2024/SRR28453434.bam"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "SRR28453434"
        assert fq_path == url

    def test_bedmethyl_url(self):
        url = "https://example.com/data/sample.bed.gz"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "sample"
        assert fq_path == url

    def test_bedmethyl_url_with_query_params(self):
        url = "https://example.com/data/sample.bed.gz?token=abc123"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "sample"
        assert fq_path == url

    def test_fastq_url_se(self):
        url = "https://example.com/reads.fastq.gz"
        seq_id, fq_path = get_seq_id_and_path(url, "SE")
        assert seq_id == "URL"
        assert fq_path == url

    def test_fastq_url_pe(self):
        r1 = "https://example.com/r1.fq.gz"
        r2 = "https://example.com/r2.fq.gz"
        seq_id, fq_path = get_seq_id_and_path(f"{r1},{r2}", "PE")
        assert seq_id == "URL"
        assert fq_path == f"{r1},{r2}"


# ---------------------------------------------------------------------------
# Analysis key / name
# ---------------------------------------------------------------------------

class TestBuildAnalysisKey:
    def test_chip(self, pombe_df):
        row = pombe_df.iloc[0]
        key = build_analysis_key(row)
        assert key == ("ChIP_broad", "WT", "H3K9me2", "Spombe")

    def test_rnaseq(self, pombe_df):
        row = pombe_df.iloc[5]
        key = build_analysis_key(row)
        assert key == ("RNAseq", "WT", "", "Spombe")


class TestBuildAnalysisName:
    def test_chip(self, pombe_df):
        row = pombe_df.iloc[0]
        name = build_analysis_name(row)
        assert name == "ChIP_broad__WT__H3K9me2__Spombe"

    def test_rnaseq(self, pombe_df):
        row = pombe_df.iloc[5]
        name = build_analysis_name(row)
        assert name == "RNAseq__WT__Spombe"


# ---------------------------------------------------------------------------
# Analysis-to-replicates
# ---------------------------------------------------------------------------

class TestBuildAnalysisToReplicates:
    def test_basic(self, pombe_df):
        a2r = build_analysis_to_replicates(pombe_df)
        # WT H3K9me2 ChIP has 2 reps
        key = ("ChIP_broad", "WT", "H3K9me2", "Spombe")
        assert key in a2r
        assert a2r[key] == ["rep1", "rep2"]

    def test_controls_excluded(self, pombe_df):
        a2r = build_analysis_to_replicates(pombe_df)
        # Input samples should not appear as analysis keys
        for key in a2r:
            assert key[2] != "WCE"  # IP_target should never be a control type

    def test_single_rep(self, pombe_df):
        a2r = build_analysis_to_replicates(pombe_df)
        key = ("ChIP_broad", "dcr1", "H3K9me2", "Spombe")
        assert key in a2r
        assert a2r[key] == ["rep1"]


# ---------------------------------------------------------------------------
# Control handling
# ---------------------------------------------------------------------------

class TestIdentifyControlSamples:
    def test_basic(self, pombe_df):
        controls = identify_control_samples(pombe_df)
        assert controls == {"WT_WCE_rep1"}


class TestGetControlSampleId:
    def test_has_control(self, pombe_df):
        ctrl = get_control_sample_id("WT_H3K9me2_rep1", pombe_df)
        assert ctrl == "WT_WCE_rep1"

    def test_no_control(self, pombe_df):
        ctrl = get_control_sample_id("WT_RNA_rep1", pombe_df)
        assert ctrl is None

    def test_is_control(self, pombe_df):
        ctrl = get_control_sample_id("WT_WCE_rep1", pombe_df)
        assert ctrl is None


# ---------------------------------------------------------------------------
# Analysis samples (non-control)
# ---------------------------------------------------------------------------

class TestGetAnalysisSamples:
    def test_controls_excluded(self, pombe_df):
        analysis = get_analysis_samples(pombe_df)
        # Controls should not appear
        assert "WT_WCE_rep1" not in analysis["Sample_ID"].values
        # dcr1_WCE_rep1 removed — no dcr1 WCE in public data; dcr1 samples use WT_WCE_rep1

    def test_deduplicated(self, pombe_df):
        analysis = get_analysis_samples(pombe_df)
        # WT H3K9me2 has 2 reps but should appear once in analysis
        chip_wt = analysis[
            (analysis["Assay"] == "ChIP_broad") &
            (analysis["IP_target"] == "H3K9me2") &
            (analysis["Levels"].str.contains("WT"))
        ]
        assert len(chip_wt) == 1


# ---------------------------------------------------------------------------
# Field lookup
# ---------------------------------------------------------------------------

class TestGetSampleField:
    def test_basic(self, pombe_df):
        val = get_sample_field("WT_H3K9me2_rep1", pombe_df, "Assay")
        assert val == "ChIP_broad"

    def test_missing(self, pombe_df):
        val = get_sample_field("nonexistent", pombe_df, "Assay")
        assert val is None


# ---------------------------------------------------------------------------
# Environment and peak type
# ---------------------------------------------------------------------------

class TestGetEnv:
    def test_chip_broad(self):
        assert get_env("ChIP_broad") == "ChIP"

    def test_rnaseq(self):
        assert get_env("RNAseq") == "RNA"

    def test_unknown(self):
        with pytest.raises(ValueError, match="Unknown assay"):
            get_env("unknown_assay")


class TestGetPeaktype:
    def test_chip_broad(self):
        assert get_peaktype("ChIP_broad") == "broad"

    def test_chip_narrow(self):
        assert get_peaktype("ChIP_narrow") == "narrow"

    def test_atac(self):
        assert get_peaktype("ATAC") == "narrow"

    def test_cut_run_broad(self):
        assert get_peaktype("CUT_RUN_broad") == "broad"

    def test_cut_run_narrow(self):
        assert get_peaktype("CUT_RUN_narrow") == "narrow"

    def test_cut_tag_broad(self):
        assert get_peaktype("CUT_TAG_broad") == "broad"

    def test_cut_tag_narrow(self):
        assert get_peaktype("CUT_TAG_narrow") == "narrow"

    def test_config_override(self):
        assert get_peaktype("ChIP_broad", {"ChIP_broad": "narrow"}) == "narrow"

    def test_no_peaktype(self):
        with pytest.raises(ValueError, match="No peak type"):
            get_peaktype("RNAseq")


class TestCUTAssayVocabulary:
    """Sample_sheet should treat CUT&RUN/CUT&Tag as IP_PEAK_ASSAYS routed to the ChIP env."""

    def test_all_four_in_valid_assays(self):
        for a in ("CUT_RUN_broad", "CUT_RUN_narrow",
                  "CUT_TAG_broad", "CUT_TAG_narrow"):
            assert a in VALID_ASSAYS

    def test_all_four_route_to_chip_env(self):
        for a in ("CUT_RUN_broad", "CUT_RUN_narrow",
                  "CUT_TAG_broad", "CUT_TAG_narrow"):
            assert ASSAY_TO_ENV[a] == "ChIP"

    def test_peaktype_suffix_matches_designation(self):
        for a in ("CUT_RUN_broad", "CUT_TAG_broad"):
            assert ASSAY_TO_PEAKTYPE[a] == "broad"
        for a in ("CUT_RUN_narrow", "CUT_TAG_narrow"):
            assert ASSAY_TO_PEAKTYPE[a] == "narrow"

    def test_in_ip_peak_assays_set(self):
        from workflow.scripts.sample_sheet import IP_PEAK_ASSAYS
        for a in ("ChIP_broad", "ChIP_narrow",
                  "CUT_RUN_broad", "CUT_RUN_narrow",
                  "CUT_TAG_broad", "CUT_TAG_narrow"):
            assert a in IP_PEAK_ASSAYS
        # Non-IP assays must not leak into the set
        for a in ("ATAC", "RNAseq", "RAMPAGE", "sRNA", "WGBS", "dmC"):
            assert a not in IP_PEAK_ASSAYS


# ---------------------------------------------------------------------------
# Compat columns
# ---------------------------------------------------------------------------

class TestAddCompatColumns:
    def test_basic_columns_exist(self, pombe_df):
        result = add_compat_columns(pombe_df)
        for col in ["data_type", "ref_genome", "replicate", "paired",
                     "line", "tissue", "sample_type", "extra_info",
                     "levels_label", "env", "sample_name",
                     "seq_id", "fastq_path"]:
            assert col in result.columns, f"Missing column: {col}"

    def test_sample_name_is_sample_id(self, pombe_df):
        result = add_compat_columns(pombe_df)
        assert (result["sample_name"] == result["Sample_ID"]).all()

    def test_env_mapping(self, pombe_df):
        result = add_compat_columns(pombe_df)
        chip_rows = result[result["Assay"] == "ChIP_broad"]
        assert (chip_rows["env"] == "ChIP").all()

    def test_line_from_levels(self, pombe_df):
        result = add_compat_columns(pombe_df)
        wt_row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert wt_row["line"] == "WT"
        assert wt_row["tissue"] == ""

    def test_sample_type_chip(self, pombe_df):
        result = add_compat_columns(pombe_df)
        wt_chip = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert wt_chip["sample_type"] == "H3K9me2"

    def test_sample_type_rna(self, pombe_df):
        result = add_compat_columns(pombe_df)
        rna_row = result[result["Sample_ID"] == "WT_RNA_rep1"].iloc[0]
        assert rna_row["sample_type"] == "RNAseq"

    def test_extra_info_chip(self, pombe_df):
        result = add_compat_columns(pombe_df)
        chip_row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert chip_row["extra_info"] == "H3K9me2"

    def test_extra_info_rna(self, pombe_df):
        result = add_compat_columns(pombe_df)
        rna_row = result[result["Sample_ID"] == "WT_RNA_rep1"].iloc[0]
        assert rna_row["extra_info"] == "N/A"

    def test_sra_seq_id(self, pombe_df):
        result = add_compat_columns(pombe_df)
        row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert row["seq_id"] == "SRR20678305"
        assert row["fastq_path"] == "SRA"

    def test_levels_label(self, pombe_df):
        result = add_compat_columns(pombe_df)
        row = result[result["Sample_ID"] == "WT_H3K9me2_rep1"].iloc[0]
        assert row["levels_label"] == "WT"


class TestReadSampleSheet:
    HEADER = "Sample_ID\tAssay\tGenome\tLevels\tReplicate_ID\tRead_files\tRead_layout\tIP_target\tControl"
    ROW_A = "A\tChIP_broad\tColCEN\tgenotype:WT\trep1\tSRR000001\tSE\tH3K9me2\tInputA"
    ROW_B = "B\tChIP_broad\tColCEN\tgenotype:WT\trep2\tSRR000002\tSE\tH3K9me2\tInputA"
    ROW_C = "InputA\tChIP_broad\tColCEN\tgenotype:WT\trep1\tSRR000003\tSE\tInput\t"

    def test_full_line_comments_skipped(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text(
            "\n".join([
                "# this is a header comment",
                self.HEADER,
                "# parked sample, do not run",
                self.ROW_A,
                self.ROW_C,
            ]) + "\n"
        )
        df = read_sample_sheet(f)
        assert sorted(df["Sample_ID"]) == ["A", "InputA"]

    def test_no_comments_unaffected(self, tmp_path):
        f = tmp_path / "samples.tsv"
        f.write_text("\n".join([self.HEADER, self.ROW_A, self.ROW_B, self.ROW_C]) + "\n")
        df = read_sample_sheet(f)
        assert sorted(df["Sample_ID"]) == ["A", "B", "InputA"]
