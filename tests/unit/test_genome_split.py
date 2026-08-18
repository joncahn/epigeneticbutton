"""Genome as a separate token: multi-genome sheets and mapped_name (issue #39).

Sample_ID identifies the library, not the alignment. Read processing stays
genome-free so one download+trim serves every reference; everything from
alignment onward carries the genome in `mapped_name`.
"""

import os
import sys

import pandas as pd
import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from sample_sheet import (  # noqa: E402
    add_compat_columns,
    explode_genomes,
    parse_genomes,
    read_sample_sheet,
)
from scripts.samplefile_validation import check_table  # noqa: E402

HEADER = ("Sample_ID\tAssay\tGenome\tLevels\tReplicate_ID\tRead_files\t"
          "Read_layout\tIP_target\tControl")


def _sheet(tmp_path, rows, name="s.tsv"):
    f = tmp_path / name
    f.write_text("\n".join([HEADER] + rows) + "\n")
    return f


class TestParseGenomes:
    @pytest.mark.parametrize("raw,expected", [
        ("ColCEN", ["ColCEN"]),
        ("B73,W22", ["B73", "W22"]),
        (" B73 , W22 ", ["B73", "W22"]),
        ("B73,,W22", ["B73", "W22"]),
    ])
    def test_splits(self, raw, expected):
        assert parse_genomes(raw) == expected


class TestExplode:
    def test_one_row_per_genome(self):
        df = pd.DataFrame([{"Sample_ID": "s1", "Genome": "B73,W22"}])
        out = explode_genomes(df)
        assert list(out["Genome"]) == ["B73", "W22"]
        assert list(out["Sample_ID"]) == ["s1", "s1"]

    def test_single_genome_untouched(self):
        df = pd.DataFrame([{"Sample_ID": "s1", "Genome": "B73"}])
        assert list(explode_genomes(df)["Genome"]) == ["B73"]

    def test_repeated_reference_rejected(self):
        # After the explode a repeat is indistinguishable from two real rows,
        # so it has to be caught here.
        df = pd.DataFrame([{"Sample_ID": "s1", "Genome": "B73,B73"}])
        with pytest.raises(ValueError, match="more than once"):
            explode_genomes(df)


class TestMappedName:
    def test_mapped_name_is_sample_and_genome(self, tmp_path):
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73,W22\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        df = add_compat_columns(read_sample_sheet(f))
        assert sorted(df["mapped_name"]) == ["s1__B73", "s1__W22"]
        # sample_name stays bare: it is the pre-alignment identity.
        assert set(df["sample_name"]) == {"s1"}

    def test_single_genome_mapped_name(self, tmp_path):
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        df = add_compat_columns(read_sample_sheet(f))
        assert list(df["mapped_name"]) == ["s1__B73"]


class TestValidation:
    def test_multi_genome_sheet_validates(self, tmp_path):
        # Same Sample_ID and same reads on both exploded rows is the whole point.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73,W22\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        check_table(read_sample_sheet(f), check_paths=False)  # no raise

    def test_genome_with_double_underscore_rejected(self, tmp_path):
        # '__' is the delimiter in '{Sample_ID}__{Genome}'.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73__v5\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        with pytest.raises(ValueError, match="must not contain '__'"):
            check_table(read_sample_sheet(f), check_paths=False)

    def test_genuine_duplicate_still_rejected(self, tmp_path):
        # Same Sample_ID AND same genome is still a real duplicate.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
            "s1\tRNAseq\tB73\tgenotype:WT\trep2\tSRR2\tSE\t\t",
        ])
        with pytest.raises(ValueError, match="duplicate Sample_ID"):
            check_table(read_sample_sheet(f), check_paths=False)

    def test_shared_reads_across_samples_still_rejected(self, tmp_path):
        # Two DIFFERENT samples pointing at one accession remains an error.
        f = _sheet(tmp_path, [
            "s1\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
            "s2\tRNAseq\tB73\tgenotype:WT\trep1\tSRR1\tSE\t\t",
        ])
        with pytest.raises(ValueError, match="is also used by"):
            check_table(read_sample_sheet(f), check_paths=False)
