"""Unit tests for the `epicc output` subcommand helpers.

The `epicc` script is a top-level executable (no .py extension), so we
load it as a module via importlib.
"""

import importlib.util
import os
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EPICC_SRC = _REPO_ROOT / "epicc"


def _load_epicc():
    # The `epicc` script has no .py extension, so spec_from_file_location
    # can't infer a loader on its own — pass SourceFileLoader explicitly.
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader("epicc_cli", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


# ---------------------------------------------------------------------------
# derive_label_from_filename
# ---------------------------------------------------------------------------

class TestDeriveLabel:
    def test_simple_filename(self):
        assert epicc.derive_label_from_filename("/data/target_genes.txt") == "target_genes"

    def test_double_extension(self):
        assert epicc.derive_label_from_filename("/data/regions.bed.gz") == "regions"

    def test_no_extension(self):
        assert epicc.derive_label_from_filename("/data/notes") == "notes"

    def test_relative(self):
        assert epicc.derive_label_from_filename("foo.bar.baz") == "foo"


# ---------------------------------------------------------------------------
# Output type registry
# ---------------------------------------------------------------------------

class TestOutputRegistry:
    def test_all_documented_types_present(self):
        # The TODO mentions these seven plot types; make sure none get
        # accidentally renamed/removed.
        expected = {
            "rnaseq-histogram", "go", "motifs", "srna-clusters",
            "heatmap", "metaplot", "browser",
        }
        assert set(epicc.OUTPUT_TYPES.keys()) == expected

    def test_target_templates_use_known_placeholders(self):
        allowed = {"output_dir", "analysis_name", "ref_genome",
                   "env", "label", "matrix"}
        import string
        for name, spec in epicc.OUTPUT_TYPES.items():
            keys = {fname for _, fname, _, _ in
                    string.Formatter().parse(spec["target"])
                    if fname}
            assert keys <= allowed, (
                f"OUTPUT_TYPES[{name!r}].target uses unknown "
                f"placeholders: {keys - allowed}"
            )


# ---------------------------------------------------------------------------
# validate_output_input
# ---------------------------------------------------------------------------

class TestValidateOutputInput:
    def test_missing_file_exits(self, tmp_path):
        missing = tmp_path / "nope.txt"
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(missing), "tsv_geneids")
        assert "not found" in str(exc.value)

    def test_empty_file_exits(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("")
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(f), "tsv_geneids")
        assert "empty" in str(exc.value)

    def test_only_comments_exits(self, tmp_path):
        f = tmp_path / "comments.txt"
        f.write_text("# just a header\n# another comment\n")
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(f), "tsv_geneids")
        assert "no data rows" in str(exc.value)

    def test_tsv_geneids_passes(self, tmp_path):
        f = tmp_path / "genes.txt"
        f.write_text("GENE1\tlabelA\nGENE2\tlabelB\n")
        epicc.validate_output_input(str(f), "tsv_geneids")  # no raise

    def test_tsv_geneids_empty_first_col_exits(self, tmp_path):
        f = tmp_path / "genes.txt"
        f.write_text("\tlabelA\n")
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(f), "tsv_geneids")
        assert "first column" in str(exc.value)

    def test_bed_passes(self, tmp_path):
        f = tmp_path / "regions.bed"
        f.write_text("chrI\t100\t200\tname1\nchrI\t500\t800\tname2\n")
        epicc.validate_output_input(str(f), "bed")

    def test_bed_too_few_cols_exits(self, tmp_path):
        f = tmp_path / "bad.bed"
        f.write_text("chrI\t100\n")
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(f), "bed")
        assert "BED" in str(exc.value)

    def test_bed_with_header_line_tolerated(self, tmp_path):
        # Non-numeric coords on first line: probably a header, not an error.
        # Snakemake's has_header checkpoint handles this downstream.
        f = tmp_path / "with_header.bed"
        f.write_text("chrom\tstart\tend\tname\nchrI\t100\t200\tn1\n")
        epicc.validate_output_input(str(f), "bed")

    def test_bed_negative_start_exits(self, tmp_path):
        f = tmp_path / "neg.bed"
        f.write_text("chrI\t-1\t200\n")
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(f), "bed")
        assert "BED coordinates" in str(exc.value)

    def test_bed_end_before_start_exits(self, tmp_path):
        f = tmp_path / "rev.bed"
        f.write_text("chrI\t500\t100\n")
        with pytest.raises(SystemExit) as exc:
            epicc.validate_output_input(str(f), "bed")

    def test_loci_format_passes(self, tmp_path):
        f = tmp_path / "loci.gff3"
        f.write_text("chr1\tShortStack\tcluster\t1\t1000\t.\t+\t.\tID=c1\n")
        epicc.validate_output_input(str(f), "loci")

    def test_gzipped_input(self, tmp_path):
        import gzip
        f = tmp_path / "regions.bed.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("chrI\t10\t100\tn1\n")
        epicc.validate_output_input(str(f), "bed")
