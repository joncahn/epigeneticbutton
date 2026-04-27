"""Unit tests for workflow/scripts/samplefile_validation.py

Focuses on the file-existence rules added for Read_files paths and genome
config file paths. SRA accessions and HTTP(S) URLs are skipped — only
local-filesystem paths are probed.
"""

import os
import sys

import pandas as pd
import pytest

# samplefile_validation imports `from scripts.sample_sheet import ...`,
# which only resolves when the parent of the workflow/ tree is on sys.path
# (mirroring how Snakemake adds workflow.basedir's parent at runtime).
_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from scripts.samplefile_validation import check_table, check_genome_config


def _row(sample_id, assay, read_files, layout="SE", **overrides):
    base = {
        "Sample_ID": sample_id,
        "Assay": assay,
        "Genome": "test_genome",
        "Levels": "genotype:WT",
        "Replicate_ID": "rep1",
        "Read_files": read_files,
        "Read_layout": layout,
        "IP_target": "",
        "Control": "",
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# Read_files: local path existence
# ---------------------------------------------------------------------------

class TestReadFilesPathExistence:
    def test_existing_local_path_passes(self, tmp_path):
        f = tmp_path / "reads.fq.gz"
        f.write_bytes(b"")
        df = pd.DataFrame([_row("s1", "RNAseq", str(f))])
        check_table(df)  # no raise

    def test_missing_local_path_raises(self, tmp_path):
        missing = tmp_path / "missing.fq.gz"
        df = pd.DataFrame([_row("s1", "RNAseq", str(missing))])
        with pytest.raises(ValueError) as excinfo:
            check_table(df)
        assert "does not exist" in str(excinfo.value)
        assert str(missing) in str(excinfo.value)

    def test_missing_pe_mate_raises(self, tmp_path):
        r1 = tmp_path / "r1.fq.gz"
        r1.write_bytes(b"")
        r2 = tmp_path / "r2.fq.gz"  # never created
        rf = f"{r1},{r2}"
        df = pd.DataFrame([_row("s1", "RNAseq", rf, layout="PE")])
        with pytest.raises(ValueError) as excinfo:
            check_table(df)
        msg = str(excinfo.value)
        assert "r2.fq.gz" in msg
        assert "does not exist" in msg

    def test_missing_merge_component_raises(self, tmp_path):
        a = tmp_path / "a.fq.gz"
        a.write_bytes(b"")
        b = tmp_path / "b.fq.gz"  # never created
        df = pd.DataFrame([_row("s1", "RNAseq", f"{a}+{b}")])
        with pytest.raises(ValueError) as excinfo:
            check_table(df)
        assert "b.fq.gz" in str(excinfo.value)

    def test_sra_accession_not_probed(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "SRR1234567")])
        check_table(df)  # no raise — SRA IDs are never probed on disk

    def test_ena_drr_not_probed(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "DRR400324", layout="PE")])
        check_table(df)

    def test_http_url_not_probed(self):
        df = pd.DataFrame([
            _row("s1", "dmC", "https://example.org/data.bam"),
        ])
        check_table(df)

    def test_check_paths_false_skips_existence(self, tmp_path):
        missing = tmp_path / "missing.fq.gz"
        df = pd.DataFrame([_row("s1", "RNAseq", str(missing))])
        # Bypass: caller signals that inputs may not be staged yet
        check_table(df, check_paths=False)


# ---------------------------------------------------------------------------
# Genome config: file existence
# ---------------------------------------------------------------------------

def _minimal_sample_df(genome="test_genome", assay="RNAseq"):
    return pd.DataFrame([_row("s1", assay, "SRR1234567")]).assign(Genome=genome)


class TestGenomeConfigPathExistence:
    def _genome_cfg(self, **fields):
        cfg = {"genomes": {"test_genome": fields}}
        return cfg

    def test_missing_fasta_raises(self, tmp_path):
        gff = tmp_path / "annot.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(
            fasta_file=str(tmp_path / "absent.fa"),
            gff_file=str(gff),
        )
        with pytest.raises(ValueError) as excinfo:
            check_genome_config(_minimal_sample_df(), cfg)
        assert "fasta_file" in str(excinfo.value)
        assert "does not exist" in str(excinfo.value)

    def test_existing_fasta_and_gff_pass(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_bytes(b"")
        gff = tmp_path / "g.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(fasta_file=str(fa), gff_file=str(gff))
        check_genome_config(_minimal_sample_df(), cfg)

    def test_auto_sentinel_skipped(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_bytes(b"")
        gff = tmp_path / "g.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(
            fasta_file=str(fa), gff_file=str(gff),
            gtf_file="<auto>",
        )
        check_genome_config(_minimal_sample_df(), cfg)

    def test_url_not_probed_for_existence(self, tmp_path):
        cfg = self._genome_cfg(
            fasta_file="https://example.org/g.fa",
            gff_file="https://example.org/g.gff3",
        )
        check_genome_config(_minimal_sample_df(), cfg)

    def test_te_file_existence(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_bytes(b"")
        gff = tmp_path / "g.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(
            fasta_file=str(fa), gff_file=str(gff),
            te_file=str(tmp_path / "absent.te.bed"),
        )
        with pytest.raises(ValueError) as excinfo:
            check_genome_config(_minimal_sample_df(), cfg)
        assert "te_file" in str(excinfo.value)

    def test_check_paths_false_skips_existence(self, tmp_path):
        cfg = self._genome_cfg(
            fasta_file=str(tmp_path / "absent.fa"),
            gff_file=str(tmp_path / "absent.gff3"),
        )
        check_genome_config(_minimal_sample_df(), cfg, check_paths=False)
