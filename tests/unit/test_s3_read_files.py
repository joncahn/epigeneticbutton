"""Unit tests for public-S3 (``s3://bucket/key``) support in Read_files.

EPICC supports authentication-free (public) S3 objects only: no signing, no
credentials, no aws CLI. A public S3 object is an ordinary HTTPS GET, so
support amounts to recognising the scheme and translating it to the public
HTTPS endpoint at the one boundary where download paths are produced
(``get_seq_id_and_path``). The download rules stay S3-unaware.
"""

import os
import sys

import pandas as pd
import pytest

_REPO_ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow", "scripts"))
# samplefile_validation imports `from scripts.sample_sheet import ...`, so the
# workflow/ dir must be importable too (mirrors test_samplefile_validation.py).
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from sample_sheet import (  # noqa: E402
    get_seq_id_and_path,
    s3_uri_to_https,
)
from scripts.samplefile_validation import (  # noqa: E402
    _is_local_path,
    _merge_component_kind,
    check_table,
)


# ---------------------------------------------------------------------------
# s3:// -> https translation
# ---------------------------------------------------------------------------

class TestS3UriToHttps:
    def test_virtual_hosted_style_for_plain_bucket(self):
        assert s3_uri_to_https("s3://bucket/key/reads.fq.gz") == \
            "https://bucket.s3.amazonaws.com/key/reads.fq.gz"

    def test_nested_key_preserved(self):
        assert s3_uri_to_https("s3://bucket/a/b/c/aln.bam") == \
            "https://bucket.s3.amazonaws.com/a/b/c/aln.bam"

    def test_dotted_bucket_uses_path_style(self):
        # The wildcard cert for *.s3.amazonaws.com does not match a dotted
        # bucket label, so virtual-hosted style would fail TLS verification.
        assert s3_uri_to_https("s3://my.dotted.bucket/k/r.fq.gz") == \
            "https://s3.amazonaws.com/my.dotted.bucket/k/r.fq.gz"

    def test_region_less_endpoint_emitted(self):
        # No region is guessed: S3 resolves/redirects the global endpoint and
        # every curl in the workflow passes --location.
        assert ".s3.amazonaws.com" in s3_uri_to_https("s3://b/k.fq.gz")
        assert "us-east-1" not in s3_uri_to_https("s3://b/k.fq.gz")

    def test_query_string_preserved(self):
        assert s3_uri_to_https("s3://b/k/r.fq.gz?v=2") == \
            "https://b.s3.amazonaws.com/k/r.fq.gz?v=2"

    @pytest.mark.parametrize("value", [
        "https://example.com/x.fq.gz",
        "/local/path/r.fq.gz",
        "SRR12345",
    ])
    def test_non_s3_returned_unchanged(self, value):
        assert s3_uri_to_https(value) == value

    @pytest.mark.parametrize("bad", ["s3://bucket", "s3://bucket/", "s3:///key"])
    def test_malformed_uri_raises(self, bad):
        with pytest.raises(ValueError, match="Malformed S3 URI"):
            s3_uri_to_https(bad)

    def test_matches_ncbi_published_url_for_real_run(self):
        # issue #47's example accession. NCBI's SDL API publishes this exact
        # HTTPS link for the run in the AWS Open Data bucket sra-pub-run-odp,
        # so our translation is pinned against a real, externally-defined URL
        # rather than only against itself. Verified fetchable (HTTP 206 on a
        # range request, no credentials).
        assert s3_uri_to_https("s3://sra-pub-run-odp/sra/SRR35264187/SRR35264187") == \
            "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR35264187/SRR35264187"


# ---------------------------------------------------------------------------
# get_seq_id_and_path: seq_id derivation + type detection + translation
# ---------------------------------------------------------------------------

class TestS3SeqIdAndPath:
    def test_fastq_se_uses_url_sentinel_and_translates(self):
        seq_id, path = get_seq_id_and_path("s3://bucket/k/reads.fq.gz", "SE")
        # "URL" is what makes the download rules take their curl branch.
        assert seq_id == "URL"
        assert path == "https://bucket.s3.amazonaws.com/k/reads.fq.gz"

    def test_fastq_pe_translates_both_mates(self):
        seq_id, path = get_seq_id_and_path(
            "s3://b/k/r1.fq.gz,s3://b/k/r2.fq.gz", "PE")
        assert seq_id == "URL"
        assert path == ("https://b.s3.amazonaws.com/k/r1.fq.gz,"
                        "https://b.s3.amazonaws.com/k/r2.fq.gz")

    def test_fastq_merge_translates_every_component(self):
        seq_id, path = get_seq_id_and_path(
            "s3://b/k/a.fq.gz+s3://b/k/b.fq.gz", "SE")
        assert seq_id == "URL"
        assert path == ("https://b.s3.amazonaws.com/k/a.fq.gz+"
                        "https://b.s3.amazonaws.com/k/b.fq.gz")

    def test_mixed_s3_and_https_merge(self):
        _, path = get_seq_id_and_path(
            "s3://b/k/a.fq.gz+https://h/b.fq.gz", "SE")
        assert path == "https://b.s3.amazonaws.com/k/a.fq.gz+https://h/b.fq.gz"

    def test_bam_seq_id_and_translation(self):
        seq_id, path = get_seq_id_and_path("s3://bucket/key/aln.bam", "SE")
        assert seq_id == "aln"
        assert path == "https://bucket.s3.amazonaws.com/key/aln.bam"

    def test_bedmethyl_seq_id_and_translation(self):
        seq_id, path = get_seq_id_and_path("s3://bucket/key/m.bed.gz", "SE")
        assert seq_id == "m"
        assert path == "https://bucket.s3.amazonaws.com/key/m.bed.gz"

    def test_query_string_does_not_defeat_type_detection(self):
        # Regression: before s3:// was recognised as remote, the query string
        # was not stripped and this raised "Unrecognized Read_files format".
        seq_id, path = get_seq_id_and_path("s3://b/k/reads.fq.gz?v=2", "SE")
        assert seq_id == "URL"
        assert path == "https://b.s3.amazonaws.com/k/reads.fq.gz?v=2"

    def test_unknown_extension_still_rejected(self):
        with pytest.raises(ValueError, match="Unrecognized Read_files format"):
            get_seq_id_and_path("s3://bucket/key/notes.txt", "SE")

    @pytest.mark.parametrize("value,layout,expected", [
        ("https://h/r.fq.gz", "SE", ("URL", "https://h/r.fq.gz")),
        ("/local/r.fq.gz", "SE", ("EXPLICIT", "/local/r.fq.gz")),
        ("SRR12345", "SE", ("SRR12345", "SRA")),
    ])
    def test_non_s3_inputs_unaffected(self, value, layout, expected):
        assert get_seq_id_and_path(value, layout) == expected


# ---------------------------------------------------------------------------
# Validation: s3:// is remote (never stat'd), and merges correctly
# ---------------------------------------------------------------------------

class TestS3Validation:
    def test_s3_is_not_a_local_path(self):
        # Otherwise validation would try to stat it on disk and fail.
        assert _is_local_path("s3://bucket/k/r.fq.gz") is False

    @pytest.mark.parametrize("value,kind", [
        ("s3://b/k/r.fq.gz", "fastq"),
        ("s3://b/k/a.bam", "bam"),
        ("s3://b/k/m.bed.gz", "bedmethyl"),
        ("s3://b/k/m.bedmethyl.gz", "bedmethyl"),
    ])
    def test_merge_component_kind(self, value, kind):
        assert _merge_component_kind(value) == kind

    @staticmethod
    def _row(sample_id, read_files, layout="SE"):
        return {
            "Sample_ID": sample_id, "Assay": "RNAseq", "Genome": "test_genome",
            "Levels": "genotype:WT", "Replicate_ID": "rep1",
            "Read_files": read_files, "Read_layout": layout,
            "IP_target": "", "Control": "",
        }

    def test_s3_fastq_accepted_without_touching_disk(self):
        df = pd.DataFrame([self._row("s1", "s3://bucket/k/reads.fq.gz")])
        check_table(df)  # check_paths defaults True — must not stat the URI

    def test_s3_fastq_merge_accepted(self):
        df = pd.DataFrame([
            self._row("s1", "s3://b/k/a.fq.gz+s3://b/k/b.fq.gz"),
        ])
        check_table(df)

    def test_s3_bam_merge_rejected(self):
        # BAM needs samtools merge, so '+'-merge is rejected regardless of
        # whether the components are local, HTTP or S3.
        df = pd.DataFrame([
            self._row("s1", "s3://b/k/a.bam+s3://b/k/b.bam"),
        ])
        with pytest.raises(ValueError, match="not supported for bam"):
            check_table(df)

    def test_malformed_s3_uri_reported_as_validation_error(self):
        df = pd.DataFrame([self._row("s1", "s3://bucket-only")])
        with pytest.raises(ValueError, match="malformed S3 URI"):
            check_table(df)


class TestSraVersionSuffix:
    """SRA "Original Format" objects carry a trailing .<digits> version suffix
    (e.g. ``..._Nanopore.bam.1``), which otherwise hides the real extension.

    Note these objects live in Requester-Pays buckets and so cannot actually be
    fetched anonymously — see the follow-up issue. Type detection is fixed here
    so the remaining gap is purely credentials, not parsing.
    """

    @pytest.mark.parametrize("name,expected_seq_id", [
        ("Nvec_F1_Gastrula_GSKxUntreated_Nanopore.bam.1",
         "Nvec_F1_Gastrula_GSKxUntreated_Nanopore"),
        ("sample.bam.1", "sample"),
        ("m.bed.gz.1", "m"),
    ])
    def test_version_suffix_reveals_real_extension(self, name, expected_seq_id):
        seq_id, path = get_seq_id_and_path(f"s3://bucket/key/{name}", "SE")
        assert seq_id == expected_seq_id
        # The download path keeps the true object name, suffix and all.
        assert path.endswith(name)

    def test_versioned_fastq_uses_url_sentinel(self):
        seq_id, path = get_seq_id_and_path("s3://bucket/k/reads.fastq.gz.2", "SE")
        assert seq_id == "URL"
        assert path.endswith("reads.fastq.gz.2")

    def test_local_path_version_suffix_also_handled(self):
        seq_id, path = get_seq_id_and_path("/data/sample.bam.1", "SE")
        assert seq_id == "sample"
        assert path == "/data/sample.bam.1"

    @pytest.mark.parametrize("name", ["notafile.1", "chunk.12", "archive.7"])
    def test_bare_numeric_suffix_still_rejected(self, name):
        # Stripping must only ever *reveal* a known extension — never invent one.
        with pytest.raises(ValueError, match="Unrecognized Read_files format"):
            get_seq_id_and_path(f"s3://bucket/key/{name}", "SE")
