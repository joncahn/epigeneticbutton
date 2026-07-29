"""Unit tests for --config merging in build_snakemake_cmd (issue #31).

When the user passes `-- --config key=val` via epicc's passthrough, Snakemake
would receive two separate --config flags and use only the last one, silently
dropping epicc's own keys (e.g. sample_file). The fix merges all --config
values into a single block.
"""

import importlib.util
import sys
import types
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EPICC_SRC = _REPO_ROOT / "epicc"


def _load_epicc():
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader("epicc_cli_cfg", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli_cfg", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli_cfg", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


def _make_args(samples=None, output_dir=None, genome_dir=None,
               keep_intermediates=None, use_node_tmpdir=False,
               options="config/epicc-options.yaml",
               profile=None, cores=None, verbose=False, quiet=False,
               workflow_profile=None):
    """Minimal args namespace for build_snakemake_cmd."""
    ns = types.SimpleNamespace(
        samples=samples,
        output_dir=output_dir,
        genome_dir=genome_dir,
        keep_intermediates=keep_intermediates,
        use_node_tmpdir=use_node_tmpdir,
        options=options,
        profile=profile,
        cores=cores,
        verbose=verbose,
        quiet=quiet,
        workflow_profile=workflow_profile,
        no_rerun_incomplete=True,   # keep test output deterministic
    )
    return ns


def _config_values(cmd):
    """Return the list of key=value tokens that follow --config in cmd."""
    try:
        idx = cmd.index("--config")
    except ValueError:
        return []
    values = []
    for tok in cmd[idx + 1:]:
        if tok.startswith("-"):
            break
        values.append(tok)
    return values


def _config_count(cmd):
    return cmd.count("--config")


class TestConfigMerge:
    def test_no_passthrough_emits_one_config_block(self, tmp_path):
        samples = tmp_path / "s.tsv"
        samples.touch()
        args = _make_args(samples=str(samples), options=str(tmp_path / "opts.yaml"))
        (tmp_path / "opts.yaml").touch()
        cmd = epicc.build_snakemake_cmd(args, extra_args=[], dry_run=True)
        assert _config_count(cmd) == 1
        assert f"sample_file={samples}" in _config_values(cmd)

    def test_passthrough_config_merged_into_single_block(self, tmp_path):
        samples = tmp_path / "s.tsv"
        samples.touch()
        (tmp_path / "opts.yaml").touch()
        args = _make_args(samples=str(samples), options=str(tmp_path / "opts.yaml"))
        extra = ["--config", "analysis_name=myrun"]
        cmd = epicc.build_snakemake_cmd(args, extra_args=extra, dry_run=True)
        assert _config_count(cmd) == 1, "must be a single --config block"
        vals = _config_values(cmd)
        assert f"sample_file={samples}" in vals, "sample_file must survive merge"
        assert "analysis_name=myrun" in vals, "user key must be present"

    def test_passthrough_config_equals_form_merged(self, tmp_path):
        samples = tmp_path / "s.tsv"
        samples.touch()
        (tmp_path / "opts.yaml").touch()
        args = _make_args(samples=str(samples), options=str(tmp_path / "opts.yaml"))
        extra = ["--config=analysis_name=myrun"]
        cmd = epicc.build_snakemake_cmd(args, extra_args=extra, dry_run=True)
        assert _config_count(cmd) == 1
        vals = _config_values(cmd)
        assert f"sample_file={samples}" in vals
        assert "analysis_name=myrun" in vals

    def test_non_config_passthrough_preserved(self, tmp_path):
        (tmp_path / "opts.yaml").touch()
        args = _make_args(options=str(tmp_path / "opts.yaml"))
        extra = ["--config", "analysis_name=x", "--dry-run", "--quiet"]
        cmd = epicc.build_snakemake_cmd(args, extra_args=extra, dry_run=False)
        assert "--dry-run" in cmd
        assert "--quiet" in cmd
        assert _config_count(cmd) == 1

    def test_multiple_passthrough_config_keys_all_merged(self, tmp_path):
        samples = tmp_path / "s.tsv"
        samples.touch()
        (tmp_path / "opts.yaml").touch()
        args = _make_args(samples=str(samples), options=str(tmp_path / "opts.yaml"))
        extra = ["--config", "analysis_name=x", "output_dir=out2"]
        cmd = epicc.build_snakemake_cmd(args, extra_args=extra, dry_run=True)
        assert _config_count(cmd) == 1
        vals = _config_values(cmd)
        assert f"sample_file={samples}" in vals
        assert "analysis_name=x" in vals
        assert "output_dir=out2" in vals
