"""Unit tests for the `epicc` snakemake version guard.

The guard exists because the workflow's own `min_version("9.0")` can never
fire against an old snakemake: version 7 rejects `--executor` while parsing
arguments, so the user gets "unrecognized arguments: --executor=slurm" and no
hint that the environment is the problem. The `epicc` script has no .py
extension, so we load it via importlib (same pattern as
test_epicc_placeholders.py).
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EPICC_SRC = _REPO_ROOT / "epicc"


def _load_epicc():
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader("epicc_cli", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


# ---------------------------------------------------------------------------
# parse_version
# ---------------------------------------------------------------------------

class TestParseVersion:
    @pytest.mark.parametrize("text,expected", [
        ("9.21.0", (9, 21, 0)),
        ("7.32.4", (7, 32, 4)),
        ("9.0.0", (9, 0, 0)),
        ("v9.21.0", (9, 21, 0)),
        ("  9.21.0  ", (9, 21, 0)),
        ("9.21.0\n", (9, 21, 0)),
    ])
    def test_plain_versions(self, text, expected):
        assert epicc.parse_version(text) == expected

    def test_two_components_pad_to_three(self):
        assert epicc.parse_version("9.1") == (9, 1, 0)

    def test_single_component_pads_to_three(self):
        # Must NOT compare as older than 9.0.0 just for being shorter.
        assert epicc.parse_version("9") == (9, 0, 0)
        assert epicc.parse_version("9") >= epicc.MIN_SNAKEMAKE

    def test_prerelease_suffix_is_truncated(self):
        assert epicc.parse_version("9.1.0rc1") == (9, 1, 0)

    def test_extra_components_ignored(self):
        assert epicc.parse_version("9.21.0.1") == (9, 21, 0)

    def test_trailing_words_ignored(self):
        assert epicc.parse_version("9.21.0 (some build)") == (9, 21, 0)

    @pytest.mark.parametrize("text", ["", "   ", "unknown", "vX.Y.Z"])
    def test_unparseable_returns_none(self, text):
        assert epicc.parse_version(text) is None


# ---------------------------------------------------------------------------
# check_snakemake_version
# ---------------------------------------------------------------------------

class _FakeProc:
    def __init__(self, stdout="", stderr=""):
        self.stdout = stdout
        self.stderr = stderr


@pytest.fixture
def fake_snakemake(monkeypatch):
    """Point the guard at a fabricated `snakemake --version` result."""
    def _install(version_text, path="/fake/env/bin/snakemake"):
        monkeypatch.setattr(epicc.shutil, "which", lambda _cmd: path)
        monkeypatch.setattr(
            epicc.subprocess, "run",
            lambda *a, **k: _FakeProc(stdout=version_text),
        )
    return _install


class TestCheckSnakemakeVersion:
    def test_current_version_passes(self, fake_snakemake):
        fake_snakemake("9.21.0")
        epicc.check_snakemake_version()   # must not raise

    def test_exact_minimum_passes(self, fake_snakemake):
        fake_snakemake("9.0.0")
        epicc.check_snakemake_version()

    def test_jons_version_is_rejected(self, fake_snakemake):
        fake_snakemake("7.32.4")
        with pytest.raises(SystemExit) as exc:
            epicc.check_snakemake_version()
        msg = str(exc.value)
        assert "7.32.4" in msg          # what was found
        assert "9.0.0" in msg           # what is required
        assert "/fake/env/bin/snakemake" in msg   # which one was used

    def test_snakemake_8_is_rejected(self, fake_snakemake):
        fake_snakemake("8.25.1")
        with pytest.raises(SystemExit):
            epicc.check_snakemake_version()

    def test_error_mentions_case_sensitivity(self, fake_snakemake):
        # The EPICC/epicc collision is the trap this guard is meant to expose.
        fake_snakemake("7.32.4")
        with pytest.raises(SystemExit) as exc:
            epicc.check_snakemake_version()
        assert "case-sensitive" in str(exc.value)

    def test_missing_snakemake_exits_with_guidance(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _cmd: None)
        with pytest.raises(SystemExit) as exc:
            epicc.check_snakemake_version()
        assert "not found on PATH" in str(exc.value)

    def test_version_read_from_stderr(self, monkeypatch):
        # Some builds print the version to stderr.
        monkeypatch.setattr(epicc.shutil, "which", lambda _cmd: "/x/snakemake")
        monkeypatch.setattr(
            epicc.subprocess, "run",
            lambda *a, **k: _FakeProc(stdout="", stderr="7.32.4"),
        )
        with pytest.raises(SystemExit):
            epicc.check_snakemake_version()

    def test_unparseable_version_warns_but_continues(self, fake_snakemake, capsys):
        fake_snakemake("some-dev-build")
        epicc.check_snakemake_version()   # must not raise
        assert "skipping the snakemake version check" in capsys.readouterr().err

    def test_subprocess_failure_exits_cleanly(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _cmd: "/x/snakemake")

        def _boom(*a, **k):
            raise OSError("exec format error")

        monkeypatch.setattr(epicc.subprocess, "run", _boom)
        with pytest.raises(SystemExit) as exc:
            epicc.check_snakemake_version()
        assert "exec format error" in str(exc.value)

    def test_timeout_exits_cleanly(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _cmd: "/x/snakemake")

        def _timeout(*a, **k):
            raise subprocess.TimeoutExpired(cmd="snakemake", timeout=120)

        monkeypatch.setattr(epicc.subprocess, "run", _timeout)
        with pytest.raises(SystemExit):
            epicc.check_snakemake_version()


# ---------------------------------------------------------------------------
# Which subcommands the guard applies to
# ---------------------------------------------------------------------------

class TestGuardScope:
    """The guard must not block subcommands that never call snakemake.

    `perf` runs dev/profile_snakemake_log.py under sys.executable, and
    `clean` does filesystem work (only its --conda-envs branch needs
    snakemake). Both have to keep working on a machine with no epicc
    environment -- e.g. profiling a log fetched from a cluster.
    """

    def test_guarded_commands(self):
        assert epicc.SNAKEMAKE_COMMANDS == {"run", "validate", "output", "unlock"}

    @pytest.mark.parametrize("command", ["perf", "clean", "init-profile"])
    def test_unguarded_commands(self, command):
        assert command not in epicc.SNAKEMAKE_COMMANDS

    def test_output_is_guarded(self):
        # cmd_output builds `["snakemake", ...]`, so it does need the check.
        assert "output" in epicc.SNAKEMAKE_COMMANDS
