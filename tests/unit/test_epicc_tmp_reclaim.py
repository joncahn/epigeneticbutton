"""Unit tests for `epicc clean --tmp` scratch-dir reclamation.

The shell.prefix in workflow/Snakefile creates a per-job scratch dir under
`<output_dir>/.tmp/${SLURM_JOB_ID:-local}.$$` and traps EXIT to remove it.
SIGKILL escapes the trap, so SLURM timeouts, OOM kills and node failures each
leak one, and nothing else ever reclaims them. Deleting scratch out from under
a *running* job would corrupt it, so the sweep has to be conservative: these
tests pin both safety gates. The `epicc` script has no .py extension, so we
load it via importlib (same pattern as test_epicc_placeholders.py).
"""

import importlib.util
import os
import subprocess
import sys
import time
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

_HOUR = 3600


def _make_dir(root, name, age_hours=0.0):
    d = root / name
    d.mkdir()
    (d / "spill.tmp").write_text("x")
    if age_hours:
        old = time.time() - (age_hours * _HOUR)
        os.utime(d, (old, old))
    return d


# ---------------------------------------------------------------------------
# parse_tmp_dir_name
# ---------------------------------------------------------------------------

class TestParseTmpDirName:
    def test_slurm_job_dir(self):
        assert epicc.parse_tmp_dir_name("3060216.12345") == ("3060216", "12345")

    def test_local_run_dir(self):
        # No SLURM job to consult; job component is None, pid is kept.
        assert epicc.parse_tmp_dir_name("local.12345") == (None, "12345")

    @pytest.mark.parametrize("name", [
        "notajob",              # no separator
        "3060216",              # no pid
        "3060216.abc",          # non-numeric pid
        "weird.name.here",      # trailing component not a pid
        "",
    ])
    def test_unrecognized_names(self, name):
        assert epicc.parse_tmp_dir_name(name) == (None, None)

    def test_dotted_job_component_is_not_guessed(self):
        # Anything we can't confidently read must be left alone.
        assert epicc.parse_tmp_dir_name("foo.bar.12345")[0] is None


# ---------------------------------------------------------------------------
# sweep_tmp_dirs -- safety gates
# ---------------------------------------------------------------------------

class TestSweepSafety:
    def test_active_slurm_job_is_never_reclaimed(self, tmp_path):
        _make_dir(tmp_path, "111.5", age_hours=99)   # old, but job is live
        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs={"111"}
        )
        assert reclaimable == []
        assert "still active" in skipped[0][1]

    def test_unknown_job_status_skips_slurm_dirs(self, tmp_path):
        # active_jobs=None means "cannot rule out a live job" -- must not be
        # read as "nothing is running".
        _make_dir(tmp_path, "222.5", age_hours=99)
        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs=None
        )
        assert reclaimable == []
        assert "status unknown" in skipped[0][1]

    def test_recent_dir_is_skipped(self, tmp_path):
        _make_dir(tmp_path, "333.5", age_hours=1)
        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs=set()
        )
        assert reclaimable == []
        assert "modified" in skipped[0][1]

    def test_array_task_id_matches_its_job(self, tmp_path):
        # squeue reports array tasks as "444_7"; the scratch dir is named from
        # SLURM_JOB_ID, so the leading component has to match.
        _make_dir(tmp_path, "444.5", age_hours=99)
        reclaimable, _ = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs={"444"}
        )
        assert reclaimable == []

    def test_unrecognized_name_is_left_alone(self, tmp_path):
        _make_dir(tmp_path, "something-else", age_hours=99)
        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs=set()
        )
        assert reclaimable == []
        assert "unrecognized" in skipped[0][1]

    def test_files_are_ignored(self, tmp_path):
        (tmp_path / "stray.txt").write_text("x")
        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs=set()
        )
        assert reclaimable == []
        assert skipped == []


# ---------------------------------------------------------------------------
# sweep_tmp_dirs -- the cases that should be reclaimed
# ---------------------------------------------------------------------------

class TestSweepReclaims:
    def test_dead_slurm_job_is_reclaimed(self, tmp_path):
        d = _make_dir(tmp_path, "555.5", age_hours=99)
        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs={"999"}
        )
        assert reclaimable == [d]
        assert skipped == []

    def test_old_local_dir_is_reclaimed(self, tmp_path):
        d = _make_dir(tmp_path, "local.4242", age_hours=99)
        reclaimable, _ = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs=set()
        )
        assert reclaimable == [d]

    def test_local_dir_ignores_unknown_job_status(self, tmp_path):
        # A local run has no SLURM job, so squeue being unavailable is
        # irrelevant -- only the age gate applies.
        d = _make_dir(tmp_path, "local.4242", age_hours=99)
        reclaimable, _ = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs=None
        )
        assert reclaimable == [d]

    def test_mixed_set_is_partitioned(self, tmp_path):
        live = _make_dir(tmp_path, "700.1", age_hours=99)
        dead = _make_dir(tmp_path, "701.1", age_hours=99)
        recent = _make_dir(tmp_path, "702.1", age_hours=0.5)
        local_old = _make_dir(tmp_path, "local.9", age_hours=99)

        reclaimable, skipped = epicc.sweep_tmp_dirs(
            tmp_path, min_age_hours=6, active_jobs={"700"}
        )
        assert set(reclaimable) == {dead, local_old}
        assert {d for d, _ in skipped} == {live, recent}

    def test_empty_root(self, tmp_path):
        assert epicc.sweep_tmp_dirs(tmp_path, 6, set()) == ([], [])


# ---------------------------------------------------------------------------
# active_slurm_job_ids
# ---------------------------------------------------------------------------

class _FakeProc:
    def __init__(self, stdout="", returncode=0):
        self.stdout = stdout
        self.stderr = ""
        self.returncode = returncode


class TestActiveSlurmJobIds:
    def test_none_when_squeue_absent(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _c: None)
        assert epicc.active_slurm_job_ids() is None

    def test_parses_ids(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _c: "/usr/bin/squeue")
        monkeypatch.setattr(
            epicc.subprocess, "run",
            lambda *a, **k: _FakeProc(stdout="111\n222\n333\n"),
        )
        assert epicc.active_slurm_job_ids() == {"111", "222", "333"}

    def test_array_tasks_reduce_to_job_id(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _c: "/usr/bin/squeue")
        monkeypatch.setattr(
            epicc.subprocess, "run",
            lambda *a, **k: _FakeProc(stdout="444_1\n444_2\n555\n"),
        )
        assert epicc.active_slurm_job_ids() == {"444", "555"}

    def test_empty_queue_is_empty_set_not_none(self, monkeypatch):
        # An empty queue is a real answer; None means "unknowable".
        monkeypatch.setattr(epicc.shutil, "which", lambda _c: "/usr/bin/squeue")
        monkeypatch.setattr(
            epicc.subprocess, "run", lambda *a, **k: _FakeProc(stdout="")
        )
        assert epicc.active_slurm_job_ids() == set()

    def test_none_on_squeue_failure(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _c: "/usr/bin/squeue")
        monkeypatch.setattr(
            epicc.subprocess, "run",
            lambda *a, **k: _FakeProc(stdout="", returncode=1),
        )
        assert epicc.active_slurm_job_ids() is None

    def test_none_on_timeout(self, monkeypatch):
        monkeypatch.setattr(epicc.shutil, "which", lambda _c: "/usr/bin/squeue")

        def _timeout(*a, **k):
            raise subprocess.TimeoutExpired(cmd="squeue", timeout=60)

        monkeypatch.setattr(epicc.subprocess, "run", _timeout)
        assert epicc.active_slurm_job_ids() is None
