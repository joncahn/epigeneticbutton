"""Shared fixtures and helpers for integration tests."""

import csv
import os
import re
import shutil
import subprocess
import tempfile
import yaml
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(__file__).parent / "data"

# Sample sheets that drive the per-assay dry-run tests. Their Read_files
# columns reference paths under tests/integration/data/mock_inputs/ that
# are not committed (binary stubs would bloat the repo); the autouse
# fixture below stages empty placeholders before any integration tests run.
_DRYRUN_SAMPLE_SHEETS = (
    "test_samples_ATAC.tsv",
    "test_samples_ChIP.tsv",
    "test_samples_CUT.tsv",
    "test_samples_RNA.tsv",
    "test_samples_sRNA.tsv",
    "test_samples_mC.tsv",
)
_SRA_RE = re.compile(r"^[SDE]RR\d+$")


def _local_paths_from_sheet(sheet_path):
    """Yield local-filesystem entries from a sample sheet's Read_files column."""
    with open(sheet_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rf = (row.get("Read_files") or "").strip()
            if not rf:
                continue
            for component in rf.split("+"):
                for token in component.split(","):
                    t = token.strip()
                    if not t or t.startswith("http://") or t.startswith("https://"):
                        continue
                    if _SRA_RE.match(t):
                        continue
                    yield t


def _stage_mock_inputs():
    """Create empty placeholders for every Read_files path that the
    dry-run sample sheets reference but which isn't already on disk.

    Snakemake's dry-run never reads the file contents, and our sample-sheet
    validation only checks for existence, so empty stubs are sufficient.
    """
    for sheet_name in _DRYRUN_SAMPLE_SHEETS:
        sheet_path = DATA_DIR / sheet_name
        if not sheet_path.exists():
            continue
        for rel in _local_paths_from_sheet(sheet_path):
            target = (REPO_ROOT / rel).resolve()
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.touch()


@pytest.fixture(scope="session", autouse=True)
def stage_mock_inputs():
    """Auto-stage missing mock Read_files placeholders for integration tests."""
    _stage_mock_inputs()
    yield


def load_output_dir(options_filename):
    """Load output_dir from a test options YAML file."""
    options_path = DATA_DIR / options_filename
    with open(options_path) as f:
        return yaml.safe_load(f).get("output_dir", "results")


def load_genome_dir(options_filename):
    """Load genome_dir from a test options YAML file."""
    options_path = DATA_DIR / options_filename
    with open(options_path) as f:
        return yaml.safe_load(f).get("genome_dir", "genomes")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def repo_root():
    """Get the repository root directory."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def snakemake_available():
    """Check if snakemake is available.

    Snakemake's startup is import-heavy; when its conda env lives on an
    NFS-mounted shared filesystem (typical on HPC), a cold first
    invocation under load can take well over a minute. The fixture is
    session-scoped, so a single timeout silently skips every dryrun
    test in the session — better to wait. Cheap shortcut: shutil.which
    confirms the binary is present without invoking it; only fall back
    to the subprocess probe if which fails.
    """
    if shutil.which("snakemake") is None:
        return False
    try:
        result = subprocess.run(
            ["snakemake", "--version"],
            capture_output=True,
            text=True,
            timeout=300,
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


# ---------------------------------------------------------------------------
# Snakemake helpers
# ---------------------------------------------------------------------------
#
# Dry-run and DAG-build are pure *planning* operations, but Snakemake keeps
# its persistence (.snakemake/) under the working directory (the repo root),
# shared across every run. A real run that is interrupted leaves stale
# .snakemake/incomplete/ markers there; if one matches an output the planned
# DAG would produce, the planning step fails:
#   - dry-run raises IncompleteFilesException;
#   - --dag additionally trips a Snakemake checkpoint-postprocess crash
#     (`'NoneType' object has no attribute 'edit_notebook'`, because no-exec
#     modes never initialise workflow.execution_settings).
# Neither failure reflects a problem with the workflow under test. We make
# both helpers hermetic by redirecting output_dir/genome_dir to an ephemeral
# temp tree: the planned DAG then targets fresh paths that no stale marker
# can match, while its rule structure (what the tests assert on) is identical.
# The requested target is rewritten into the temp output_dir so it resolves.


def _config_dirs(options_file):
    """Return (output_dir, genome_dir) declared in a Snakemake options YAML."""
    with open(options_file) as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg.get("output_dir", "results"), cfg.get("genome_dir", "genomes")


def _retarget(target, real_out, tmp_out):
    """Rewrite a path target's output_dir prefix to the hermetic temp dir.

    Phony targets (rule names like ``map_only``) and any target not under the
    configured output_dir pass through unchanged.
    """
    if not target:
        return target
    prefix = real_out.rstrip("/")
    if target == prefix or target.startswith(prefix + "/"):
        return tmp_out + target[len(prefix):]
    return target


def _split_config(extra_args):
    """Split extra_args into (config KEY=VALUE assignments, remaining args).

    Snakemake honors a single ``--config`` group; a second one would override
    the first. Since we inject our own output_dir/genome_dir assignments, any
    caller-supplied ``--config`` entries must be merged into one group rather
    than passed as a competing flag.
    """
    config_pairs = []
    other = []
    it = iter(extra_args or [])
    for arg in it:
        if arg == "--config":
            for nxt in it:  # collect KEY=VALUE tokens until the next flag
                if nxt.startswith("-"):
                    other.append(nxt)
                    break
                config_pairs.append(nxt)
        else:
            other.append(arg)
    return config_pairs, other


def _run_planning(mode_args, repo_root, options_file, target,
                  extra_args=None, timeout=120):
    """Run a no-exec Snakemake planning command in a hermetic output tree."""
    real_out, _ = _config_dirs(options_file)
    config_pairs, other = _split_config(extra_args)
    with tempfile.TemporaryDirectory(prefix="epicc_plan_") as tmp:
        tmp_out = os.path.join(tmp, "results")
        tmp_genomes = os.path.join(tmp, "genomes")
        cmd = [
            "snakemake",
            *mode_args,
            # No-exec planning needs no working-dir lock; --nolock avoids
            # spurious LockExceptions when many planning runs hit the shared
            # repo-root .snakemake/ back-to-back (e.g. NFS lock-release lag).
            "--nolock",
            "--configfile", options_file,
            "--config", f"output_dir={tmp_out}", f"genome_dir={tmp_genomes}",
            *config_pairs,
            "--cores", "1",
        ]

        cmd.extend(other)

        # Use -- separator to prevent Snakemake 9 from interpreting targets as options
        if target:
            cmd.append("--")
            cmd.append(_retarget(target, real_out, tmp_out))

        return subprocess.run(
            cmd,
            cwd=str(repo_root),
            capture_output=True,
            text=True,
            timeout=timeout,
        )


def run_snakemake_dryrun(repo_root, options_file, target=None,
                         extra_args=None, timeout=120):
    """Run snakemake in dry-run mode with the given config."""
    return _run_planning(
        ["--dry-run", "--quiet", "progress"],
        repo_root, options_file, target,
        extra_args=extra_args, timeout=timeout,
    )


def run_snakemake_dag(repo_root, options_file, target=None, timeout=120):
    """Generate the Snakemake DAG."""
    return _run_planning(
        ["--dag"],
        repo_root, options_file, target,
        timeout=timeout,
    )
