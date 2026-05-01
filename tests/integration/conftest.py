"""Shared fixtures and helpers for integration tests."""

import csv
import re
import shutil
import subprocess
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
    "test_samples_dmc.tsv",
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

def run_snakemake_dryrun(repo_root, options_file, target=None,
                         extra_args=None, timeout=120):
    """Run snakemake in dry-run mode with the given config."""
    cmd = [
        "snakemake",
        "--dry-run",
        "--configfile", options_file,
        "--cores", "1",
        "--quiet", "progress",
    ]

    if extra_args:
        cmd.extend(extra_args)

    # Use -- separator to prevent Snakemake 9 from interpreting targets as options
    if target:
        cmd.append("--")
        cmd.append(target)

    return subprocess.run(
        cmd,
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def run_snakemake_dag(repo_root, options_file, target=None, timeout=120):
    """Generate the Snakemake DAG."""
    cmd = [
        "snakemake",
        "--dag",
        "--configfile", options_file,
        "--cores", "1",
    ]

    # Use -- separator to prevent Snakemake 9 from interpreting targets as options
    if target:
        cmd.append("--")
        cmd.append(target)

    return subprocess.run(
        cmd,
        cwd=str(repo_root),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
