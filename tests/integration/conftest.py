"""Shared fixtures and helpers for integration tests."""

import subprocess
import yaml
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(__file__).parent / "data"


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
    """Check if snakemake is available."""
    try:
        result = subprocess.run(
            ["snakemake", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
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
