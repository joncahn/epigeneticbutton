# Integration Tests Quick Start Guide

## TL;DR

Run dmC (direct methylation) integration tests:

```bash
# From repository root
pytest tests/integration/test_dmc_dryrun.py -v
```

## What These Tests Do

Integration tests verify that the Snakemake pipeline can correctly build the DAG (Directed Acyclic Graph) for dmC (direct methylation) workflows **without executing any rules**. They use Snakemake's `--dry-run` mode.

## Requirements

```bash
# Install snakemake
conda install -c conda-forge -c bioconda snakemake

# OR
pip install snakemake
```

## Test Structure

```
tests/integration/
├── README.md                           # Full documentation
├── QUICKSTART.md                       # This file
├── test_dmc_dryrun.py                  # dmC integration tests (25 tests)
└── data/
    ├── test_samples_dmc.tsv           # Mock sample metadata
    └── test_config_dmc.yaml           # Test configuration
```

## Running Tests

### All dmC tests
```bash
pytest tests/integration/test_dmc_dryrun.py -v
```

### Specific test class
```bash
pytest tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow -v
```

### Specific test
```bash
pytest tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow::test_dmc_modbam_dryrun_succeeds -v
```

### With verbose output (see Snakemake output)
```bash
pytest tests/integration/test_dmc_dryrun.py -v -s
```

## Test Coverage

The 25 tests cover:

1. **dmC modBAM workflow** (7 tests)
   - Dry-run success
   - Correct rule triggering (get_modbam, align_modbam, modkit_pileup)
   - Exclusion of Bismark rules
   - All contexts (CG, CHG, CHH)

2. **bedMethyl workflow** (3 tests)
   - Pre-computed bedMethyl handling
   - Skipping of alignment steps
   - Correct rule selection

3. **DMR analysis** (2 tests)
   - DMRcaller for DMR calling
   - Merged replicate handling

4. **DAG structure** (3 tests)
   - DAG generation
   - Rule dependencies
   - dmC rule inclusion

5. **Wildcard resolution** (2 tests)
   - dmC sample wildcards
   - bedMethyl sample wildcards

6. **Error handling** (2 tests)
   - Invalid reference genome
   - Invalid context

7. **Pipeline targets** (2 tests)
   - all_mc rule
   - dmC output inclusion

8. **Context BEDs** (2 tests)
   - Generation of CG/CHG/CHH BEDs
   - Dependency verification

9. **Multiple replicates** (2 tests)
   - Replicate processing
   - Merged DMR analysis

## Manual Dry-Run Testing

Test individual targets manually:

```bash
# dmC modBAM sample
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw

# bedMethyl sample
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__root__bedMethyl__rep1__test_genome__CG.bw

# DMR analysis
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/DMRs/summary__mC__WT__leaf__dmC__test_genome__vs__mC__mutant__leaf__dmC__test_genome__DMRs.txt
```

## Generate DAG Visualization

```bash
snakemake --dag \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw \
    | dot -Tpng > dag_dmc.png
```

## Expected Results

### If Snakemake is Installed

```
tests/integration/test_dmc_dryrun.py::TestDmcDryRunBasic::test_snakemake_installed PASSED
tests/integration/test_dmc_dryrun.py::TestDmcDryRunBasic::test_config_file_exists PASSED
tests/integration/test_dmc_dryrun.py::TestDmcDryRunBasic::test_sample_file_exists PASSED
...
========================= 25 passed in X.XXs =========================
```

### If Snakemake is Not Installed

```
tests/integration/test_dmc_dryrun.py::TestDmcDryRunBasic::test_snakemake_installed SKIPPED
tests/integration/test_dmc_dryrun.py::TestDmcDryRunBasic::test_config_file_exists PASSED
tests/integration/test_dmc_dryrun.py::TestDmcDryRunBasic::test_sample_file_exists PASSED
...
========================= 2 passed, 23 skipped in X.XXs =========================
```

Tests automatically skip if Snakemake is not available.

## What Tests DON'T Require

- Actual reference genomes
- Real sequencing data (BAM/bedMethyl files)
- Conda environments to be built
- Rules to be executed
- Software tools (modkit, samtools, etc.)

## What Tests DO Require

- Snakemake installed and in PATH
- Python 3.8+
- pytest

## Troubleshooting

### "Snakemake not installed"
```bash
# Check installation
which snakemake
snakemake --version

# Install if missing
conda install -c conda-forge -c bioconda snakemake
```

### "Config file not found"
```bash
# Run from repository root
cd /path/to/epigeneticbutton
pytest tests/integration/test_dmc_dryrun.py -v
```

### Tests timeout
Increase timeout in test file (default 60s):
```python
timeout=120  # in subprocess.run() calls
```

## Adding New Tests

1. Add test samples to `data/test_samples_dmc.tsv`
2. Update config if needed: `data/test_config_dmc.yaml`
3. Add test methods to `test_dmc_dryrun.py`

Example:
```python
def test_new_feature(self, snakemake_available, repo_root, test_config):
    if not snakemake_available:
        pytest.skip("Snakemake not installed")

    target = "results/new_feature/output.txt"
    result = run_snakemake_dryrun(repo_root, test_config, target)
    assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
```

## CI/CD Integration

See `.github_workflow_example.yml` for GitHub Actions integration.

## More Information

- Full documentation: `README.md`
- Unit tests: `../unit/README.md`
- Pipeline documentation: `../../CLAUDE.md`
