# Integration Tests for EpigeneticButton

This directory contains integration tests for the EpigeneticButton Snakemake pipeline. These tests verify that the pipeline rules and workflows are correctly structured and can build valid DAGs without executing the actual analysis.

## Overview

Integration tests use Snakemake's `--dry-run` mode to:
- Validate that the DAG can be constructed
- Verify wildcard resolution works correctly
- Ensure rule dependencies are properly defined
- Test that sample-specific workflows are triggered correctly
- Confirm that incompatible rules are not triggered

These tests do **not** require:
- Actual reference genomes
- Real sequencing data
- Conda environments to be built
- Rules to be executed

## Test Organization

### dmC (Direct Methylation) Integration Tests (`test_dmc_dryrun.py`)

Tests for direct methylation sequencing workflow (e.g., from Oxford Nanopore Technologies).

**Test Classes:**

1. `TestDmcDryRunBasic` - Basic setup validation
   - Snakemake installation check
   - Config and sample file existence

2. `TestDmcModBAMWorkflow` - dmC modBAM input workflow
   - Dry-run success for modBAM samples
   - Correct dmC rules are triggered (get_modbam, align_modbam, modkit_pileup, etc.)
   - Bismark rules are NOT triggered
   - All methylation contexts (CG, CHG, CHH) are generated

3. `TestBedMethylWorkflow` - Pre-computed bedMethyl input workflow
   - Dry-run success for bedMethyl samples
   - Correct bedMethyl rules are triggered (get_bedmethyl, copy_bedmethyl_for_pileup)
   - Alignment rules are skipped (no need to align pre-computed data)

4. `TestDmcDMRWorkflow` - DMR analysis for dmC samples
   - Dry-run success for DMR analysis
   - DMRcaller is used by default for DMR calling

5. `TestDAGStructure` - DAG generation and structure
   - DAG can be generated successfully
   - DAG contains expected dmC rules
   - Rule dependencies are correct

6. `TestWildcardResolution` - Wildcard resolution
   - dmC sample wildcards resolve correctly
   - bedMethyl sample wildcards resolve correctly

7. `TestErrorHandling` - Error cases
   - Missing reference genome fails gracefully
   - Invalid methylation context fails

8. `TestAllMCTarget` - Pipeline completion target
   - all_mc rule works with dmC samples
   - dmC-specific outputs are included

9. `TestContextBedGeneration` - Context BED file generation
   - CG/CHG/CHH context BED files are generated
   - Context BED files are dependencies for splitting

10. `TestMultipleReplicates` - Multiple replicate handling
    - Multiple dmC replicates are processed
    - Merged replicates work for DMR calling

## Test Data

Test data files are located in `tests/integration/data/`:

### `test_samples_dmc.tsv`
Mock sample metadata file with:
- 2 dmC modBAM samples (WT leaf, rep1 and rep2)
- 1 bedMethyl sample (WT root, rep1)
- 2 dmC modBAM samples for mutant (leaf, rep1 and rep2)

Format matches the pipeline's 9-column TSV format:
```
data_type  line  tissue  sample_type  replicate  seq_id  fastq_path  paired  ref_genome
```

### `test_config_dmc.yaml`
Minimal test configuration that:
- Points to test sample file
- Defines mock reference genome (test_genome)
- Sets dmC methylation parameters
- Uses minimal resource allocations
- Enables full analysis with DMR calling

## Running the Tests

### Prerequisites

1. Snakemake must be installed and in PATH:
   ```bash
   conda install -c conda-forge -c bioconda snakemake
   # OR
   pip install snakemake
   ```

2. Pytest should be installed:
   ```bash
   pip install pytest
   ```

### Run All Integration Tests

From the repository root:

```bash
# Run all integration tests
pytest tests/integration/ -v

# Run only dmC dry-run tests
pytest tests/integration/test_dmc_dryrun.py -v

# Run specific test class
pytest tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow -v

# Run specific test
pytest tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow::test_dmc_modbam_dryrun_succeeds -v
```

### Run with Integration Test Marker

Integration tests are marked with `@pytest.mark.integration`:

```bash
# Run only integration tests (if other test types exist)
pytest -m integration -v

# Run all tests except integration tests
pytest -m "not integration" -v
```

### Verbose Output

For debugging, use verbose mode to see Snakemake output:

```bash
pytest tests/integration/test_dmc_dryrun.py -v -s
```

### Skip if Snakemake Not Installed

Tests will automatically skip if Snakemake is not available:

```bash
pytest tests/integration/test_dmc_dryrun.py -v
# Output: SKIPPED [1] Snakemake not installed or not in PATH
```

## Manual Dry-Run Testing

You can also manually test the dry-run:

```bash
# Test dmC modBAM sample bigwig generation
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw

# Test bedMethyl sample
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__root__bedMethyl__rep1__test_genome__CG.bw

# Test DMR analysis
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/DMRs/summary__mC__WT__leaf__dmC__test_genome__vs__mC__mutant__leaf__dmC__test_genome__DMRs.txt

# Generate DAG visualization
snakemake --dag \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw \
    | dot -Tpng > dag_dmc.png
```

## What the Tests Verify

### dmC Workflow Correctness

1. **Rule Selection**
   - dmC samples trigger dmC-specific rules (modkit, etc.)
   - dmC samples do NOT trigger Bismark rules
   - bedMethyl samples skip alignment steps

2. **Wildcard Resolution**
   - Sample names are parsed correctly
   - File paths are constructed with correct wildcards
   - All methylation contexts (CG, CHG, CHH) are handled

3. **Dependencies**
   - Context BED files are generated before splitting
   - Alignment happens before pileup
   - Pileup happens before context splitting
   - Context splitting happens before bigwig generation

4. **DMR Analysis**
   - dmC samples use DMRcaller for DMR calling by default
   - Merged replicates are used correctly
   - Sample comparisons are properly structured

### Error Handling

Tests verify that invalid inputs fail appropriately:
- Non-existent reference genomes
- Invalid methylation contexts
- Missing sample files

## Interpreting Test Results

### Success

```
tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow::test_dmc_modbam_dryrun_succeeds PASSED
```

The Snakemake DAG was successfully built for the target, indicating correct rule structure.

### Failure

```
tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow::test_dmc_modbam_includes_expected_rules FAILED
AssertionError: Expected dmC rule 'modkit_pileup' not found in dry-run output
```

Indicates a problem with rule triggering or naming. Check:
1. Rule name in `workflow/rules/mC.smk`
2. Input functions for rule selection
3. Wildcard constraints

### Skipped

```
tests/integration/test_dmc_dryrun.py::TestDmcModBAMWorkflow::test_dmc_modbam_dryrun_succeeds SKIPPED
```

Test was skipped (usually because Snakemake is not installed). This is expected in environments without Snakemake.

## Adding New Integration Tests

To add new integration tests:

1. **Add test samples** to `test_samples_dmc.tsv` or create a new sample file

2. **Update test config** if needed in `test_config_dmc.yaml` or create a new config

3. **Add test methods** to existing test classes or create new classes:

```python
class TestNewFeature:
    """Test new feature dry-run."""

    @pytest.fixture
    def feature_target(self):
        """Return target for new feature."""
        return "results/new_feature/output.txt"

    def test_feature_dryrun(self, snakemake_available, repo_root, test_config, feature_target):
        """Test that dry-run succeeds for new feature."""
        if not snakemake_available:
            pytest.skip("Snakemake not installed")

        result = run_snakemake_dryrun(repo_root, test_config, feature_target)
        assert result.returncode == 0, f"Dry-run failed: {result.stderr}"
```

4. **Test both success and failure cases** to ensure robust validation

## Troubleshooting

### "Snakemake not installed" but it is installed

Ensure Snakemake is in your PATH:
```bash
which snakemake
snakemake --version
```

If using conda, activate the appropriate environment before running tests.

### Tests timeout

Increase the timeout in test functions (default is 60 seconds):
```python
result = subprocess.run(
    cmd,
    timeout=120  # Increase to 120 seconds
)
```

### "Config file not found"

Verify you're running pytest from the repository root:
```bash
cd /path/to/epigeneticbutton
pytest tests/integration/test_dmc_dryrun.py -v
```

### Dry-run fails with "AmbiguousRuleException"

This indicates overlapping rule outputs. Check:
1. Rule output patterns in `workflow/rules/mC.smk`
2. Wildcard constraints
3. Input functions that select rules

## CI/CD Integration

These tests can be integrated into GitHub Actions or other CI systems:

```yaml
name: Integration Tests

on: [push, pull_request]

jobs:
  integration-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: conda-incubator/setup-miniconda@v2
        with:
          python-version: 3.9
          channels: conda-forge,bioconda
      - name: Install dependencies
        run: |
          conda install snakemake
          pip install pytest
      - name: Run integration tests
        run: pytest tests/integration/ -v
```

## Additional Resources

- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [Pytest Documentation](https://docs.pytest.org/)
- [ONT Modkit Documentation](https://nanoporetech.github.io/modkit/)
- [EpigeneticButton Pipeline Documentation](../../CLAUDE.md)
