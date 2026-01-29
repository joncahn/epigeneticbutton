# Files Created for dmC Integration Tests

This document lists all files created for the dmC (direct methylation) dry-run integration tests.

## Test Files

### 1. Test Sample Metadata
**File:** `tests/integration/data/test_samples_dmc.tsv`

Mock sample metadata file with 5 dmC samples:
- 2 dmC modBAM samples (WT leaf, rep1 & rep2)
- 1 bedMethyl sample (WT root, rep1)
- 2 dmC modBAM samples (mutant leaf, rep1 & rep2)

Format: 9-column TSV matching pipeline requirements
- Columns: data_type, line, tissue, sample_type, replicate, seq_id, fastq_path, paired, ref_genome
- Uses mock paths (no actual files needed)
- Uses test reference genome "test_genome"

### 2. Test Configuration
**File:** `tests/integration/data/test_config_dmc.yaml`

Minimal configuration file for dry-run testing:
- Points to test sample file
- Defines mock reference genome (test_genome)
- Includes dmC methylation parameters
- Sets minimal resource allocations
- Enables full analysis with DMR calling
- No actual files required (uses mock paths)

### 3. Integration Test Suite
**File:** `tests/integration/test_dmc_dryrun.py`

Comprehensive pytest test suite with 25 tests organized into 10 test classes:

**Test Classes:**
1. `TestDmcDryRunBasic` (3 tests) - Basic setup validation
2. `TestDmcModBAMWorkflow` (4 tests) - dmC modBAM workflow
3. `TestBedMethylWorkflow` (3 tests) - bedMethyl input workflow
4. `TestDmcDMRWorkflow` (2 tests) - DMR analysis
5. `TestDAGStructure` (3 tests) - DAG generation and structure
6. `TestWildcardResolution` (2 tests) - Wildcard resolution
7. `TestErrorHandling` (2 tests) - Error cases
8. `TestAllMCTarget` (2 tests) - Pipeline completion target
9. `TestContextBedGeneration` (2 tests) - Context BED generation
10. `TestMultipleReplicates` (2 tests) - Multiple replicate handling

**Key Features:**
- Uses Snakemake `--dry-run` mode
- Automatically skips if Snakemake not installed
- Verifies rule selection (dmC rules included, Bismark rules excluded)
- Tests wildcard resolution
- Validates DAG structure and dependencies
- Checks error handling for invalid inputs

## Documentation Files

### 4. Comprehensive README
**File:** `tests/integration/README.md`

Full documentation covering:
- Overview of integration testing approach
- Detailed description of all test classes
- Test data file descriptions
- Running tests (various options)
- Manual dry-run testing commands
- DAG visualization
- What tests verify
- Interpreting test results
- Adding new tests
- Troubleshooting guide
- CI/CD integration
- Additional resources

### 5. Quick Start Guide
**File:** `tests/integration/QUICKSTART.md`

Quick reference for:
- TL;DR commands
- Requirements
- Test structure
- Running tests (examples)
- Test coverage summary
- Manual dry-run commands
- Expected results
- Troubleshooting
- Adding new tests

### 6. GitHub Actions Example
**File:** `tests/integration/.github_workflow_example.yml`

Example CI/CD workflow showing:
- Matrix testing (multiple Python/Snakemake versions)
- Conda environment setup
- Integration test execution
- Coverage reporting
- DAG visualization artifact upload
- Test result summary

### 7. This File
**File:** `tests/integration/FILES_CREATED.md`

Summary of all created files with descriptions.

## File Tree

```
tests/integration/
├── README.md                           # Full documentation (detailed)
├── QUICKSTART.md                       # Quick start guide (concise)
├── FILES_CREATED.md                    # This file (summary)
├── .github_workflow_example.yml       # CI/CD example
├── test_dmc_dryrun.py                 # Test suite (25 tests)
└── data/
    ├── test_samples_dmc.tsv           # Mock sample metadata (5 samples)
    └── test_config_dmc.yaml           # Test configuration
```

## Usage

### Quick Test
```bash
pytest tests/integration/test_dmc_dryrun.py -v
```

### Collect Tests Only
```bash
pytest tests/integration/test_dmc_dryrun.py --collect-only
```

### Manual Dry-Run
```bash
snakemake --dry-run \
    --configfile tests/integration/data/test_config_dmc.yaml \
    results/mC/tracks/mC__WT__leaf__dmC__rep1__test_genome__CG.bw
```

## File Sizes

```bash
# Check file sizes
du -h tests/integration/data/test_samples_dmc.tsv    # ~400 bytes
du -h tests/integration/data/test_config_dmc.yaml    # ~3.5 KB
du -h tests/integration/test_dmc_dryrun.py           # ~21 KB
du -h tests/integration/README.md                     # ~19 KB
du -h tests/integration/QUICKSTART.md                 # ~5 KB
```

## Dependencies

**Required:**
- Python 3.8+
- pytest
- Snakemake (for tests to run; will skip if not installed)

**Not Required:**
- Actual reference genomes
- Real sequencing data
- Conda environments
- Bioinformatics tools (modkit, samtools, etc.)

## Testing Coverage

**What is Tested:**
- DAG construction and validation
- Rule selection based on sample type
- Wildcard resolution
- Rule dependencies
- dmC-specific workflow paths
- bedMethyl-specific workflow paths
- DMR analysis workflow
- Error handling for invalid inputs
- Multiple replicate handling
- Methylation context (CG/CHG/CHH) handling

**What is NOT Tested:**
- Actual rule execution
- Tool functionality (modkit, samtools, etc.)
- Data processing correctness
- Output file formats
- Performance/resource usage

## Integration with Existing Tests

These integration tests complement the existing unit tests:

**Unit Tests** (`tests/unit/`):
- Test individual functions
- Test helper functions (parse_sample_name, is_dmc_sample, etc.)
- Test validation scripts
- Fast execution
- No external dependencies

**Integration Tests** (`tests/integration/`):
- Test workflow assembly
- Test rule interactions
- Test Snakemake DAG construction
- Medium execution time
- Requires Snakemake

## Maintenance

**When to Update These Tests:**

1. **New dmC rules added** - Add tests to verify rule triggering
2. **Sample naming changes** - Update test samples and wildcard tests
3. **Config structure changes** - Update test config
4. **New output targets** - Add tests for new targets
5. **Error handling changes** - Update error tests

## Related Files

**Existing Files (Not Modified):**
- `tests/conftest.py` - Shared pytest fixtures (already includes dmC fixtures)
- `pytest.ini` - Pytest configuration (already has integration marker)
- `tests/unit/test_mC_helpers.py` - Unit tests for mC helper functions
- `tests/unit/test_validate_dmc_input.py` - Unit tests for dmC validation

**Pipeline Files (Tested):**
- `workflow/Snakefile` - Main pipeline orchestrator
- `workflow/rules/mC.smk` - Methylation analysis rules (including dmC)
- `config/config.yaml` - Example configuration

## License

These test files follow the same license as the main EpigeneticButton repository.

## Contact

For questions or issues with these tests, please open an issue in the repository.
