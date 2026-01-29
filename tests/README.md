# EpigeneticButton Test Suite

Comprehensive test suite for the EpigeneticButton (EPICC) bioinformatics pipeline, with emphasis on the direct methylation (dmC) analysis feature for native base modification detection from long-read sequencing.

## Quick Start

```bash
# Install test dependencies
pip install -r tests/requirements-test.txt

# Run all unit tests
./tests/run_tests.sh

# Run with coverage report
./tests/run_tests.sh --cov

# Run specific tests
pytest tests/unit/test_validate_dmc_input.py -v
```

## Directory Structure

```
tests/
├── README.md                           # This file
├── conftest.py                         # Shared pytest fixtures
├── requirements-test.txt               # Testing dependencies
├── run_tests.sh                        # Convenient test runner script
│
├── unit/                               # Unit tests
│   ├── README.md                       # Unit test documentation
│   ├── test_validate_dmc_input.py     # Tests for dmC input validation
│   └── test_mC_helpers.py             # Tests for mC.smk helper functions
│
├── integration/                        # Integration tests (future)
│   └── README.md                       # Integration test documentation
│
└── data/                               # Mock data for testing
    ├── README.md                       # Data file documentation
    ├── sample.chrom.sizes              # Mock chromosome sizes
    ├── sample_valid.bedmethyl          # Valid bedMethyl example
    ├── sample_invalid_coords.bedmethyl # Invalid coordinates example
    └── sample_invalid_percent.bedmethyl # Invalid percent example
```

## Test Categories

### Unit Tests (`tests/unit/`)

Fast, isolated tests for individual functions and modules. No external dependencies required.

**Coverage:**
- Input validation for modBAM files (MM/ML tag detection, alignment validation)
- Input validation for bedMethyl files (format, coordinates, coverage values)
- Automatic format detection (modBAM vs bedMethyl)
- Sample name parsing and type detection
- dmC vs Bismark sample identification
- Parameter selection logic for different methylation assay types

**Run unit tests:**
```bash
pytest tests/unit/ -v
```

### Integration Tests (`tests/integration/`)

Tests for complete workflows and pipelines. May require external tools (samtools, modkit).

**Status:** Placeholder directory for future integration tests

## Test Files

### validate_dmc_input.py Tests

`tests/unit/test_validate_dmc_input.py` - 40+ tests covering:

- **Valid inputs:** 11-column bedMethyl, 10-column bedMethyl, gzipped files
- **Invalid inputs:** Wrong column counts, bad coordinates, invalid values
- **Edge cases:** Empty files, headers, chromosome validation
- **modBAM validation:** MM/ML tag checking, reference alignment
- **CLI behavior:** Argument parsing, exit codes

### mC.smk Helper Tests

`tests/unit/test_mC_helpers.py` - 30+ tests covering:

- **parse_sample_name:** Standard names, ChIP groups, TF names
- **is_dmc_sample:** dmC detection, Bismark vs dmC
- **get_dmc_input_type:** modBAM vs bedMethyl detection
- **parameters_for_mc:** Parameter set selection for all sample types
- **Integration scenarios:** Complete workflows, mixed sample types

## Running Tests

### Basic Usage

```bash
# All tests with verbose output
pytest tests/unit/ -v

# Single test file
pytest tests/unit/test_validate_dmc_input.py

# Single test class
pytest tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl

# Single test function
pytest tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl::test_valid_bedmethyl_11_columns
```

### With Coverage

```bash
# Generate coverage report
pytest tests/unit/ --cov=workflow/scripts --cov-report=html --cov-report=term

# View HTML report
firefox htmlcov/index.html
```

### Using the Test Runner Script

```bash
# Basic run
./tests/run_tests.sh

# With coverage
./tests/run_tests.sh --cov

# Verbose output
./tests/run_tests.sh -v

# Run specific tests
./tests/run_tests.sh -k test_validate_bedmethyl
```

### Advanced Options

```bash
# Stop on first failure
pytest tests/unit/ -x

# Run last failed tests only
pytest tests/unit/ --lf

# Run tests in parallel (requires pytest-xdist)
pytest tests/unit/ -n auto

# Show test durations
pytest tests/unit/ --durations=10

# Drop into debugger on failure
pytest tests/unit/ --pdb
```

## Test Fixtures

Shared fixtures are defined in `tests/conftest.py`:

- `temp_dir` - Temporary directory for test files
- `sample_chrom_sizes` - Mock chromosome sizes file
- `valid_bedmethyl_content` - Valid bedMethyl data (11 columns)
- `valid_bedmethyl_10col_content` - Valid bedMethyl data (10 columns)
- `bedmethyl_with_header` - bedMethyl with track header
- `mock_bam_header` - Mock BAM file header
- `mock_bam_read_with_mm_ml` - BAM read with methylation tags
- `mock_bam_read_without_mm` - BAM read without methylation tags
- `dmc_sample_names` - Dictionary of sample name examples

## Writing Tests

### Test Structure

Follow the Arrange-Act-Assert pattern:

```python
def test_example(self, fixture):
    # Arrange - set up test data
    input_data = prepare_input()

    # Act - call the function under test
    result = function_under_test(input_data)

    # Assert - verify the expected outcome
    assert result == expected_value
```

### Naming Conventions

- Test files: `test_<module>.py`
- Test classes: `Test<Feature>`
- Test functions: `test_<function>_<scenario>_<expected>`

Example: `test_validate_bedmethyl_invalid_coordinates`

### Test Categories with Markers

Mark tests by category:

```python
import pytest

@pytest.mark.unit
def test_fast_unit_test():
    pass

@pytest.mark.slow
def test_slow_integration():
    pass

@pytest.mark.requires_samtools
def test_with_external_tool():
    pass
```

Run specific categories:
```bash
pytest -m unit          # Only unit tests
pytest -m "not slow"    # Skip slow tests
```

## Continuous Integration

### GitHub Actions Example

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, '3.10']

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}

    - name: Install dependencies
      run: |
        pip install -r tests/requirements-test.txt

    - name: Run tests with coverage
      run: |
        pytest tests/unit/ --cov=workflow/scripts --cov-report=xml

    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml
```

## Test Coverage Goals

Target coverage percentages:

- **validate_dmc_input.py**: >90%
- **mC.smk helper functions**: >95%

Check current coverage:
```bash
pytest tests/unit/ --cov=workflow/scripts/validate_dmc_input --cov-report=term-missing
```

## Troubleshooting

### Import Errors

```bash
# Run from repository root
cd /grid/martienssen/home/eernst/src/epigeneticbutton

# Set PYTHONPATH if needed
PYTHONPATH=. pytest tests/unit/
```

### Missing Dependencies

```bash
# Install test dependencies
pip install -r tests/requirements-test.txt

# Verify pytest is installed
pytest --version
```

### Test Failures

1. Run with verbose output: `pytest -v`
2. Check specific test: `pytest tests/unit/test_file.py::test_name -v`
3. Use debugger: `pytest --pdb`
4. Check fixture setup in `conftest.py`

## Best Practices

1. **Write tests first** (TDD approach when adding new features)
2. **Keep tests fast** (unit tests should run in seconds)
3. **One assertion per test** (or closely related assertions)
4. **Use descriptive names** (test name should explain what is tested)
5. **Mock external dependencies** (files, network, subprocess)
6. **Clean up after tests** (use fixtures with cleanup/teardown)
7. **Test edge cases** (empty, null, boundary values, errors)
8. **Maintain high coverage** (aim for >90% on new code)

## Contributing

When adding new features:

1. Write tests for new functionality
2. Ensure all existing tests pass
3. Maintain or improve coverage
4. Update test documentation
5. Add new fixtures to `conftest.py` if needed

## Resources

- [pytest documentation](https://docs.pytest.org/)
- [pytest fixtures](https://docs.pytest.org/en/stable/fixture.html)
- [unittest.mock](https://docs.python.org/3/library/unittest.mock.html)
- [Coverage.py](https://coverage.readthedocs.io/)

## Support

For issues or questions about tests:

1. Check test output and error messages
2. Review relevant test documentation
3. Check existing issues in repository
4. Create new issue with test failure details
