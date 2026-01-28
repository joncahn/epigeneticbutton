# Test Suite Quick Start Guide

Get up and running with the dmC (direct methylation) test suite in 5 minutes.

## Prerequisites

- Python 3.8 or higher
- pip package manager
- Access to the repository

## Installation (1 minute)

```bash
# Navigate to repository root
cd /grid/martienssen/home/eernst/src/epigeneticbutton

# Install test dependencies
pip install -r tests/requirements-test.txt
```

This installs:
- pytest (testing framework)
- pytest-cov (coverage reporting)
- pytest-mock (mocking utilities)
- pytest-xdist (parallel execution)
- pytest-sugar (better output)

## Run Tests (30 seconds)

### Option 1: Use the test runner script (recommended)

```bash
./tests/run_tests.sh
```

### Option 2: Run pytest directly

```bash
pytest tests/unit/ -v
```

## Expected Output

```
============================= test session starts ==============================
platform linux -- Python 3.9.7, pytest-7.4.0, pluggy-1.2.0
rootdir: /grid/martienssen/home/eernst/src/epigeneticbutton
configfile: pytest.ini
testpaths: tests
plugins: cov-4.1.0, mock-3.11.1, xdist-3.3.1
collected 72 items

tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl::test_valid_bedmethyl_11_columns PASSED [  1%]
tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl::test_valid_bedmethyl_10_columns PASSED [  2%]
...
tests/unit/test_mC_helpers.py::TestIntegrationScenarios::test_mixed_sample_types PASSED [100%]

============================== 72 passed in 2.45s ===============================
```

## Common Commands

### Run with coverage report

```bash
./tests/run_tests.sh --cov
```

Output shows coverage percentage for each file:
```
Name                                    Stmts   Miss  Cover   Missing
---------------------------------------------------------------------
workflow/scripts/validate_dmc_input.py    127      5    96%   45-47, 89, 126
---------------------------------------------------------------------
TOTAL                                     127      5    96%

Coverage HTML report: htmlcov/index.html
```

### Run specific test file

```bash
pytest tests/unit/test_validate_dmc_input.py -v
```

### Run specific test

```bash
pytest tests/unit/test_validate_dmc_input.py::TestValidateBedMethyl::test_valid_bedmethyl_11_columns -v
```

### Run tests matching a pattern

```bash
pytest tests/unit/ -k "bedmethyl" -v     # Only bedmethyl tests
pytest tests/unit/ -k "dmc_sample" -v    # Only dmc_sample tests
```

## Understanding Test Results

### ✓ PASSED - Test succeeded
```
test_validate_dmc_input.py::test_valid_bedmethyl_11_columns PASSED
```

### ✗ FAILED - Test failed
```
test_validate_dmc_input.py::test_valid_bedmethyl_11_columns FAILED

AssertionError: assert False is True
```

### s SKIPPED - Test was skipped
```
test_validate_dmc_input.py::test_requires_samtools SKIPPED
```

## Troubleshooting

### "pytest: command not found"

```bash
# Make sure you've installed the test dependencies
pip install -r tests/requirements-test.txt

# Verify installation
pytest --version
```

### Import errors

```bash
# Make sure you're in the repository root
cd /grid/martienssen/home/eernst/src/epigeneticbutton

# Set PYTHONPATH if needed
PYTHONPATH=. pytest tests/unit/
```

### "No module named 'validate_dmc_input'"

The test file imports from `workflow/scripts/`. Make sure you're running from the repository root.

## Next Steps

### View detailed documentation

```bash
# Main test suite documentation
cat tests/README.md

# Unit test documentation
cat tests/unit/README.md

# Implementation summary
cat tests/TEST_SUITE_SUMMARY.md
```

### View HTML coverage report

```bash
# Generate coverage
./tests/run_tests.sh --cov

# Open in browser
firefox htmlcov/index.html
```

### Run tests in parallel (faster)

```bash
pytest tests/unit/ -n auto
```

## Test Suite Structure

```
tests/
├── conftest.py              # Shared fixtures
├── requirements-test.txt    # Dependencies
├── run_tests.sh            # Test runner script
│
├── unit/                   # Unit tests (72 tests)
│   ├── test_validate_dmc_input.py  # 42 tests
│   └── test_mC_helpers.py          # 30 tests
│
└── data/                   # Mock data
    ├── sample.chrom.sizes
    └── sample_valid.bedmethyl
```

## What's Being Tested

### validate_dmc_input.py (42 tests)
- bedMethyl file format validation
- modBAM file validation
- Chromosome reference checking
- Error handling and edge cases

### mC.smk helpers (30 tests)
- `is_dmc_sample()` - dmC sample detection
- `get_dmc_input_type()` - Input type detection
- `parameters_for_mc()` - Parameter selection
- Sample name parsing

## Writing Your First Test

1. Open an existing test file
2. Copy a test function
3. Modify it for your use case
4. Run it!

Example:
```python
def test_my_new_feature(self, temp_dir):
    """Test description here."""
    # Arrange
    input_file = temp_dir / "test.bed"
    input_file.write_text("Chr1\t100\t200\n")

    # Act
    result = my_function(str(input_file))

    # Assert
    assert result is True
```

## Getting Help

1. Check error messages carefully
2. Read the documentation in `tests/README.md`
3. Look at existing tests for examples
4. Run with `-v` for verbose output
5. Use `--pdb` to drop into debugger on failures

## Tips

- Tests should run fast (< 5 seconds total)
- Use fixtures for reusable test data
- Mock external dependencies (subprocess, file I/O)
- Write descriptive test names
- One logical assertion per test
- Keep tests independent

## Success Criteria

You've successfully set up the test suite when:

- ✓ All 72 tests pass
- ✓ Coverage report is generated
- ✓ Tests complete in < 5 seconds
- ✓ No import errors
- ✓ HTML coverage report opens

## Resources

- [pytest documentation](https://docs.pytest.org/)
- [Repository test documentation](tests/README.md)
- [Test implementation summary](tests/TEST_SUITE_SUMMARY.md)

---

**Questions?** Check the main README or review the test documentation.

**Happy testing!** 🧪
