# Test Data Files

This directory contains mock data files used for unit testing the direct methylation (dmC) feature.

## Files

### Reference Files

- `sample.chrom.sizes` - Mock chromosome sizes file for Arabidopsis thaliana (ColCEN reference)

### Valid Test Files

- `sample_valid.bedmethyl` - Valid bedMethyl format file with 11 columns
  - Includes track header and comment line
  - Contains methylation data for multiple chromosomes
  - All coordinates and values are valid

### Invalid Test Files

- `sample_invalid_coords.bedmethyl` - Invalid coordinates (end position before start)
- `sample_invalid_percent.bedmethyl` - Invalid percent modified value (>100)

## BedMethyl Format

The bedMethyl format used by ONT modkit has 11 columns:

1. chrom - Chromosome name
2. start - Start position (0-based)
3. end - End position
4. name - Modified base code (e.g., "m" for 5mC)
5. score - Score value (0-1000)
6. strand - Strand (+ or -)
7. thickStart - Start of thick region (same as start)
8. thickEnd - End of thick region (same as end)
9. itemRgb - Color value (RGB format)
10. coverage - Number of valid reads (Nvalid)
11. percent_modified - Percent of reads showing modification (0-100)

## Usage in Tests

These files are used by:
- `tests/unit/test_validate_dmc_input.py` - For testing bedMethyl validation
- Integration tests that need realistic sample data

## Adding New Test Data

When adding new test data:
1. Use the same format as existing files
2. Document the purpose in this README
3. Keep files small (< 100 lines) for fast testing
4. Name files descriptively (e.g., `sample_invalid_<reason>.bedmethyl`)
