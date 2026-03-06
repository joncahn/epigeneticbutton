# EPICC Test Suite

## Quick Start

```bash
conda activate epicc

# Unit tests (fast, <1s)
pytest tests/unit/ -v

# Integration dry-run (S. pombe, ~43s)
pytest tests/integration/test_pombe_dryrun.py -v

# Full pombe validation (dry-run + pipeline execution + output checks)
scripts/validate_pombe.sh --all
```

See `dev/docs/README.md` for environment setup (pytest, samtools, and other dev dependencies).

## Directory Structure

```
tests/
├── README.md                          # This file
├── conftest.py                        # Shared pytest fixtures
│
├── unit/                              # Fast unit tests (110+, <1s)
│   ├── test_sample_sheet.py           #   Sample sheet parsing and validation
│   ├── test_rule_commands.py          #   Samtools pipelines on synthetic SAM data
│   ├── test_validate_dmc_input.py     #   dmC input validation (modBAM, bedMethyl)
│   └── test_mC_helpers.py             #   mC.smk helper functions
│
├── integration/                       # Snakemake dry-run and post-run tests
│   ├── test_pombe_dryrun.py           #   S. pombe DAG validation (31 tests, ~43s)
│   ├── test_pombe_postrun.py          #   Post-run output checks (29 tests)
│   ├── test_dmc_dryrun.py             #   dmC rule routing (modBAM vs bedMethyl)
│   └── data/
│       ├── test_samples_pombe.tsv     #   17-sample S. pombe (4 assays)
│       ├── test_options_pombe.yaml
│       ├── Spombe/                    #   S. pombe reference genome
│       ├── test_samples_chr21.tsv     #   33-sample human chr21 (7 assays)
│       ├── test_options_chr21.yaml
│       ├── human_chr21_design.md      #   Dataset design rationale
│       └── hg38_chr21/               #   Human chr21 reference and data prep
│           ├── prepare_reference.sh   #     Download hg38 + extract chr21
│           └── prep_manifest.tsv      #     Manifest for subset_test_data.sh
```

## Unit Tests

| Test file | What it tests | External deps |
|-----------|---------------|---------------|
| `test_sample_sheet.py` | Sample sheet parsing, validation, analysis name construction, factor/level handling | None |
| `test_rule_commands.py` | Samtools filter/sort/dedup pipelines on synthetic SAM data | `samtools` (skips if absent) |
| `test_validate_dmc_input.py` | modBAM validation (MM/ML tags), bedMethyl format checking, auto-detection | None |
| `test_mC_helpers.py` | `parse_sample_name`, `is_dmc_sample`, `parameters_for_mc`, dmC vs Bismark routing | None |

```bash
pytest tests/unit/ -v
```

## Integration Tests

### S. pombe (primary)

18 samples, 4 assays (ChIP_broad, ChIP_narrow, RNAseq, sRNA). Reference genome is small enough to store in the repo.

| Test | What it validates | Runtime |
|------|-------------------|---------|
| `test_pombe_dryrun.py` | DAG construction, rule selection, wildcard resolution, control linking, replicate handling, checkpoints (31 tests) | ~43s |
| `test_pombe_postrun.py` | Output file existence, BAM/bigwig integrity, peak files, DEG tables, stats reports (29 tests) | ~5s |
| `test_dmc_dryrun.py` | dmC-specific rule routing: modBAM vs bedMethyl inputs, Bismark exclusion, DMR analysis | ~30s |

```bash
# Dry-run only
pytest tests/integration/test_pombe_dryrun.py -v

# Full validation cycle
scripts/validate_pombe.sh --all    # dry-run → full run → output checks
```

### Human chr21 (multi-omics)

33 samples covering all 7 assay types: ChIP_broad, ChIP_narrow, RNAseq, RAMPAGE, sRNA, WGBS, and dmC. EN-TEx stomach vs colon comparison with donors as biological replicates, plus ONT Open Data HG002 for dmC. All data subsetted to chromosome 21.

See `tests/integration/data/human_chr21_design.md` for the complete design rationale, ENCODE accessions, and known considerations.

**Data preparation** (one-time, requires SLURM and data prep tools — see `dev/docs/README.md` Section 3):

```bash
# 1. Download full hg38 reference (~3 GB) and extract chr21
bash tests/integration/data/hg38_chr21/prepare_reference.sh

# 2. Subset SRA data to chr21 (submits ~36 SLURM jobs)
bash scripts/subset_test_data.sh \
    tests/integration/data/hg38_chr21/prep_manifest.tsv chr21 \
    --outdir tests/integration/data/hg38_chr21/fastqs

# 3. Dry-run validation (works without data prep)
snakemake --configfile tests/integration/data/test_options_chr21.yaml --dry-run
```

## Running Tests

```bash
conda activate epicc

# All unit tests
pytest tests/unit/ -v

# All integration tests
pytest tests/integration/ -v

# Stop on first failure
pytest tests/unit/ -x

# Run last failed tests only
pytest tests/unit/ --lf

# Single test class or function
pytest tests/unit/test_sample_sheet.py::TestLevelsToLabel -v

# Coverage report
pytest tests/unit/ --cov=workflow/scripts --cov-report=term-missing
```

## Profiling Pipeline Runs

After a completed pipeline run:

```bash
python dev/profile_snakemake_log.py --latest                # markdown to stdout
python dev/profile_snakemake_log.py --latest --html r.html  # HTML with Gantt chart
```
