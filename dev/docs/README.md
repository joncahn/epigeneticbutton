# Developer Guide

Setup and conventions for developing and testing EPICC.

## Environment setup

One command sets up a complete dev environment — Snakemake, testing tools, aligners for test data preparation, and the GitHub CLI:

```bash
conda create -n epicc -y -c conda-forge -c bioconda \
    --file config/epicc-env.txt \
    pytest samtools gh \
    bowtie2 star bismark pigz sra-tools awscli
conda activate epicc
```

| Package | Purpose |
|---------|---------|
| `config/epicc-env.txt` | Snakemake + executor plugins (base pipeline requirement) |
| `pytest` | Test framework (`pytest tests/` commands) |
| `samtools` | Unit tests on synthetic SAM data + BAM subsetting for test data prep |
| `gh` | GitHub CLI for PRs, issues, CI interaction |
| `bowtie2` | Alignment for ChIP-seq / sRNA-seq test data subsetting |
| `star` | Alignment for RNA-seq / RAMPAGE test data subsetting |
| `bismark` | Alignment for WGBS/EMseq test data subsetting |
| `pigz` | Parallel gzip for downloaded FASTQs |
| `sra-tools` | `fasterq-dump` for SRA downloads |
| `awscli` | S3 downloads for ONT Open Data (dmC modBAM) |

The aligner and SRA packages are only needed for preparing new test datasets with `scripts/subset_test_data.sh`. For routine development (editing code and running existing tests), just Snakemake + pytest + samtools is sufficient.

### Optional tools

Not needed for routine development but useful for specific tasks. Install as needed:

```bash
conda install -n epicc -y -c conda-forge -c bioconda gffread infernal seqkit bedtools
```

| Package | Purpose |
|---------|---------|
| `gffread` | GFF/GTF conversion for preparing new test genome annotations |
| `infernal` | Structural RNA depletion databases for new reference genomes (see `Help/Structural_RNAs_Rfam.md`) |
| `seqkit` | FASTA/FASTQ inspection during test data preparation |
| `bedtools` | Ad-hoc genome arithmetic (available at runtime in conda sub-envs) |

## Running tests

All test commands assume the `epicc` env is active.

```bash
# Unit tests (fast, no data downloads)
pytest tests/unit/ -v

# Integration dry-run (validates Snakemake DAG resolution with S. pombe test data)
pytest tests/integration/test_pombe_dryrun.py -v

# mC dry-run (validates bisulfite/dmC workflow DAG across assay types)
pytest tests/integration/test_mC_dryrun.py -v

# Full validation (dry-run + pipeline execution + output checks)
scripts/validate_pombe.sh --all
```

See `tests/README.md` for detailed test documentation and `CLAUDE.md` for the full test command reference.

## Profiling pipeline runs

`dev/profile_snakemake_log.py` parses Snakemake log files and reports per-rule
timing, phase summaries, slowest individual jobs, and parallelism stats.

```bash
python dev/profile_snakemake_log.py --latest                # markdown to stdout
python dev/profile_snakemake_log.py --latest --html r.html  # self-contained HTML with Gantt chart
python dev/profile_snakemake_log.py .snakemake/log/<file>.snakemake.log
```

## Project layout (dev-relevant)

```
dev/
    docs/
        README.md             # This file
        design-decisions.md   # Architectural decisions and rationale
        sample-sheet-spec.md  # Sample sheet format specification
        benchmarks.md         # Pipeline benchmarking notes
        TODO.md               # Active roadmap and backlog
    plans/                    # Archived implementation plans (gitignored)
    profile_snakemake_log.py  # Snakemake run profiler (see above)
tests/
    unit/                     # Fast unit tests (200+)
    integration/              # Snakemake dry-run and post-run tests
        data/                 # Test data, sample sheets, and design docs
            pombe_design.md   # S. pombe test case design
            colcen_design.md  # A. thaliana ColCEN test case design
            hg38_chr21_design.md  # Human chr21 test case design
scripts/
    validate_pombe.sh         # S. pombe integration test orchestrator
    subset_test_data.sh       # SLURM-based chromosome subsetting for test data prep
    migrate_sample_sheet.py   # Old → new sample sheet format migration
tools/
    epicc-builder.html        # Sample sheet and config builder app
```

## Conventions

- **Commits**: Imperative mood, concise subject line. No co-author attribution lines.
- **Tests**: Write tests for new logic. Unit tests go in `tests/unit/test_<module>.py`. Use `@pytest.mark.slow` for anything that takes more than a few seconds.
- **Design decisions**: When making a significant architectural choice, add an entry to `dev/docs/design-decisions.md`.
- **Conda environments**: Analysis tool environments live in `workflow/envs/`. See `dev/docs/design-decisions.md` (Conda Environment Strategy) for the rationale behind the current 5-env structure.
