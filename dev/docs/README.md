# Developer Guide

Setup and conventions for developing and testing EPICC.

## Environment setup

### 1. Base environment (required)

The base `epicc` env provides Snakemake and its executor plugins. This is sufficient to run the pipeline (analysis tools are installed automatically by Snakemake into per-rule conda environments).

```bash
conda create -n epicc -y --file config/epicc-env.txt
conda activate epicc
```

### 2. Developer extras

For development and testing, install these additional packages into the same environment:

```bash
conda install -n epicc -y -c conda-forge -c bioconda \
    pytest \
    samtools \
    gh
```

| Package | Purpose |
|---------|---------|
| `pytest` | Test framework. Required for all `pytest tests/` commands. |
| `samtools` | Required by `tests/unit/test_rule_commands.py` (executes real samtools pipelines on synthetic SAM data). Tests are skipped if samtools is not found. |
| `gh` | GitHub CLI. Useful for creating PRs, viewing issues, and CI interaction from the command line. Not required for tests. |

Or as a single create command that sets up a ready-to-go dev environment:

```bash
conda create -n epicc -y -c conda-forge -c bioconda \
    --file config/epicc-env.txt \
    pytest samtools gh
```

### 3. Test data preparation (for subsetting new test datasets)

The `scripts/subset_test_data.sh` script downloads SRA data, aligns to a full reference, subsets reads to a target region, and regenerates FASTQs. It requires aligners and SRA tools in the `epicc` env:

```bash
conda install -n epicc -y -c bioconda -c conda-forge \
    bowtie2 star bismark samtools pigz sra-tools awscli
```

| Package | Purpose |
|---------|---------|
| `bowtie2` | Alignment for ChIP-seq and sRNA-seq samples |
| `star` | Alignment for RNA-seq and RAMPAGE samples |
| `bismark` | Alignment for WGBS/EMseq samples |
| `samtools` | BAM subsetting, sorting, indexing, FASTQ extraction |
| `pigz` | Parallel gzip compression of downloaded FASTQs |
| `sra-tools` | `fasterq-dump` for downloading from SRA |
| `awscli` | S3 downloads for ONT Open Data (dmC modBAM) |

See `tests/integration/data/hg38_chr21/prep_manifest.tsv` for the data prep manifest.

### 4. Optional tools

These are not needed for routine development but are useful for specific tasks:

| Package | Purpose |
|---------|---------|
| `gffread` | GFF/GTF format conversion. Needed when preparing new test genome annotations (e.g. deriving GTF from GFF3). |
| `infernal` | RNA homology search. Used to build structural RNA depletion databases for new reference genomes (see `Help/Structural_RNAs_Rfam.md`). |
| `seqkit` | Sequence toolkit. Handy for inspecting/subsetting FASTA/FASTQ during test data preparation. |
| `bedtools` | Genome arithmetic. Available in the base `epicc` env's conda environments at runtime, but useful to have in the dev env for ad-hoc inspection. |

Install any of these as needed:

```bash
conda install -n epicc -y -c conda-forge -c bioconda gffread infernal seqkit bedtools
```

## Running tests

All test commands assume the `epicc` env is active.

```bash
# Unit tests (fast, no data downloads)
pytest tests/unit/ -v

# Integration dry-run (validates Snakemake DAG resolution with S. pombe test data)
pytest tests/integration/test_pombe_dryrun.py -v

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
        TODO.md               # Active roadmap and backlog
    plans/                    # Archived implementation plans (gitignored)
    profile_snakemake_log.py  # Snakemake run profiler (see above)
tests/
    unit/                     # Fast unit tests
    integration/              # Snakemake dry-run and post-run tests
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
