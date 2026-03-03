# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EpigeneticButton (EPICC - Epigenetic Pipeline for Integrative Chromatin Characterization) is a Snakemake-based bioinformatics pipeline for analyzing and integrating epigenomics datasets: ChIP-seq, RNA-seq, small RNA-seq, bisulfite methylC-seq, and direct methylation from long-read sequencing.

## Running the Pipeline

```bash
# Install environment
conda create -n epicc -y --file config/epicc-env.txt
conda activate epicc

# Run locally
snakemake --use-conda --conda-frontend conda --cores 12

# Run on SLURM cluster
snakemake --profile profiles/slurm

# Run in background with logging
snakemake --profile profiles/slurm > epigeneticbutton.log 2>&1 &

# Generate DAG to validate configuration
snakemake --dag | dot -Tpng > dag.png
```

### Intermediate Targets

```bash
# Mapping only (no analysis)
snakemake --profile profiles/slurm map_only

# ChIP coverage bigwigs only
snakemake --profile profiles/slurm coverage_chip
```

## Architecture

### Snakemake Structure

- `workflow/Snakefile` - Main orchestrator; loads config, parses sample metadata, includes rule modules
- `workflow/rules/` - Modular rule files by data type:
  - `environment_setup.smk` - Reference genome preparation (indexing, annotation processing)
  - `sample_download.smk` - SRA download and FASTQ processing
  - `ChIPseq.smk` - Histone/TF ChIP mapping, peak calling (MACS2), IDR analysis
  - `RNAseq.smk` - STAR alignment, differential expression (edgeR)
  - `mC.smk` - Bismark alignment, methylation calling, DMR analysis (DMRcaller)
  - `smallRNA.smk` - ShortStack analysis, structural RNA filtering
  - `combined_analysis.smk` - Cross-datatype heatmaps, metaplots, browsers (deeptools)
- `workflow/scripts/` - R scripts for statistical analysis and plotting
- `workflow/envs/` - Conda environment YAML files per analysis type

### Sample Sheet and Naming

Sample metadata is defined in a TSV file with 9 columns:

`Sample_ID | Assay | Genome | Levels | Replicate_ID | Read_files | Read_layout | IP_target | Control`

- **Sample_ID**: Unique identifier, used as filesystem name. Must be unique and filesystem-safe (no `__`, `/`, whitespace).
- **Assay**: Controlled vocabulary: `ChIP_broad`, `ChIP_narrow`, `ATAC`, `RNAseq`, `RAMPAGE`, `sRNA`, `WGBS`, `EMseq`, `dmC`
- **Genome**: Reference genome name (e.g. `Spombe`, `ColCEN`)
- **Levels**: Comma-separated `factor:level` pairs (e.g. `genotype:WT,tissue:root`). All samples must have the same factors.
- **Replicate_ID**: Replicate identifier (e.g. `rep1`, `rep2`)
- **Read_files**: SRA accession (`SRR12345`), local path, or `+`-separated for merging multiple inputs
- **Read_layout**: `SE` or `PE`
- **IP_target**: Required for ChIP assays (e.g. `H3K9me2`, `WCE`, `Input`). Blank for others.
- **Control**: Sample_ID of the control sample (e.g. WCE or Input for ChIP). No chaining.

Per-replicate files use `Sample_ID` directly (e.g. `final__WT_H3K9me2_rep1.bam`). Analysis-level (merged replicate) files use a derived name: `{Assay}__{levels_label}__{IP_target}__{Genome}` (e.g. `ChIP_broad__WT__H3K9me2__Spombe`).

Peak type is determined by Assay: `ChIP_broad` → broad peaks (histone marks), `ChIP_narrow` → narrow peaks (transcription factors, H3K4me3, etc.). Both share the `ChIP` env (`results/ChIP/`).

Central sample-sheet logic lives in `workflow/scripts/sample_sheet.py`.

### Tools

- `tools/epicc-builder.html` - Self-contained HTML5 app for building sample sheets and config files. Tabulator-based editor with validation, dynamic factor columns, per-cell examples, and YAML config export. Open directly in a browser.

### Configuration

- `config/config.yaml` - Main configuration (paths, parameters, resource allocation)
- `config/all_samples.tsv` - Sample metadata (see above)
- `profiles/slurm/config.yaml` - SLURM executor settings

### Output Structure

Results go to `results/{env}/` where env is one of: `ChIP`, `ATAC`, `RNA`, `sRNA`, `mC`, `combined`. Each contains `chkpts/` (checkpoint files for pipeline logic), `logs/`, `tracks/` (bigwigs), and analysis-specific subdirectories.

Reference genomes are prepared in `genomes/{ref_genome}/`.

## Testing

```bash
# Unit tests (sample_sheet.py logic, samtools rule commands)
pytest tests/unit/ -v

# Integration tests — S. pombe dry-run (no cluster needed)
pytest tests/integration/test_pombe_dryrun.py -v

# Full pombe validation (dry-run + postrun checks)
scripts/validate_pombe.sh --all
```

- `tests/unit/test_sample_sheet.py` - Tests for `sample_sheet.py` functions
- `tests/unit/test_rule_commands.py` - Tests samtools pipelines on synthetic SAM data (requires samtools)
- `tests/integration/test_pombe_dryrun.py` - Snakemake dry-run validation with S. pombe test data
- `tests/integration/test_pombe_postrun.py` - Post-run output checks (requires completed pipeline run)
- `tests/integration/data/test_samples_pombe.tsv` - S. pombe test sample sheet (18 samples, 4 assays)

## Key Implementation Details

- Requires Snakemake 9.0+
- Control samples are linked explicitly via the `Control` column (a valid Sample_ID). Controls can be any sample (WCE, H3, IgG, etc.).
- Peak types are determined by Assay (`ChIP_broad` → broad, `ChIP_narrow` → narrow). No regex config needed.
- Env mapping: `ChIP_broad`/`ChIP_narrow` → `ChIP`, `ATAC` → `ATAC`, `RNAseq`/`RAMPAGE` → `RNA`, `sRNA` → `sRNA`, `WGBS`/`EMseq`/`dmC` → `mC`
- RNA-seq strandedness is configurable per protocol (RNAseq vs RAMPAGE)
- Checkpoint files in `results/*/chkpts/` control re-running analyses; delete to force rerun
