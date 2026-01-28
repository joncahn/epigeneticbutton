# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EpigeneticButton (EPICC - Epigenetic Pipeline for Integrative Chromatin Characterization) is a Snakemake-based bioinformatics pipeline for analyzing and integrating epigenomics datasets: ChIP-seq, RNA-seq, small RNA-seq, bisulfite methylC-seq, and direct methylation from long-read sequencing.

## Running the Pipeline

```bash
# Install environment
conda create -n smk9 -y --file config/smk9.txt
conda activate smk9

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

### Sample Naming Convention

Samples are identified by a compound name: `{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}` (double underscore separators).

Data types: `ChIP`, `ChIP_<group>`, `TF_<name>`, `RNAseq`, `sRNA`, `mC`

Sample types for mC data_type:
- Bisulfite sequencing: `mC`, `WGBS`, `Pico`, `EMseq` (processed via Bismark)
- Direct methylation: `dmC` (native base modifications from Oxford Nanopore/PacBio; processed via modkit)

### Configuration

- `config/config.yaml` - Main configuration (paths, parameters, resource allocation)
- `config/all_samples.tsv` - Sample metadata (9 columns: data_type, line, tissue, sample_type, replicate, seq_id, fastq_path, paired, ref_genome)
- `profiles/slurm/config.yaml` - SLURM executor settings

### Output Structure

Results go to `results/{env}/` where env is one of: `ChIP`, `TF`, `RNA`, `sRNA`, `mC`, `combined`. Each contains `chkpts/` (checkpoint files for pipeline logic), `logs/`, `tracks/` (bigwigs), and analysis-specific subdirectories.

Reference genomes are prepared in `genomes/{ref_genome}/`.

## Key Implementation Details

- Requires Snakemake 9.0+
- ChIP Input samples must have `sample_type: Input` regardless of actual control type (H3, IgG, etc.)
- TF ChIP uses `TF_<name>` data_type to link IP with corresponding Input
- Peak types (narrow/broad) are determined by regex patterns in `chip_callpeaks.peaktype` config
- RNA-seq strandedness is configurable per protocol (RNAseq vs RAMPAGE)
- Checkpoint files in `results/*/chkpts/` control re-running analyses; delete to force rerun
