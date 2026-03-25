# CLAUDE.md

EPICC (Epigenetic Pipeline for Integrative Chromatin Characterization) — Snakemake pipeline for ChIP-seq, RNA-seq, small RNA-seq, bisulfite methylC-seq, and direct methylation from long-read sequencing.

## Running the Pipeline

```bash
conda create -n epicc -y --file config/epicc-env.txt && conda activate epicc
snakemake --use-conda --conda-frontend conda --cores 12      # local
snakemake --profile profiles/slurm                            # SLURM cluster
```

## Architecture

### Snakemake Structure

- `workflow/Snakefile` - Main orchestrator; loads config, parses sample metadata, includes rule modules
- `workflow/rules/` - Modular rule files by data type:
  - `environment_setup.smk` - Reference genome preparation (indexing, annotation processing)
  - `sample_download.smk` - SRA download and FASTQ processing
  - `ChIPseq.smk` - Histone/TF ChIP mapping, peak calling (MACS2), IDR analysis
  - `ATACseq.smk` - ATAC-seq mapping and peak calling (shares ChIP conda env)
  - `RNAseq.smk` - STAR alignment, differential expression (edgeR)
  - `mC.smk` - Bismark alignment, methylation calling, DMR analysis (DMRcaller)
  - `smallRNA.smk` - ShortStack analysis, structural RNA filtering
  - `combined_analysis.smk` - Cross-datatype heatmaps, metaplots, browsers (deeptools)
- `workflow/scripts/` - R and Python scripts for statistical analysis, plotting, and sample-sheet utilities
- `workflow/envs/` - Conda environment YAML files per analysis type

### Sample Sheet and Naming

Sample metadata is defined in a TSV file with 9 columns:

`Sample_ID | Assay | Genome | Levels | Replicate_ID | Read_files | Read_layout | IP_target | Control`

- **Sample_ID**: Unique identifier, used as filesystem name. Must be unique and filesystem-safe (no `__`, `/`, whitespace).
- **Assay**: Controlled vocabulary: `ChIP_broad`, `ChIP_narrow`, `ATAC`, `RNAseq`, `RAMPAGE`, `sRNA`, `WGBS`, `WGBS_nd`, `PBAT`, `EMseq`, `dmC`
- **Genome**: Reference genome name (e.g. `Spombe`, `ColCEN`)
- **Levels**: Comma-separated `factor:level` pairs (e.g. `genotype:WT,tissue:root`). All samples must have the same factors.
- **Replicate_ID**: Replicate identifier (e.g. `rep1`, `rep2`)
- **Read_files**: SRA accession (`SRR12345`), local path, HTTP(S) URL, or `+`-separated for merging multiple inputs
- **Read_layout**: `SE` or `PE`
- **IP_target**: Required for ChIP assays (e.g. `H3K9me2`, `WCE`, `Input`). Blank for others.
- **Control**: Sample_ID of the control sample (e.g. WCE or Input for ChIP). No chaining.

Per-replicate files use `Sample_ID` directly (e.g. `final__WT_H3K9me2_rep1.bam`). Analysis-level (merged replicate) files use a derived name: `{Assay}__{levels_label}__{IP_target}__{Genome}` (e.g. `ChIP_broad__WT__H3K9me2__Spombe`).

Peak type is determined by Assay: `ChIP_broad` → broad peaks (histone marks), `ChIP_narrow` → narrow peaks (transcription factors, H3K4me3, etc.). Both share the `ChIP` env (`results/ChIP/`).

Central sample-sheet logic lives in `workflow/scripts/sample_sheet.py`.

### Tools

- `tools/epicc-builder.html` - Self-contained HTML5 app for building sample sheets and options files. Tabulator-based editor with validation, dynamic factor columns, per-cell examples, and YAML options export. Open directly in a browser.
- `dev/profile_snakemake_log.py` - Snakemake log profiler. Parses `.snakemake/log/*.snakemake.log` files and reports per-rule timing, phase summary, slowest jobs, and parallelism stats. Supports markdown (stdout) and self-contained HTML with Gantt timeline chart.
  ```bash
  python dev/profile_snakemake_log.py --latest              # markdown to stdout
  python dev/profile_snakemake_log.py --latest --html r.html # HTML report
  python dev/profile_snakemake_log.py path/to/log.snakemake.log
  ```
- `scripts/subset_test_data.sh` - SLURM-based test data preparation. Downloads SRA data, aligns to a full reference, and subsets reads mapping to a target region (e.g. chr21) for use as integration test fixtures. Self-resubmitting controller with three phases (index, per-sample, gather).

### Configuration

- `config/epicc-options.yaml` - Main options file (paths, parameters, resource allocation)
  - `repo_folder` is optional; auto-detected from `workflow.basedir` at runtime (one level above the Snakefile). Override explicitly only when the repo is accessed from a non-standard path.
  - Reference genomes are namespaced under `genomes:`, each entry containing annotation file paths and species-level parameters (e.g. `genus`, `species`, `ncbi_taxid`)
  - `gtf_file`, `genomesize`, `star_index`, `ncbi_taxid`, and `structural_rna_fafile` are auto-computed at runtime (GTF derived from GFF via gffread, genome stats from FASTA, TaxId from NCBI Datasets CLI, structural RNA FASTA via Infernal/Rfam); user-provided values in the options file override the computed values
  - GO database name is auto-derived as `org.<G><species>_<GenomeName>.eg.db` (e.g. `org.Athaliana_ColCEN.eg.db`)
  - Access pattern in rule files: `config["genomes"][ref_genome][field]`
  - Old bare-key format (genome blocks as top-level keys + separate species blocks) is auto-migrated at startup with a deprecation warning
- `config/all_samples.tsv` - Sample metadata (see Sample Sheet section above)
- `profiles/slurm/config.yaml` - SLURM executor settings

## Testing

```bash
# Unit tests (sample_sheet.py logic, samtools/mC helpers, dmC input validation)
pytest tests/unit/ -v

# Integration tests — S. pombe dry-run (no cluster needed)
pytest tests/integration/test_pombe_dryrun.py -v

# Full pombe validation (dry-run + postrun checks)
scripts/validate_pombe.sh --all
```

- `tests/unit/test_sample_sheet.py` - Tests for `sample_sheet.py` functions
- `tests/unit/test_rule_commands.py` - Tests samtools pipelines on synthetic SAM data (requires samtools)
- `tests/unit/test_mC_helpers.py` - Tests mC rule helper functions (`is_dmc_sample`, `parameters_for_mc`)
- `tests/unit/test_validate_dmc_input.py` - Tests dmC input validation (modBAM MM/ML tags, bedMethyl format)
- `tests/integration/test_pombe_dryrun.py` - Snakemake dry-run validation with S. pombe test data
- `tests/integration/test_pombe_postrun.py` - Post-run output checks (requires completed pipeline run)
- `tests/integration/data/test_samples_pombe.tsv` - S. pombe test sample sheet (17 samples, 4 assays)
- `tests/integration/data/test_samples_hg38_chr21.tsv` - Human chr21 test sample sheet (33 samples, all 6 assay types)
- `tests/integration/data/test_samples_colcen.tsv` - A. thaliana ColCEN test sample sheet (38 samples: ChIP, ATAC, EMseq, PBAT, dmC)
- `tests/integration/data/test_options_colcen.yaml` - ColCEN test options (GitHub URLs for genome, SRA/URL inputs)

## Key Details

- Snakemake 9.0+. Results go to `{output_dir}/{env}/` (`ChIP`, `ATAC`, `RNA`, `sRNA`, `mC`, `combined`); genomes to `{genome_dir}/{ref_genome}/`. Both directories default to `results` and `genomes` respectively, configurable via `output_dir`/`genome_dir` config keys or `epicc --output-dir`/`--genome-dir` CLI flags.
- In Python context (input/output/params), paths use `RESULTS_DIR`/`GENOMES_DIR` variables. In shell blocks, paths use `{config[output_dir]}`/`{config[genome_dir]}` Snakemake substitution.
- Env mapping: `ChIP_broad`/`ChIP_narrow` → `ChIP`, `ATAC` → `ATAC`, `RNAseq`/`RAMPAGE` → `RNA`, `sRNA` → `sRNA`, `WGBS`/`WGBS_nd`/`PBAT`/`EMseq`/`dmC` → `mC`
- Checkpoint files in `{output_dir}/*/chkpts/` control re-running analyses; delete to force rerun.
- Read_files supports HTTP(S) URLs for FASTQ, BAM, and bedMethyl inputs. Genome config fields (`fasta_file`, `gff_file`, `te_file`) also accept URLs — downloaded automatically via curl at rule execution time.
- `te_file` accepts `.bed(.gz)` (pass-through) or `.gff3(.gz)` (auto-converted to BED6 using the GFF3 `ID=` attribute as the name column).
