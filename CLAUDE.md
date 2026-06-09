# CLAUDE.md

EPICC (Epigenetic Pipeline for Integrative Chromatin Characterization) — Snakemake pipeline for ChIP-seq, RNA-seq, small RNA-seq, bisulfite methylC-seq, and direct methylation from long-read sequencing.

## Running the Pipeline

```bash
conda create -n epicc -y --file config/epicc-env.txt && conda activate epicc
epicc run --samples your_samples.tsv --cores 12              # local
epicc run --samples your_samples.tsv --profile profiles/slurm # SLURM cluster
epicc validate --samples your_samples.tsv                     # config check + dry-run
epicc validate --build-envs --samples your_samples.tsv        # pre-create rule conda envs
                                                              # (run from a login/dev node
                                                              # before sbatch-wrapped runs
                                                              # on clusters where conda
                                                              # env create fails inside a
                                                              # SLURM allocation)
```

`epicc` is a thin argparse wrapper around snakemake; subcommands are `run`, `validate`, `output`, `unlock`, `perf`, `clean`, `init-profile`. Anything after `--` is forwarded verbatim to snakemake. Raw `snakemake ...` invocations still work but bypass `epicc`'s placeholder-detection and TMPDIR routing.

`epicc init-profile TYPE [NAME]` copies a bundled cluster profile template (`slurm` or `uge`) to `~/.config/snakemake/NAME/config.yaml` — the XDG path Snakemake searches when `--profile NAME` is used. Recommended setup for conda-install users and the standard path when sharing a profile across multiple projects. Use `--force` to overwrite an existing profile.

`run` and `validate` accept named overrides for commonly-tuned settings (all override the options YAML):
- `--chip-aligner {bowtie2,chromap}`, `--atac-aligner`, `--chip/atac-mapping-strategy`
- `--[no-]full-analysis`, `--[no-]te-analysis`, `--[no-]go`, `--[no-]rna-deduplicate`, `--methylation-contexts CG,CHG,CHH`
- `--dmr-caller {metilene,dmrcaller}`, `--dmr-min-diff SPEC` (`auto` / float / `CG=f,CHG=f`), `--dmr-min-cytosines`, `--dmr-p-value`, `--dmr-sigma-n`
- `--chip-callpeaks-params STR`, `--cut-broad-caller {epic2,seacr,macs2}`, `--cut-narrow-caller {seacr,macs2}`

Nested settings (`dmr_thresholds.*`, `chip_callpeaks.params`, `cut_callpeaks.*`) are applied via a temp merged config file; sibling keys not specified on the CLI are preserved from the options YAML.

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
- **Assay**: Controlled vocabulary: `ChIP_broad`, `ChIP_narrow`, `CUT_RUN_broad`, `CUT_RUN_narrow`, `CUT_TAG_broad`, `CUT_TAG_narrow`, `ATAC`, `RNAseq`, `RAMPAGE`, `sRNA`, `WGBS`, `WGBS_nd`, `PBAT`, `EMseq`, `dmC`
- **Genome**: Reference genome name (e.g. `Spombe`, `ColCEN`)
- **Levels**: Comma-separated `factor:level` pairs (e.g. `genotype:WT,tissue:root`). All samples must have the same factors.
- **Replicate_ID**: Replicate identifier (e.g. `rep1`, `rep2`)
- **Read_files**: SRA accession (`SRR12345`), local path, HTTP(S) URL, or `+`-separated for merging multiple inputs (`+`-merge supported for SRA accessions and FASTQ files, local or URL; BAM/bedMethyl must be merged upstream)
- **Read_layout**: `SE` or `PE`
- **IP_target**: Required for ChIP assays (e.g. `H3K9me2`, `WCE`, `Input`). Blank for others.
- **Control**: Sample_ID of the control sample (e.g. WCE or Input for ChIP). No chaining.

Control-row replicate merging keys on `(Levels, IP_target, Genome)` only — the `Assay` value on a control row is decorative for merging purposes, so a single biological Input/IgG/WCE serving multiple IP types (broad + narrow, ChIP + CUT&RUN) merges correctly regardless of how individual rep rows are labeled. See `build_control_merge_key` in `workflow/scripts/sample_sheet.py`.

Per-replicate files use `Sample_ID` directly (e.g. `final__WT_H3K9me2_rep1.bam`). Analysis-level (merged replicate) files use a derived name: `{Assay}__{levels_label}__{IP_target}__{Genome}` (e.g. `ChIP_broad__WT__H3K9me2__Spombe`).

Peak type is determined by Assay: `ChIP_broad`/`CUT_RUN_broad`/`CUT_TAG_broad` → broad peaks (histone marks), `ChIP_narrow`/`CUT_RUN_narrow`/`CUT_TAG_narrow` → narrow peaks (transcription factors, H3K4me3, etc.). All six "IP-with-peaks" assays share the `ChIP` env (`results/ChIP/`). Default peak callers: ChIP* → MACS2; CUT&* `_broad` → epic2; CUT&* `_narrow` → SEACR. Override via `cut_callpeaks.{broad,narrow}_caller` (`epic2`, `seacr`, `macs2`).

Central sample-sheet logic lives in `workflow/scripts/sample_sheet.py`.

### Tools

- `tools/epicc-builder.html` - Self-contained HTML5 app for building sample sheets and options files. Tabulator-based editor with validation, dynamic factor columns, per-cell examples, and YAML options export. Open directly in a browser.
- `dev/profile_snakemake_log.py` - Snakemake log profiler. Parses `.snakemake/log/*.snakemake.log` files and reports per-rule timing, phase summary, slowest jobs, and parallelism stats. Supports markdown (stdout) and self-contained HTML with Gantt timeline chart. Default mode auto-aggregates all logs from the same resumed run, identified by `output_dir` + `analysis_name` wildcards extracted from log content (not timestamps).
  ```bash
  python dev/profile_snakemake_log.py                       # aggregate latest run
  python dev/profile_snakemake_log.py --html r.html         # HTML report
  python dev/profile_snakemake_log.py --single              # newest log only
  python dev/profile_snakemake_log.py path/to/log.snakemake.log
  ```
- `scripts/subset_test_data.sh` - SLURM-based test data preparation. Downloads SRA data, aligns to a full reference, and subsets reads mapping to a target region (e.g. chr21) for use as integration test fixtures. Self-resubmitting controller with three phases (index, per-sample, gather).

### Configuration

- `config/epicc-options.yaml` - Main options file (paths, parameters, resource allocation)
  - `repo_folder` is optional; auto-detected from `workflow.basedir` at runtime (one level above the Snakefile). Override explicitly only when the repo is accessed from a non-standard path.
  - Reference genomes are namespaced under `genomes:`, each entry containing annotation file paths and species-level parameters (e.g. `genus`, `species`, `ncbi_taxid`)
  - `gtf_file`, `genomesize`, `star_index`, `ncbi_taxid`, and `structural_rna_fafile` are auto-computed at runtime (GTF derived from GFF via gffread, genome stats from FASTA, TaxId from NCBI Datasets CLI, structural RNA FASTA via Infernal/Rfam); user-provided values in the options file override the computed values
  - GO database name is auto-derived as `org.<G><species>.eg.db` (e.g. `org.Athaliana.eg.db`) — matches AnnotationForge's strict `org.<G><species>.eg.db` package-name format. To keep multiple reference genomes of the same species (e.g. ColCEN + TAIR10) from colliding in the conda env's shared R library, each genome's GO package is installed into and loaded from `genomes/<refgenome>/GO/` via per-call `lib=` / prepended `.libPaths()`.
  - Access pattern in rule files: `config["genomes"][ref_genome][field]`
  - Old bare-key format (genome blocks as top-level keys + separate species blocks) is auto-migrated at startup with a deprecation warning
  - `methylation_contexts` (default `["CG", "CHG", "CHH"]`) gates per-context mC analysis: bigwigs, DMR calls, and PCA plots are produced only for listed contexts. Set to `["CG"]` for animal genomes where non-CpG methylation is negligible. Subcontexts (CAG/CAA/...) not currently supported.
  - `dmr_thresholds:` block configures DMR calling parameters for both callers. Shared: `min_diff` (per-context, default CG=0.3/CHG=0.2/CHH=0.1), `min_cytosines` (default 5). DMRcaller-only: `bin_size`, `p_value`, `min_gap`, `min_size`, `min_reads`. Metilene-only: `max_cpgs`, `valley` (per-context), `maxseg` (per-context). Defaults are Arabidopsis-leaf-tuned; adjust for other organisms/tissues.
  - `chip_aligner` / `atac_aligner` (default `"bowtie2"`) — set to `"chromap"` for ~10x faster mapping. chromap is appropriate for SE-only datasets; PE data shows only 11–28% properly-paired reads with chromap vs 99%+ with bowtie2 (chromap maps ends independently, so one mate is MAPQ-filtered when repeat paralogs are present). When `chip_mapping_strategy` / `atac_mapping_strategy` is `repeat` or `repeatall`, bowtie2 is used automatically regardless of the aligner setting.
  - `use_node_tmpdir` (default `false`) toggles TMPDIR routing — see Key Details below.
- `config/example_samples.tsv` - Documented sample-sheet template (copy and edit; pass to epicc via `--samples`)
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
- Env mapping: `ChIP_broad`/`ChIP_narrow`/`CUT_RUN_broad`/`CUT_RUN_narrow`/`CUT_TAG_broad`/`CUT_TAG_narrow` → `ChIP`, `ATAC` → `ATAC`, `RNAseq`/`RAMPAGE` → `RNA`, `sRNA` → `sRNA`, `WGBS`/`WGBS_nd`/`PBAT`/`EMseq`/`dmC` → `mC`
- Checkpoint files in `{output_dir}/*/chkpts/` control re-running analyses; delete to force rerun.
- Read_files supports HTTP(S) URLs for FASTQ, BAM, and bedMethyl inputs. Genome config fields (`fasta_file`, `gff_file`, `te_file`) also accept URLs — downloaded automatically via curl at rule execution time.
- `te_file` accepts `.bed(.gz)` (pass-through) or `.gff3(.gz)` (auto-converted to BED6 using the GFF3 `ID=` attribute as the name column).
- TMPDIR routing: `workflow/Snakefile` registers a `shell.prefix` that points TMPDIR at `{output_dir}/.tmp/${SLURM_JOB_ID:-$$}/` per job, with a `trap rm` cleanup on EXIT. Tools that spill through TMPDIR (sort, samtools, STAR, fasterq-dump, deeptools) write to the project filesystem instead of cluster `/tmp` — important on sites where `/tmp` is a RAM-tmpfs sized off `mem_mb` (e.g. SLURM `JobContainerType=job_container/tmpfs`). Opt out via `use_node_tmpdir: true` in the options file or `epicc run --use-node-tmpdir`. Bismark is the one mapper with its own `--temp_dir` pointed under `results/mC/mapped/<sample>/`, so it's unaffected by the TMPDIR override either way.
- dmC (ONT/PacBio) flow uses `modkit pileup` → bedMethyl → CX_report rather than DMRcaller's `readONTbam`. Reasons: (1) modkit is the de facto standard ONT methylation caller with calibrated probability filtering and broader mode support (5mC/5hmC/dual); (2) bedMethyl is a portable intermediate (IGV, methylpy, bsmooth can all consume it); (3) routing dmC through the same CX_report format as Bismark means the entire downstream pipeline (bigwigs, DMRs, PCA) has one code path; (4) `readONTbam`'s defaults are mammal-centric (`BSgenome.Hsapiens.UCSC.hg38`, CG-only, single 0.5 probability threshold) and would need plant BSgenome packages per genome and three separate calls per replicate for our 3-context use case.
- DMR pair-context jobs (`call_DMRs_for_pair_context`) run `computeDMRsReplicates` with `parallel=FALSE` and `cores=1L` — serial per-chromosome iteration. Parallel mode (`SnowParam(type="SOCK")` across chromosomes) consistently crashed BiocParallel's reducer on heterogeneous plant data with `Error: BiocParallel errors / wrong args for environment subassignment`, masking the underlying worker R error. A diagnostic with identical input + `parallel=FALSE` completes cleanly (verified PBAT__Col0_seedling vs WGBS__Col0_seedling CHG, 2026-05-18, 2147 DMRs in 56 min). Trade-off: high-coverage 4-replicate pairs (e.g. EMseq×EMseq on ColCEN-scale plant data) take 1–4 h serial vs ~10–20 min parallel-when-it-worked; bins method is still fast enough for this to be acceptable. The pooled fallback (`computeDMRs` for N<2 pairs) keeps parallel SnowParam since it doesn't trigger the bug.
- The shipped `profiles/slurm/config.yaml` sets `retries: 1` to absorb both transient SLURM control-plane glitches and the rare SnowParam port collisions that can still hit the pooled-fallback DMR path. Site profiles inheriting from it should keep `retries` at 1 or higher.
