# Big Refactor TODO

## Features

### (Differential) splicing analysis for RNAseq

* [ ] Jon: heavy pipeline. Here is what was used in my MBD paper:
To look for novel splicing changes that occurred within the mC reader mutants, the reads mapped withSTAR (Dobin et al. 2013) were processed by StringTie and merged together (Pertea et al. 2015) into amaster novel transcriptome comprising splicing events from TAIR10 and ones uniquely identified withinthis study. Then reads underwent lightweight alignment using Salmon version 1.4.0 (Patro et al. 2017) against the novel transcriptome from StringTie. Novel and known transcripts belonging to the same genewere analysed for splicing events by SUPPA2 (Trincado et al. 2018). Differential alternative splicing (DAS)was calculated for each event based on abundance of transcripts with and without inclusion of those eventsby SUPPA2 (Trincado et al. 2018).

### Generic pre-computed bedMethyl input support

* [ ] Support bedMethyl as a generic pre-computed methylation input format for any mC assay (WGBS, EMseq, dmC), not just dmC. Currently bedMethyl handling is exclusively gated behind dmC wildcard constraints (`_DMC_WC`). The underlying scripts (`validate_dmc_input.py`, `bedmethyl_to_cx_report.py`) are already assay-agnostic — the limitation is purely in rule routing.
  * **Current state**: A sample with `Assay: WGBS` and `Read_files: /path/to/precomputed.bed.gz` fails because (1) excluded from dmC rules by `_DMC_WC`, (2) `define_cx_report_input()` routes it to Bismark expecting FASTQs.
  * **Approach**: Add `bedMethyl` as a valid Assay in `VALID_ASSAYS` (maps to `mC` env). Modify `define_cx_report_input()` and wildcard constraints to route bedMethyl samples to the conversion pipeline. Reuse existing `get_dmc_input` → `copy_bedmethyl_input` → `convert_bedmethyl_to_cx_report` chain.
  * **Alternative**: Add an optional `Input_format` column (FASTQ/BAM/bedMethyl) to separate assay type from input format. More flexible but larger change.

### ATAC-seq input sample support

* [ ] Support calling ATAC peaks with Input

## Documentation

* [ ] Integrate README + Read the docs + epicc-builder app for concerted changes

## UI/UX

### refactor sample sheet

#### New format

* [x] New sample sheet format streamlines and clarifies input specifications. Controls can now be arbitrary samples (WCE for yeast ChIP, RNA-seq for RAMPAGE, etc.) referenced by Sample_ID rather than a boolean flag. **Done**: New format with 9 columns (Sample_ID, Assay, Genome, Levels, Replicate_ID, Read_files, Read_layout, IP_target, Control). Central logic in `workflow/scripts/sample_sheet.py`. Migration script at `scripts/migrate_sample_sheet.py`.
  * [x] Update Snakefile sample-sheet parsing to read the new columns and build sample metadata accordingly. **Done**: `read_sample_sheet()` + `add_compat_columns()` in Snakefile.
  * [x] Update rule files: automatically generated filenames will change (e.g. "Input" -> the control's Sample_ID). **Done**: All 8 rule files migrated (ChIPseq, RNAseq, ATACseq, smallRNA, mC, combined_analysis, sample_download, environment_setup). TF env eliminated.
  * [x] Update documentation (README, Read the Docs, example sample sheets) to reflect the new format. **Done**: README.md and CLAUDE.md updated.
  * [x] Update test sample sheets and test code. **Done**: All 4 test sample sheets converted (pombe, colcen, chr5, dmc). Unit tests for sample_sheet.py (49 tests).
  * [x] Validate with a full run of the S. pombe integration test. **Done**: 257 pipeline steps, all completed successfully.

  | Sample_ID | Assay | Genome | Levels | Replicate_ID | Read_files | Read_layout | IP_target | Control |
  |-----------|-------|--------|--------|--------------|------------|-------------|-----------|---------|
  | [freetext] | [ChIP_broad, ChIP_narrow, ATAC, RNAseq, RAMPAGE, sRNA, WGBS, dmC, EMseq] | [freetext] | [freetext] | [freetext] | FASTQ SE:[/path/to/file/name.r1.fq], FASTQ PE:[/path/to/file/name.r1.fq,/path/to/file/name.r2.fq], BAM SE or PE: [/path/to/file/name.bam], SRA: [SRRxxxxx], SRA merge multiple:   [SRRxxxxx+SRRxxxxx+SRRxxxxx] | [SE or PE] | [freetext name of IP target or control, e.g. H3K9me2, or Input, WCE, etc.] required for ChIP_broad/ChIP_narrow | [valid sample ID] |
  
  Sample_ID: a name that uniquely identifies this sample. Will be used to track the sample internally, and can be used to assign controls to ChIP_broad, ChIP_narrow, and RAMPAGE Assays. We will not enforce any format, other than uniqueness, but the epicc-builder app should suggest a concise ID (see epicc-builder specification).
  Assay: controlled vocabulary, replaces data_type/sample_type and provides the menu of accepted assay   types for analysis.
  Genome: Reference genome name
  Levels: Comma-separated list of factor:level pairs describing the experimental conditions for this sample. Uses the syntax `factor_name:level_value` (e.g. `genotype:WT,tissue:root`). Factors can be any experimental variable — genotype (B73, Mo17), tissue (root, leaf), temperature (37deg, 24deg), time point (T0, T1), etc. If multiple factors are specified, multifactorial comparisons will be performed. Factor names are not currently used in pipeline logic (only levels are used for constructing comparisons and plot labels), but they improve readability and enable epicc-builder to parse and validate entries progressively. **All samples must have the same number of factor:level pairs, including controls.** Controls should still specify meaningful levels where possible (e.g. `genotype:WT,tissue:root`).
  Replicate_ID: identifies biological or technical replicates (e.g. rep1, rep2, repA, repB, 1, 2). Replicates are treated independently in analysis and merged only for specific downstream steps like peak calling. Note: `+` in Read_files is for concatenating multi-file inputs from the same library (e.g. multiple lanes or runs) before any processing — this is distinct from replicate handling.
  Read_files: Path to FASTQ, BAM files, or bare SRA IDs. In this last case, read files will be downloaded from SRA. For paired-end FASTQs with separate mate files, use a comma to separate the R1 and R2 (/path/to/file.R1.fastq.gz,/path/to/file.R2.fastq.gz). If multiple read files or SRA IDs are separated by a "+", they will be merged before any processing.
  Read_layout: controlled vocabulary, single-end or paired-end sequencing (SE or PE)
  IP_target: Required for all ChIP_broad and ChIP_narrow samples, including controls. Describes what was pulled down — e.g. `H3K9me2` for an IP sample, or `Input`, `WCE`, `IgG` for a control sample.
  Control: Must be a valid Sample_ID. Used only for normalizing ChIP_broad, ChIP_narrow, and RAMPAGE samples. Validation rules: (1) the referenced Sample_ID must exist in the sheet, (2) a control sample must not itself reference another control (no chaining), (3) it is an error to specify a Control for non-ChIP/RAMPAGE assays, (4) multiple IP samples may share the same control.

##### Input validation

* [x] Per-field validation rules are defined in [`dev/docs/sample-sheet-spec.md`](sample-sheet-spec.md) (the canonical specification) and implemented in `workflow/scripts/samplefile_validation.py`. The epicc-builder app should implement the same rules. See the spec file for the full list of per-field constraints, cross-field checks, and derived name definitions.

#### epicc-builder

* [ ] Currently, epicc-builder (as referenced in the README) is a standalone web app hosted remotely, and currently broken. We should develop a new version of epicc-builder as a self-contained HTML5/javascript app that helps users create a valid sample sheet with a tabular GUI. It will be deployed as a single HTML file that can be opened offline in any modern browser.

  **Prerequisite research**: Identify a suitable JS library for table drawing and widgets (e.g. Handsontable, AG Grid, or a lightweight alternative that can be bundled into a single HTML file).

  **I/O**:
  * [ ] Output format: TSV with the column header as the first row, matching the pipeline's expected input format.
  * [ ] Import an existing sample sheet (TSV) for editing.

  **UX**:
  * [ ] Accessibility (to e.g. screen readers) is an important concern.
  * [ ] When focus (keyboard or mouse) is on a column, the description of that column should be shown to the user via a dynamic text display.
  * [ ] Hamburger menu on each sample row allows for common actions like "Add a replicate", "Insert duplicate below", "Remove sample", etc.
  * [ ] Sample rows are reorderable.
  * [ ] Controlled-vocabulary fields are represented as a drop-down menu.
  * [ ] Example text is shown for the freetext fields until a user enters information. For example, epicc-builder has the opportunity to progressively suggest unique sample IDs as a user fills out the other fields for that sample. The user can edit the suggestion.

  **Validation** (see [`dev/docs/sample-sheet-spec.md`](sample-sheet-spec.md) for the canonical rules):
  * [ ] Perform the same input validation as the pipeline code, and give users feedback through diagnostic messages.
  * [ ] Continuously evaluate user input as the table is filled out, opportunistically assigning defaults to sample column entries when there is sufficient input to do so.

### config.yaml

* [ ] one thing we can think of as well, maybe for the future big reworking, is that parameters in the config file that are sub-settings in the yaml cannot be fed directly into snakemake command line. So for things that people might want to customize on the fly to try different things (plots mostly, but parameters for peak calling and DMRs for example), it could be good to have them as single entries.

### custom adapter handling

* [ ] Sequencing adapters could vary on a per-library. Maybe there should be an optional sample file column for custom adapters and we remove the global params from config.yaml. If we use skewer for trimming, auto-detection of most standard adapters is built-in if I’m not mistaken.

### species-specific parameters

* [ ] Species (as in Species-dependent parameters in config.yaml) should probably be defined as their binomial like Zea_mays to avoid collisions. Could we just get rid of that section altogether and stick ncbiID and go_database along with the params for each reference genome? We can just compute the genome size of the reference and not bother with it in the config, same with —genomeSAindexNbases, no?

### Eliminate redundant requirement for GTF

* [ ] Currently, users must supply both a GFF annotation file and GTF transcript annotation file. We provide instructions for deriving the latter from the former in the README.md, but we should instead simply try to create the GTF without requiring it from the user. If GTF creation fails we can raise a clear error message and ask the user to supply one explicitly and re-run.

### Consider using Infernal workflow for building structural_rna_depletion FASTA database

* [ ] Current suggested approach is cumbersome and results in a FASTA database that isn't ref genome specific. If we instead just run Infernal (with lots of threads and paralellized by chromosome if on the cluster) and filter overlapping hits, determine an e-value threshold, we could incorporate this into the pipeline and/or provide a subcommand to perform this task in the event that the user doesn't specify a file.

### Explicitly handle repeats vs coding gene annotations?

Do we currently have a way to specify whether one or the other or both should be used in the analysis?

## Plotting

* [ ] See if we can improve browser plot sample label readability

* [ ] In plotting peak stats, for now only the first 2 reps are used (empty if not, and idr only between these 2). Would be best to allow for a flexible output where all reps are shown, and all pairwise idr too. Need refactoring the way stats are compiled.

* [ ] Consider adding correlation matrix of coverage or pairwise plots between selected samples for more QC output

## Codebase Hygiene

* [ ] Rename shared routines to generic names, e.g. merging_chip_replicates → merging_bam_replicates

* [ ] Should we use the Snakemake Wrapper Repository? Looks very actively maintained:
    <https://github.com/snakemake/snakemake-wrappers>

* [ ] Improve logging system (naming, concatenating, and cleaning if chosen)

* [ ] Input checks for different files, including extra output, e.g. browser target file with bed+label=string(not starting with -)+binsize=Integer(min1)+optional (coordinates+width)

* [ ] Merging rules to call peaks for ChIP and TSS for RAMPAGE (both macs2), for merging regions (clusters, peaks, TSS)?

## Performance/Resource Usage

### Data acquisition and preparation

* [ ] Switch to direct fastq.gz downloads from ENA for download speed, better transitory disk space usage? Maybe add alternative fastq_path=ENA, or try ENA first and fall back to SRA. Look into storing SRA downloads as compressed FASTQs to avoid writing huge uncompressed data to disk, and the post-hoc wait for compression.

* [ ] pigz appears to be limited to 4 threads for at least local runs, which can bottleneck the pipeline when there are few samples, but with a high read volume.

* [ ] it looks like mate file compression for PE SRA accessions after fasterq-dump happens serially - the R1 file must complete before R2. I don't think there's any reason for this constraint provided sufficient resources are available. I noticed this on a local run, not sure if it's true on the cluster.

### Disk usage

* [ ] See if we can refactor the conda envs to decrease required disk space - current config uses ~32GiB. Are there any low-hanging dependency fruits that can be removed? Would consolidation of at least some of the environments save disk space by eliminating core package redundancy?

* [ ] Should snakemake --conda-cleanup-envs be added to the pipeline to get rid of old envs? Can this be run inside of the pipeline snakemake run?

* [ ] add option to keep all intermediate files, default to using pipelining and cleanup to avoid storing large intermediates like processed FASTQs, BAMs, etc. Also includes intermediate files for plotting (tracks, heatmap parameters, ...)

* [ ] ChIPseq.smk - we should not be writing SAM files to disk. Wasteful of both network storage I/O (slow) and disk space.

### Read trimming

* [ ] consider switching to something faster than cutadapt, like fastp or skewer, which supports automatic Illumina adapter detection (while still allowing explicit overrides in the config.yaml)

### Read Mapping

* [ ] (ChIP/ATAC):  look at adding option to use [Chromap](https://github.com/haowenz/chromap) for ~10X speedup (and possibly set as default), consider supporting different sensitivity levels if possible as with bt2.

* [ ] (WGBS):  consider switching to [bwa-meth](https://github.com/brentp/bwa-meth)

### Local

* [ ] Check snakemake log times for the S. pombe integration test to profile each stage of the DAG.

### Cluster

* [ ] Reconsider the resource request delineations, and maybe add more fine-grained/task-specific options like proc-intensive/himem-lowproc/mapping, etc.

* [ ] Examine wall clock times on Elzar and estimate reasonable requests with slop - most steps should be O(n+k) for sequence inputs, e.g. trimming and mapping, but downstream analysis might differ.

* [ ] Make better use of temporary storage on the cluster nodes to reduce NFS I/O bottlenecks and minimize temporary disk usage bloat

* [ ] Small one: mem_mb and tmp_mb should be changed to mem_mib and tmp_mib to meet expectations of binary byte counting.

* [ ] Add time for all rules + in the profiles config

* [ ] Resolve slurm issues with QOSMaxSubmitJobPerUserLimit reached sometimes (when it should be limited to 16 in the profile (specific to CSHL cluster, but could be helpful for other environments in case it' a shared bug)

## Testing

### Pico and EMseq

* [ ] Check that Pico and EMseq work (following adapter trimming improvement mentioned above)

### Schizosaccharomyces pombe test case

* [x] Add S. pombe integration test for faster development, user installation validation, and local single-host execution as well as cluster execution. **Done**: 18 samples (11 ChIP, 4 RNA-seq, 3 sRNA), 259 pipeline steps, ~1h 11m on gemmule with 56 threads. See `tests/integration/data/test_config_pombe.yaml`.

* [x] Gather all necessary genome reference resources (fasta, gff, gtf) from [Pombase.org](https://www.pombase.org/monthly_releases/2026/pombase-2026-02-01/). Derive an appropriate test config and test samplefile. **Done**: PomBase Feb 2026 FASTA/GFF3, gffread-derived GTF, Infernal/Rfam-15.0 structural RNA FASTA (261 loci). Files in `tests/integration/data/Spombe/`.

* [x] Let’s use Hyun Soo Kim’s ChIP-seq (H3K9me2, H3K9me3, sRNA) + Ekwall lab H3K4me3 + Chang/Rct1 WCE + Martienssen 2025 RNA-seq:
  * Kim et al. 2024 GSE156069: H3K9me2/me3 WT+dcr1 PE ChIP, sRNA SE
  * Ekwall lab GSE280066: H3K4me3 SE ChIP + Input
  * Chang et al. 2017 GSE97746: WCE control PE
  * Martienssen 2025 GSE278839: RNA-seq WT+dcr1 SE

* [x] Search only R.A. Martienssen publication datasets for additional RNA-seq and ChIP libraries. **Done**: used 3 Martienssen lab datasets + 1 Ekwall lab dataset.

* [x] Any necessary data caching for the development of this test case should be done in the untracked test-data-prep directory. **Done**: structural RNA build intermediates in `test-data-prep/pombe-structural-rna/`.

### Complete A. thaliana ColCEN Chr5 test case

* [ ] Add a more complete A. thaliana ColCEN Chr5 test case, using tests/integration/test_samples_chr5.tsv, tests/integration/test_samples_colcen.tsv, and the test data we have already prepared at test-data-prep/ as sources. The idea is to create a Chr5 test subset of all of the samples currently used in test_samples_colcen.tsv. Make sure we subset any input fastqs and BAMs to contain only reads mapped to Chr5. This may require alignment to the full ColCEN genome first, and then a samtools view to subset.

### H. sapiens test case

* [ ] Add H. sapiens Chr21 test case.

## Known Issues

* [ ] For now, ChIPseq replicates are only properly merged if same paired information (all PE or all SE). Not sure what happens if both PE and SE reps are available with the same line+tissue name. Corner case to check.

## BUGS

* [ ] PlotPCA can fail if no dimensions found. check npz results before starting PCA? (how in bash?)
