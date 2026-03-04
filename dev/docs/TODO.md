# Big Refactor TODO

## Documentation

* [ ] Integrate README + Read the docs + epicc-builder app for concerted changes

* [x] Update README and CLAUDE.md to suggest creating a conda environment named epicc instead of smk9. Rename the config file epicc-env.txt. **Done**: Renamed `config/smk9.txt` → `config/epicc-env.txt`, updated all references in README.md, CLAUDE.md, validate_pombe.sh, tests/unit/README.md, test_rule_commands.py.

* [ ] we should add a list of recommendations for which of ChIP_broad/narrow to use based on the histone mark

## Analysis

### (Differential) splicing analysis for RNAseq

* [ ] heavy pipeline. Here is what was used in my MBD paper:
To look for novel splicing changes that occurred within the mC reader mutants, the reads were mapped with STAR (Dobin et al. 2013) were processed by StringTie and merged together (Pertea et al. 2015) into a master novel transcriptome comprising splicing events from TAIR10 and ones uniquely identified withinthis study. Then reads underwent lightweight alignment using Salmon version 1.4.0 (Patro et al. 2017) against the novel transcriptome from StringTie. Novel and known transcripts belonging to the same genewere analysed for splicing events by SUPPA2 (Trincado et al. 2018). Differential alternative splicing (DAS)was calculated for each event based on abundance of transcripts with and without inclusion of those eventsby SUPPA2 (Trincado et al. 2018).

### Generic pre-computed bedMethyl input support

* [ ] Support bedMethyl as a generic pre-computed methylation input format for any mC assay (WGBS, EMseq, dmC), not just dmC. Currently bedMethyl handling is exclusively gated behind dmC wildcard constraints (`_DMC_WC`). The underlying scripts (`validate_dmc_input.py`, `bedmethyl_to_cx_report.py`) are already assay-agnostic — the limitation is purely in rule routing.
  * **Current state**: A sample with `Assay: WGBS` and `Read_files: /path/to/precomputed.bed.gz` fails because (1) excluded from dmC rules by `_DMC_WC`, (2) `define_cx_report_input()` routes it to Bismark expecting FASTQs.
  * **Approach**: Add `bedMethyl` as a valid Assay in `VALID_ASSAYS` (maps to `mC` env). Modify `define_cx_report_input()` and wildcard constraints to route bedMethyl samples to the conversion pipeline. Reuse existing `get_dmc_input` → `copy_bedmethyl_input` → `convert_bedmethyl_to_cx_report` chain.
  * **Alternative**: Add an optional `Input_format` column (FASTQ/BAM/bedMethyl) to separate assay type from input format. More flexible but larger change.

### ATAC-seq input sample support

* [ ] Support calling ATAC peaks with Input

## UI/UX

### Handling IDR failures

* [ ] How can we better handle IDR failures like in the case of the S. pombe test case? It shouldn't be a failure in the pipeline sense - it's an analytical outcome, and probably here reflects the biology of S. pombe. Should we 

### Configuration fileset

* [x] Change all references to the "config file" to refer to it as the "options file". This is currently config.yaml. It should become epicc-options.yaml. This will be an extensive rename - Make sure this change is uniformly applied throughout the workflow, documentation, test cases, and the builder. Keep the config/ directory name as-is. A full run configuration is actually the composite of a sample sheet and the information contained in the options file. **Done**: Renamed `config/config.yaml` → `config/epicc-options.yaml`, all test config files → `test_options_*.yaml`. Updated all references across workflow, rules, scripts, tests, builder, and documentation.

* [ ] Parameters in the options file that are sub-settings in the yaml cannot be fed directly into snakemake command line. So for things that people might want to customize on the fly to try different things (plots mostly, but parameters for peak calling and DMRs for example), it could be good to have them as single entries.

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
  Read_files: Path to FASTA, FASTQ, BAM files, or bare SRA IDs. In this last case, read files will be downloaded from SRA. For paired-end FASTQs with separate mate files, use a comma to separate the R1 and R2 (/path/to/file.R1.fastq.gz,/path/to/file.R2.fastq.gz). If multiple read files or SRA IDs are separated by a "+", they will be merged before any processing.
  Read_layout: controlled vocabulary, single-end or paired-end sequencing (SE or PE)
  IP_target: Required for all ChIP_broad and ChIP_narrow samples, including controls. Describes what was pulled down — e.g. `H3K9me2` for an IP sample, or `Input`, `WCE`, `IgG` for a control sample.
  Control: Must be a valid Sample_ID. Used only for normalizing ChIP_broad, ChIP_narrow, and RAMPAGE samples. Validation rules: (1) the referenced Sample_ID must exist in the sheet, (2) a control sample must not itself reference another control (no chaining), (3) it is an error to specify a Control for non-ChIP/RAMPAGE assays, (4) multiple IP samples may share the same control.

##### Input validation

* [x] Per-field validation rules are defined in [`dev/docs/sample-sheet-spec.md`](sample-sheet-spec.md) (the canonical specification) and implemented in `workflow/scripts/samplefile_validation.py`. The epicc-builder app should implement the same rules. See the spec file for the full list of per-field constraints, cross-field checks, and derived name definitions.

* [x] Add a check to ensure multiple samples don't share the same inputs files or SRA accessions. This is almost certainly a user data entry error, can't think of a reasonable use case. Should be validated both in the builder and in the workflow. **Done**: Cross-row duplicate check in `samplefile_validation.py` and `epicc-builder.html`. Documented in `sample-sheet-spec.md`.

* [x] Make fatal the non-fatal warnings for PE with a single (non-BAM, non-SRA) path (i.e. just one FASTQ) and SE with more than one path passed with comma separation. **Done**: Changed from warnings to errors in both `samplefile_validation.py` and `epicc-builder.html`. Updated `sample-sheet-spec.md` severity.

* [x] epicc-builder: don't allow the user to generate an invalid sample sheet. The export button should be inactive if the sheet is in a validation failed state, and we should instruct the user to fix errors in the sheet to enable export (could be a hover tool tip, let's try to be screenreader friendly). **Done**: Export button disabled when errors exist, with tooltip and `aria-disabled` attribute.

* [x] epicc-builder: IP_target should not be editable if Assay is not chip_broad or chip_narrow. **Done**: Conditional `editable` function on IP_target column checks `CHIP_ASSAYS.has(assay)`.

* [x] epicc-builder: Unless it has been edited by the user, the Sample_ID suggestion should be updated automatically and continuously as the user changes other fields for the sample. **Done**: `_sid_user_edited` flag tracks manual edits; `applySuggestions()` auto-updates cell value when not user-edited; clearing Sample_ID resumes auto-suggestion.

#### Harmonize analysis_samplefiles

* [x] Analysis samplefiles need to be updated to follow the conventions of the new input samplefile format. **Done**: `analysis_samplefile_*.tsv` now uses new-format columns (Assay, Genome, Levels, IP_target, env, sample_name).

#### Remove old format detection in samplefile_validation.py

* [x] Old format TSVs didn't always have headers, so we can't rely on column names for detection. Let's just get rid of the detection routine and make sure the README.md informs users that they can run migrate_sample_sheet.py if they're using the old format. **Done**: Removed old-format detection from `samplefile_validation.py` and `sample_sheet.py`. Removed `OLD_COLNAMES` constant (migration script defines its own). Added usage example to README.md and sample-sheet-spec.md.

* [x] Ensure migrate_sample_sheet.py is up to date and remember to update it if any changes are made to the new sample sheet format. **Done**: Verified — handles all 9 new columns, all assay types, Control linking, and Sample_ID generation.

#### Derived Names

* [x] For samples without {IP_target}, we should avoid "____" in the filenames. Since many samples could lack an IP_target, best if we can do this by just eliminating the term and spacer, avoiding a placeholder like __NA__. This may already be the case, in which case we should just document the fact in sample-sheet-spec.md. **Done**: `build_analysis_name()` now filters empty parts before joining. E.g. `RNAseq__WT_cell__Spombe` (not `RNAseq__WT_cell____Spombe`). Documented in sample-sheet-spec.md.

#### epicc-builder

* [x] Currently, epicc-builder (as referenced in the README) is a standalone web app hosted remotely, and currently broken. We should develop a new version of epicc-builder as a self-contained HTML5/javascript app that helps users create a valid sample sheet with a tabular GUI. It will be deployed as a single HTML file that can be opened offline in any modern browser. **Done**: `tools/epicc-builder.html` (Tabulator 6.x, ~1165 lines), symlinked at repo root.

##### [x] Implementation Stage 1: Sample sheet preparation helper

**Prerequisite research**: Identify a suitable JS library for table drawing andwidgets (e.g. Handsontable, AG Grid, or a lightweight alternative that can bebundled into a single HTML file).  **Decision reached: use [tabulator] <https:/tabulator.info/>**.

**I/O**:

* [x] Output format: TSV with the column header as the first row, matching thepipeline's expected input format. **Done**.
* [x] Import an existing sample sheet (TSV) for editing. **Done**.

**UX**:

* [x] Accessibility (to e.g. screen readers) is an important concern. **Done**: aria labels, roles, keyboard-contained Tab navigation.
* [x] When focus (keyboard or mouse) is on a column, the description of that columnshould be shown to the user via a dynamic text display. **Done**: description panel updates on header/cell click.
* [x] Hamburger menu on each sample row allows for common actions like "Add areplicate", "Insert duplicate below", "Remove sample", etc. **Done**: context menu with Add replicate, Duplicate, Move up/down, Remove.
* [x] Sample rows are reorderable. **Done**: Tabulator movableRows + Move up/down in context menu.
* [x] Example text is shown for the freetext fields until a user enters information.For example, epicc-builder has the opportunity to progressively suggest uniquesample IDs as a user fills out the other fields for that sample. The user can editthe suggestion. **Done**: per-row example data (dimmed, individually cleared on edit), progressive Sample_ID suggestion with auto-accept.
* [x] Controlled-vocabulary fields are represented as a drop-down menu. **Done**: Assay, Read_layout, Control as list editors.
* [x] The drop-down for Control should be filtered based on the Assay of the currentrow to show only compatible types. For ChIP_broad and ChIP_narrow, both are validAssay sources for Control. For RAMPAGE, only RNAseq Assays are valid Controls. Ifcontrols are available, instead of an empty dropdown list, show an unselectablemessage item "None available" or similar. **Done**: CONTROL_SOURCE_ASSAYS constant drives filtering.
* [x] Replace underscores in the header row with spaces for display purposes only. **Done**.
* [x] Add a selection widget to each row to support multi-row selection and actions like delete and duplicate. **Done**: Checkbox column with `rowSelection` formatter, `selectableRows: true`, "Delete Selected" and "Duplicate Selected" toolbar buttons.

**Validation** (see [`dev/docs/sample-sheet-spec.md`](sample-sheet-spec.md) for thecanonical rules):

* [x] Perform the same input validation as the pipeline code, and give usersfeedback through diagnostic messages. **Done**: full validation engine ported from samplefile_validation.py, diagnostics panel.
* [x] Continuously evaluate user input as the table is filled out, opportunisticallyassigning defaults to sample column entries when there is sufficient input to do so. **Done**: editing-aware validation scheduling, auto IP_target clearing on Assay change.
* [x] Cells with unresolved warnings should show an indicator of this in the UI.Possibilities are yellow cell border or small yellow circle or exclamation mark inthe cell. Same for errors, but with red. **Done**: red left border for errors, yellow for warnings.
* [x] When a user enters an IP_target, we should compare with the list of knownmarks, and if the user has entered a ChIP_broad or ChIP_narrow value that conflictswith the recommendation in the known marks, we should warn. **Done**: `KNOWN_MARKS` constant with regex patterns from `config/epicc-options.yaml` peaktype mappings. Warning emitted during validation when assay conflicts with recommendation.

##### [x] Implementation Stage 2: Config file helper and concordance with original hosted app

* [x] Augment the sample sheet builder with another interactive section just above allowing for the creation, editing, removal of the factors that will be used in Levels in the sample file. This can be something like an New Factor text field and button. Committed factors will then appear as tiles showing their name with an "x" inside to enable removal. Then, in the sample sheet builder, instead of showing the Levels as a table column, we should automatically populate/update the table GUI with columns having the names of the factors, so the users will just enter the level value there. The TSV will keep the same factor:level format for the Levels column, so epicc-builder should handle this transformation on import/export. **Done**: Factor panel with add/remove tiles, dynamic `_factor_*` columns replacing Levels, import/export round-trips correctly.

* [x] Implement all currently unimplemented functionality from [the original epicc-builder](https://epicc-builder.streamlit.app/), but keep our updated sample sheet builder. **Done**: Config file form with all sections from original Streamlit app (required params, per-genome/species, output/input options, motifs, advanced ChIP/ATAC/RNA/mC/sRNA/plotting), YAML import/export via js-yaml@4.
  * [x] We will need to alter the app presentation so we can show different sections (sample sheet and epicc configuration file form settings). This could be implemented with upper nav bar with "EPICC" as the home item in the menu, followed by Sample Sheet and Config File. Clicking EPICC resolves to the first menu item (Sample Sheet). Links to documentation and the EPICC/epigeneticbutton [Github repo](https://github.com/joncahn/epigeneticbutton) should also appear in this main nav menu. **Done**: Top nav bar with EPICC brand, Sample Sheet / Config File tabs, Docs / GitHub links.
  * [x] The dedicated sample file check button probably isn't necessary since we're doing real time validation. **Done**: Real-time validation replaces the button.
  * [x] Keep all of the jokey prompts in the options file form. **Done**: All jokey prompts from the original Streamlit app preserved (full_analysis, te_analysis, GO, trimmed_fastqs, aligned_bams, motifs, nextflex, structural RNA, stranded heatmaps, browser TE, etc.).

##### [ ] Deployment

* [x] Host the new epicc-builder on Github Pages so users who are primarily working with a remote installation of epicc don't necessarily have to download the HTML file and run it locally. **Done**: GitHub Actions workflow `.github/workflows/deploy-builder.yml` deploys `tools/epicc-builder.html` as `index.html` on push to main. Requires enabling GitHub Pages (Source: GitHub Actions) in repo Settings > Pages.

##### [ ] Tweaks

* [ ] Show help text for column headers on mouse over of headers specifically, instead of just on click. Retain current show-on-click behavior for header row and value cells.

* [ ] Add an explanation of what Factors are to the right of the +Add button. Something like: "Experimental variables used for grouping samples for comparative analysis. All rows must have the same factors, and may have different levels (values)."

* [ ] Promote the Reference Genomes section out of the current Config file form to a separate top-level menu tab appearing just to the right of the EPICC brand (and becoming the new link destination for EPICC). The new layout becomes Reference Genomes | Samples | Options. The data from the References page will constrain which reference Genomes are available still be used in generating the Options file.

##### [ ] Bugs

* [ ] Sample sheet example should validate, be exportable.

* [ ] Help text for the factor levels columns should be better. Maybe should say: Levels   of factor {factor}. Stored as comma-separated factor:level pairs (e.g. genotype:WT,tissue:root) in the sample sheet TSV 'Levels' column. All rows must have the same factors. Level values form the 'levels_label' in analysis names.

### custom adapter handling

* [ ] Sequencing adapters could vary on a per-library. Maybe there should be an optional sample file column for custom adapters and we remove the global params from epicc-options.yaml. If we use skewer for trimming, auto-detection of most standard adapters is built-in if I’m not mistaken.

### species-specific parameters

* [ ] Species (as in Species-dependent parameters in epicc-options.yaml) should probably be defined as their binomial like Zea_mays to avoid collisions. Could we just get rid of that section altogether and stick ncbiID and go_database along with the params for each reference genome? We can just compute the genome size of the reference and not bother with it in the config, same with —genomeSAindexNbases, no?

### Eliminate redundant requirement for GTF

* [ ] Currently, users must supply both a GFF annotation file and GTF transcript annotation file. We provide instructions for deriving the latter from the former in the README.md, but we should instead simply try to create the GTF without requiring it from the user. If GTF creation fails we can raise a clear error message and ask the user to supply one explicitly and re-run.

### Consider using Infernal workflow for building structural_rna_depletion FASTA database

* [ ] Current suggested approach is cumbersome and results in a FASTA database that isn't ref genome specific. If we instead just run Infernal (with lots of threads and paralellized by chromosome if on the cluster) and filter overlapping hits, determine an e-value threshold, we could incorporate this into the pipeline and/or provide a subcommand to perform this task in the event that the user doesn't specify a file.

### Explicitly handle repeats vs coding gene annotations?

* [ ] Do we currently have a way to specify whether one or the other or both should be used in the analysis?

### Provide a front-end CLI executable script

* [ ] Currently, users must call snakemake to run the pipeline. Instead, we should build a project-specific executable wrapper script (called epicc) that will expose the necessary command line parameters for runtime configuration of the pipeline and calling snakemake. It will have user-friendly help text, and could provide more flexibility for reporting, logging, intermediate file cleanup, and output/working folder designation (to run in a directory other than the one in which it is invoked). We should be able to pass through arbitrary arguments to snakemake as well (maybe through a --smk argument to the wrapper).

## Plotting

* [ ] See if we can improve browser plot sample label readability

* [ ] In plotting peak stats, for now only the first 2 reps are used (empty if not, and idr only between these 2). Would be best to allow for a flexible output where all reps are shown, and all pairwise idr too. Need refactoring the way stats are compiled.

* [ ] Consider adding correlation matrix of coverage or pairwise plots between selected samples for more QC output

## Codebase Hygiene

* [ ] Rename shared routines to generic names, e.g. merging_chip_replicates → merging_bam_replicates

* [ ] Improve logging system (naming, concatenating, and cleaning if chosen)

* [ ] Input checks for different files, including extra output, e.g. browser target file with bed+label=string(not starting with -)+binsize=Integer(min1)+optional (coordinates+width)

* [ ] Merging rules to call peaks for ChIP and TSS for RAMPAGE (both macs2), for merging regions (clusters, peaks, TSS)?

## Performance/Resource Usage

### Full code review of snakemake rules

* [ ] Do a full code review of the full workflow in parallel with --simplify, and improve performance and resource usage where possible for all sample type analysis rules and the combined analysis.

### Data acquisition and preparation

* [x] Switch to direct fastq.gz downloads from ENA for download speed, better transitory disk space usage? Maybe add alternative fastq_path=ENA, or try ENA first and fall back to SRA. Look into storing SRA downloads as compressed FASTQs to avoid writing huge uncompressed data to disk, and the post-hoc wait for compression. **Done**: ENA-first downloads via `workflow/scripts/ena_download.sh` with automatic fallback to fasterq-dump. ENA provides pre-compressed `.fastq.gz`, eliminating uncompressed intermediates. `fasterq-dump --temp "$TMPDIR"` uses SLURM scratch instead of `/tmp`.

* [x] it looks like mate file compression for PE SRA accessions after fasterq-dump happens serially - the R1 file must complete before R2. I don't think there's any reason for this constraint provided sufficient resources are available. I noticed this on a local run, not sure if it's true on the cluster. **Done**: PE mate compression now runs in parallel (background `&` + `wait`), each using half the available threads.

* [x] pigz appears to be limited to 4 threads for at least local runs, which can bottleneck the pipeline when there are few samples, but with a high read volume. **Done**: New `download` resource preset with 8 threads (up from 4 in `heavy`). Download rules (`get_fastq_pe`, `get_fastq_se`) now use this preset.

### Disk usage

* [x] See if we can refactor the conda envs to decrease required disk space - current config uses ~32GiB. Are there any low-hanging dependency fruits that can be removed? Would consolidation of the majority of the analysis packages into a single environment while maintaining packages with problematic dependencies or dependency version conflicts be an optimal solution to save disk space and initial installation time by eliminating core package redundancy? **Done**: Consolidated 8 → 5 envs. Merged `epibutton_upset` into `epibutton` (removed install_github hack), `epibutton_srna` into `epibutton_rnaseq` (shared R/Bioconductor stack), `epibutton_dmc` into `epibutton_mc` (independent aligners, no conflicts). Estimated ~3-4.5 GB savings.

* [x] Should snakemake --conda-cleanup-envs be added to the pipeline to get rid of old envs? Can this be run inside of the pipeline snakemake run? **Done**: Added `conda-cleanup-pkgs: tarballs` to SLURM profile. Documented `snakemake --sdm conda --conda-cleanup-envs` as a maintenance tip in README.md (standalone command, not inside pipeline runs).

* [x] add option to keep all intermediate files, default to using pipelining and cleanup to avoid storing large intermediates like processed FASTQs, BAMs, etc. Also includes intermediate files for plotting (tracks, heatmap parameters, ...) **Done**: Added `keep_intermediates` tiered option (`none`/`standard`/`custom`/`all`) with per-category booleans (`keep_trimmed_fastqs`, `keep_final_bams`, `keep_merged_bams`, `keep_shifted_bams`, `keep_cx_reports`). `maybe_temp()` helper in Snakefile conditionally wraps outputs with `temp()`. epicc-builder exposes tier selector and custom toggles.

* [x] Review to ensure all non-retainable (via the granular retention options) intermediate files are only ever written to temp storage. **Done**: Audited all rule files; applied unconditional `temp()` to ATAC shifted BED, sRNA clean FASTQs, RNA DEG intermediates (samples/counts/RData), and GO database temp files. Added `keep_dmc_intermediates` config option for dmC aligned modBAMs and bedMethyl pileups.

* [x] Refactor merge rules to pipe `samtools merge -u` into `samtools sort`, eliminating temp intermediate BAMs. **Done**: ChIPseq `merging_chip_replicates` and RNAseq `merging_rna_replicates` now pipe directly, removing `temp_merge`/`temp` output declarations.

* [x] ChIPseq.smk - we should not be writing SAM files to disk. Wasteful of both network storage I/O (slow) and disk space. **Done**: Merged bowtie2_map + filter_chip rules into single piped rules (filter_chip_pe, filter_chip_se). Bowtie2 stdout pipes directly into samtools view → fixmate → sort, eliminating the SAM-to-disk step and reducing intermediate BAMs from 3 to 1.

### Read trimming

* [x] Explore faster options than cutadapt, like fastp or skewer, which supports automatic Illumina adapter detection (while still allowing explicit overrides in the options file). There are a number of different use cases for read trimming - let's investigate all to make sure a substitution would cover all of them. **Done**: Replaced cutadapt with fastp across the pipeline. Config restructured: `trimming_quality` (CLI strings) replaced with tool-agnostic `quality_threshold`, `min_read_length`, `trim_front` keys. Standard Illumina adapters set to `"auto"` (fastp auto-detection); non-standard adapters (NextFlex, Nextera) kept explicit. Trimming metrics changed from `.txt` to `.json`; HTML QC reports added. All 6 downstream stats rules updated to parse fastp JSON. epicc-builder config form updated. TF keys removed.

### Read Mapping

* [x] (ChIP/ATAC):  look at adding option to use [Chromap](https://github.com/haowenz/chromap) for ~10X speedup (and set as default), consider supporting different sensitivity levels if possible as with bt2. **Done**: Chromap added as default aligner for ChIP/ATAC (`chip_aligner`/`atac_aligner` config keys). Auto-falls back to bowtie2 for `repeat`/`repeatall` strategies (chromap lacks `-k 100` multi-mapping). MAPQ filter extracted from config instead of hardcoded. Dual-format metrics parsing in stats rules. epicc-builder updated with aligner selection UI.

### Local

* [x] Check snakemake log times for the S. pombe integration test to profile each stage of the DAG. **Done**: See `dev/docs/profiling-pombe.md`. Top bottleneck is `plotFingerprint` (45% of CPU time, ~14 min/job). Full pipeline completes in ~42 min on 56 cores (267/282 steps; remaining 15 are data/network failures, not code bugs).

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

* [x] Add S. pombe integration test for faster development, user installation validation, and local single-host execution as well as cluster execution. **Done**: 18 samples (11 ChIP, 4 RNA-seq, 3 sRNA), 259 pipeline steps, ~1h 11m on gemmule with 56 threads. See `tests/integration/data/test_options_pombe.yaml`.

* [x] Gather all necessary genome reference resources (fasta, gff, gtf) from [Pombase.org](https://www.pombase.org/monthly_releases/2026/pombase-2026-02-01/). Derive an appropriate test config and test samplefile. **Done**: PomBase Feb 2026 FASTA/GFF3, gffread-derived GTF, Infernal/Rfam-15.0 structural RNA FASTA (261 loci). Files in `tests/integration/data/Spombe/`.

* [x] Let’s use Hyun Soo Kim’s ChIP-seq (H3K9me2, H3K9me3, sRNA) + Ekwall lab H3K4me3 + Chang/Rct1 WCE + Martienssen 2025 RNA-seq:
  * Kim et al. 2024 GSE156069: H3K9me2/me3 WT+dcr1 PE ChIP, sRNA SE
  * Ekwall lab GSE280066: H3K4me3 SE ChIP + Input
  * Chang et al. 2017 GSE97746: WCE control PE
  * Martienssen 2025 GSE278839: RNA-seq WT+dcr1 SE

* [x] Search only R.A. Martienssen publication datasets for additional RNA-seq and ChIP libraries. **Done**: used 3 Martienssen lab datasets + 1 Ekwall lab dataset.

* [x] Any necessary data caching for the development of this test case should be done in the untracked test-data-prep directory. **Done**: structural RNA build intermediates in `test-data-prep/pombe-structural-rna/`.

* [x] Add targeted rule tests executing real samtools commands on synthetic SAM data (filter_chip_pe, filter_chip_se, filter_rna_se, bamCoverage). **Done**: `tests/unit/test_rule_commands.py` (14 tests, 12 pass + 2 skip when bamCoverage unavailable).

* [x] Add post-run validation tests checking all pipeline outputs for existence, format integrity, and correct structure. **Done**: `tests/integration/test_pombe_postrun.py` (29 tests, marked `@slow`, auto-skips if no completed run).

* [x] Add validate-pombe orchestration script and Claude Code skill. **Done**: `scripts/validate_pombe.sh` (`--dry`/`--full`/`--check`/`--all`), `.claude/commands/validate-pombe.md` (`/validate-pombe` skill).

* [x] Remove tissue:cell from test_samples_pombe sample sheet - no reason to have this factor as it is not differentiated among the samples. **Done**: Reduced Levels to single `genotype:X` factor, updated all Sample_IDs, analysis names, and test assertions across 5 files.

* [x] Fix pombe test sample metadata: rename Ekwall `veg` → `WT` (HU3112 is WT), fix Kim/Chang WCE controls from `IP_target: Input` → `WCE` (SRR5445712 is whole-cell extract, not Input). **Done**: Updated test_samples_pombe.tsv, test_pombe_dryrun.py, test_pombe_postrun.py, test_sample_sheet.py.

### Complete A. thaliana ColCEN Chr5 test case

* [ ] Add a more complete A. thaliana ColCEN Chr5 test case, using tests/integration/test_samples_chr5.tsv, tests/integration/test_samples_colcen.tsv, and the test data we have already prepared at test-data-prep/ as sources. The idea is to create a Chr5 test subset of all of the samples currently used in test_samples_colcen.tsv. Make sure we subset any input fastqs and BAMs to contain only reads mapped to Chr5. This may require alignment to the full ColCEN genome first, and then a samtools view to subset.

### H. sapiens test case

* [ ] Add H. sapiens Chr21 test case.

## Known Unknowns

* [ ] For now, ChIPseq replicates are only properly merged if same paired information (all PE or all SE). Not sure what happens if both PE and SE reps are available with the same line+tissue name. Corner case to check.

## Known Issues/Bugs

* [ ] PlotPCA can fail if no dimensions found. check npz results before starting PCA? (how in bash?)
