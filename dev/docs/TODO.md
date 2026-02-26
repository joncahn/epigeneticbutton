# Big Refactor TODO

## Documentation

* [ ] Integrate README + Read the docs + epicc-builder app for concerted changes

## UI/UX

### refactor sample sheet

* [ ] New sample sheet format streamlines and clarifies input specifications. This is a major change that will require a deep review of the codebase, especially the snakemake rules. Automatically generated filenames will change, e.g. "Input" -> "Control", as controls could be arbitrary samples (WCE for yeast ChIP, RNA-seq for RAMPAGE, etc.). There will also be significant changes to documentation to reflect this new format, and to example files and test code. Finally, this change should be validated with a full run of the S. pombe integration test.

  | Assay   | Genome   | Condition1 | Condition2 | Replicate_ID | Read_files | Read_layout | Group_ID |   Control | IP_target |
  |---------|----------|------------|------------|-------|-------------|------------|-------------|  ---------|-----------|
  | [ChIP_broad, ChIP_narrow, ATAC, RNAseq, RAMPAGE, sRNA, WGBS, dmC, EMseq] | [freetext] | [freetext] |   [freetext] | [freetext] | FASTQ SE:[/path/to/file/name.r1.fq], FASTQ PE:[/path/to/file/name.r1.fq,/path/to/file/name.r2.fq], BAM SE or PE: [/path/to/file/name.bam], SRA: [SRRxxxxx], SRA merge multiple:   [SRRxxxxx+SRRxxxxx+SRRxxxxx] | [SE or PE] | [freetext] | [boolean] | [Input or freetext name of IP   target, e.g. H3K9me2] required for ChIP_broad/ChIP_narrow |
  
  Assay: controlled vocabulary, replaces data_type/sample_type and provides the menu of accepted assay   types for analysis.
  Genome: Reference genome name
  Condition1: An experimental variable like sample line/genotype (e.g. B73), tissue type, environmental   parameter value (e.g. 37deg), time point (e.g. T0), etc. Condition1 and Condition2 generalize the previous line/tissue fields.
  Condition2: (optional) another experimental variable. If provided, will be analyzed combinatorially   with Condition1.
  Replicate_ID: name your replicates (e.g. rep1, rep2, repA, repB, 1, 2). If multiple samples share the   same Replicate_ID they will be merged as technical replicates in the analysis.
  Read_files: Path to FASTQ, BAM files, or bare SRA IDs. In this last case, read files will be downloaded from SRA. For paired-end FASTQs with separate mate files, use a comma to separate the R1 and R2 (/path/to/file.R1.fastq.gz,/path/to/file.R2.fastq.gz).
  Read_layout: controlled vocabulary, single-end or paired-end sequencing (SE or PE)
  Group_ID: groups samples together for normalization with a control sample
  Control: [boolean] - indicate whether this sample is the control sample for all samples with this   Group_ID. We should validate that there is just one control per Group_ID. Controls are currently used   only for ChIP_broad, ChIP_narrow, and RAMPAGE samples.
  IP_target: Required only for ChIP_broad and ChIP_narrow samples.
  
  * All free text fields should have input validation for safety, and path validation and SRA regex validation should be applied to Read_files.

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

* [ ] (ChIP/ATAC):  look at adding option to use Chromap for ~10X speedup (and possibly set as default), consider supporting different sensitivity levels if possible as with bt2.

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

### Schizosaccharomyces pombe test case

* [x] Add S. pombe integration test for faster development, user installation validation, and local single-host execution as well as cluster execution. **Done**: 18 samples (11 ChIP, 4 RNA-seq, 3 sRNA), 259 pipeline steps, ~1h 11m on gemmule with 56 threads. See `tests/integration/data/test_config_pombe.yaml`.

* [x] Gather all necessary genome reference resources (fasta, gff, gtf) from [Pombase.org](https://www.pombase.org/monthly_releases/2026/pombase-2026-02-01/). Derive an appropriate test config and test samplefile. **Done**: PomBase Feb 2026 FASTA/GFF3, gffread-derived GTF, Infernal/Rfam-15.0 structural RNA FASTA (261 loci). Files in `tests/integration/data/Spombe/`.

* [x] Let’s use Hyun Soo Kim’s ChIP-seq (H3K9me2, H3K9me3, sRNA) + Ekwall lab H3K4me3 + Chang/Rct1 WCE + Martienssen 2025 RNA-seq:
    - Kim et al. 2024 GSE156069: H3K9me2/me3 WT+dcr1 PE ChIP, sRNA SE
    - Ekwall lab GSE280066: H3K4me3 SE ChIP + Input
    - Chang et al. 2017 GSE97746: WCE control PE
    - Martienssen 2025 GSE278839: RNA-seq WT+dcr1 SE

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