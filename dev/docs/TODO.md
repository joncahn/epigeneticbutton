# Big Refactor TODO

## Documentation

* Integrate README + Read the docs + epicc-builder app for concerted changes

## UI

* config.yaml:  one thing we can think of as well, maybe for the future big reworking, is that parameters in the config file that are sub-settings in the yaml cannot be fed directly into snakemake command line. So for things that people might want to customize on the fly to try different things (plots mostly, but parameters for peak calling and DMRs for example), it could be good to have them as single entries.

| Assay   | Genome   | Condition1 | Condition2 | Label | ReplicateID | Read_files | Read_layout | Control | IP_target |
|---------|----------|------------|------------|-------|-------------|------------|-------------|---------|-----------|
| [ChIP_broad, ChIP_narrow, ATAC, RNAseq, RAMPAGE, sRNA, WGBS, dmC, EMseq] | [freetext] | [freetext] | [freetext] | [freetext] default = [Assay]\_[Genome]\_[Condition1]\_[Condition2] | [freetext] | FASTQ SE:[/path/to/file/name.r1.fq], FASTQ PE:[/path/to/file/name.r1.fq,/path/to/file/name.r2.fq], BAM SE or PE: [/path/to/file/name.bam], SRA: [SRRxxxxx], SRA merge multiple: [SRRxxxxx+SRRxxxxx+SRRxxxxx] | [SE or PE] | [label of Input] required for ChIP/TF | [Input or name of target, e.g. H3K9me2] required for ChIP/TF |

* refactor sample sheet, especially fastq_path - this actually already not strictly a fastq_path - we search for bams and bedmethyls here. Need to think more about whether a one file per line approach to the sample file makes more sense than the current approach. See: sample-sheet-refactor-recommendations.md
  * explicit placeholder value for empty data to improve user experience with sample sheet, e.g. “-” synonym for “” in TSV
  * data_type and sample_type columns are the only targets of this refactor, all others remain the same for now
* Sequencing adapters could vary on a per-library. Maybe there should be an optional sample file column for custom adapters and we remove the global params from config.yaml. If we use skewer for trimming, auto-detection of most standard adapters is built-in if I’m not mistaken.
* Species (as in Species-dependent parameters in config.yaml) should probably be defined as their binomial like Zea_mays to avoid collisions. Could we just get rid of that section altogether and stick ncbiID and go_database along with the params for each reference genome? We can just compute the genome size of the reference and not bother with it in the config, same with —genomeSAindexNbases, no?

## Plotting

* See if we can improve browser plot sample label readability
* In plotting peak stats, for now only the first 2 reps are used (empty if not, and idr only between these 2). Would be best to allow for a flexible output where all reps are shown, and all pairwise idr too. Need refactoring the way stats are compiled.
* Consider adding correlation matrix of coverage or pairwise plots between selected samples for more QC output

## Codebase Hygiene

* Rename shared routines to generic names, e.g. merging_chip_replicates → merging_bam_replicates
* Should we use the Snakemake Wrapper Repository? Looks very actively maintained:
    <https://github.com/snakemake/snakemake-wrappers>
* Improve logging system (naming, concatenating, and cleaning if chosen)
* Input checks for different files, including extra output (e.g. browser target file with bed+label=string(not starting with -)+binsize=Integer(min1)+optional (coordinates+width)
* Merging rules to call peaks for ChIP  and TSS for RAMPAGE (both macs2), for merging regions (clusters, peaks, TSS)?

## Performance/Resource Usage

* Read trimming: consider switching to something faster than cutadapt, like fastp or skewer, which supports automatic Illumina adapter detection (while still allowing explicit overrides in the config.yaml)
* Intermediate data retention: add option to keep all intermediate files, default to using pipelining and cleanup to avoid storing large intermediates like processed FASTQs, BAMs, etc. Also includes intermediate files for plotting (tracks, heatmap parameters, ...)
* Mapping (ChIP/ATAC):  look at adding option to use Chromap for ~10X speedup (and possibly set as default), consider supporting different sensitivity levels if possible as with bt2.
* Reconsider the resource request delineations, and maybe add more fine-grained/task-specific options like proc-intensive/himem-lowproc/mapping, etc.
* See if we can refactor the conda envs to decrease required disk space - current config uses ~32GiB
* Examine wall clock times on Elzar and estimate reasonable requests with slop - most steps should be O(n) for sequence inputs, e.g. trimming and mapping, but downstream analysis might differ.
* ChIPseq.smk - we should not be writing SAM files to disk. Wasteful of both network storage I/O (slow) and disk space.
* Switch to direct fastq.gz downloads from ENA for download speed, better transitory disk space usage? Maybe add alternative fastq_path=ENA
* Make better use of temporary storage on the cluster nodes to reduce NFS I/O bottlenecks and minimize temporary disk usage bloat
* Small one: mem_mb and tmp_mb should be changed to mem_mib and tmp_mib to meet expectations of binary byte counting.
* Add time for all rules + in the profiles config
* Resolve slurm issues with QOSMaxSubmitJobPerUserLimit reached sometimes (when it should be limited to 16 in the profile (specific to CSHL cluster, but could be helpful for other environments in case it' a shared bug)

## Testing

* Add pombe test case for faster development (and user installation validation) and local machine run instead of cluster - let’s use some Hyun Soo ChIP-seq (H3K9me2, H3K4me3, PolII, Inputs), RNA-seq, sRNA-seq, pombe ATAC?, long read gDNA (methylation neg ctrl)?
    <https://www.nature.com/articles/s41467-024-53417-9#data-availability> <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156069>
* Add test coverage for pre-existing and new features.
* For now, ChIPseq replicates are only properly merged if same paired information (all PE or all SE). Not sure what happens if both PE and SE reps are available with the same line+tissue name. Corner case to check.

## BUGS

* PlotPCA can fail if no dimensions found. check npz results before starting PCA? (how in bash?)