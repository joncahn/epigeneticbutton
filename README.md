# Epigenetic Pipeline for Integrative Chromatin Characterization (EPICC) 
# aka Epigenetic Button

A Snakemake-based pipeline for analyzing and integrating various types of (epi)genomics datasets, including histone and transcription factor ChIP-seq, RNA-seq, RAMPAGE, small RNA-seq, and methylC-seq.

## Complete documentation: [Read the docs](https://epicc-documentation.readthedocs.io/en/latest/)

## Overview

EpigeneticButton is a comprehensive pipeline that processes and analyzes multiple types of genomics data. It provides an automated workflow for:
- Data preprocessing and quality control
- Read mapping and alignment
- Peak calling and differential expression analysis
- Data integration and visualization

## Features

- **Multiple Data Types Support**:
  - Histone ChIP-seq
  - Transcription Factor ChIP-seq
  - RNA-seq
  - small RNA-seq
  - MethylC-seq (mC) - Bisulfite sequencing via Bismark
  - Direct Methylation (dmC) - Long-read native methylation (ONT, PacBio)
  - RAMPAGE *\*in development*

- **Automated Analysis**:
  - Reference genome preparation
  - Sample-specific processing
  - Data type-specific analysis
  - Combined analysis across samples
  - Quality control and reporting
  - Additional output options such as heatmaps, metaplots and browsers

- **Flexible Configuration**:
  - [App](https://epicc-builder.streamlit.app/) to validate configuration options
  - Customizable mapping parameters
  - Configurable analysis options
  - Resource management
  - Parallel processing

## Installation

1. Clone the repository:
```bash
git clone https://github.com/joncahn/epigeneticbutton.git
```
or for ssh connection
```bash
git clone git@github.com:joncahn/epigeneticbutton.git
```
```bash
cd epigeneticbutton
```

2. Install snakemake and other required packages from the depency file:
```bash
conda create -n smk9 -y --file config/smk9.txt
conda activate smk9
```
If you want to run the pipeline on a different platform than locally or slurm, you will need to also install the corresponding snakemake-executor-plugin

## Usage

### Configuration

For new users, it is recommended to use the configuration app to validate your sample metadata file and choose analysis options:\
https://epicc-builder.streamlit.app/

1. Prepare your sample metadata file (default to `config/all_samples.tsv`) with the required columns below (see Input requirements for more details specific to each data-type):
   - `data_type`: Type of data [RNAseq | ChIP_* | TF_* | mC | sRNA] (RAMPAGE under development)
   - `line`: Sample line (e.g. B73)
   - `tissue`: Tissue type
   - `sample_type`: Sample identifier
   - `replicate`: Replicate ID
   - `seq_id`: Sequencing ID; use the corresponding SRR####### if downloading from SRA
   - `fastq_path`: Path to FASTQ files; if downloading from SRA, use "SRA" 
   - `paired`: [PE | SE]
   - `ref_genome`: Reference genome name

2. Update `config/config.yaml` with your paths and parameters:
   - Sample file: this is the full path to the file detailed above which contain your samples metadata. 
   - Reference genome files: for each reference genome in the sample file (last column), enter the full path of a fasta file, a gene gff file, and a gene gtf file (See [below](#common-to-all-types-of-samples) for more details)
   - Analysis parameters / options
   - Species-specific parameters
   - Resources allocation
   
3. If you are running the pipeline on a different platform than CSHL slurm cluster, you will likely need to adjust the rule-specific resource parameters at the bottom of the `config/config.yaml` and the config file for your cluster scheduler (`profiles/slurm/config.yaml` for SLURM or create a new profile for your scheduler). In slurm, the default is to start 16 jobs maximum in parallel. Keep in mind that units in the cluster file are in MB.

4. Default options: 
   - Full analysis: By default, a full analysis is performed form raw data to analysis plots. Change `full_analysis` in the config file ([see below](#main-output-options)).
   - Limited QC output: By default, some QC options are not performed to limit the time and amount of output files. Change `QC_option` in the config file ([see below](#main-output-options)).
   - No Gene Ontology analysis: Due to the difficulty in automating building a GO database, this option is OFF by default. Change `GO` option in the config file. Please refer to Additional output options #2 below and [Help GO](Help/Gene_Ontology.md) before setting it to `true` as it requires 2 other files. These files are available for Arabidopsis thaliana (Tair10 / ColCEN assembly) and Maize B73 (v5 or NAM assembly) in the `data` folder.
   - No TE analysis: By default, no analysis on transposable elements is performed. Change `te_analysis` in the config file ([see below](#main-output-options)).
   - For ChIP-seq: the default mapping parameters are bowtie2 `--end-to-end` default parameters. Other options are available in the config file `chip_mapping_option` ([see below](#chip-mapping-parameters)).
   - For sRNA-seq: the default is not based on Netflex v3 library preparation. If your data was made with this kit, an additional deduplication and read trimming is required. To turn it ON, change the `Netflex_v3_deduplication` in the config file. See [Known issues #3](#known-potential-issues) below if you have mixed libraries.
   - For sRNA-seq: the default is not to filter structural RNAs prior to shortstack analysis. Change `structural_rna_depletion` in the config file.  While this step is recommended for small interfering RNA analysis, it requires a pre-build database of fasta files. Please refer to the [Help structural RNAs](Help/Structural_RNAs_Rfam.md) before setting it to `true`. This file is available for Maize in the `data` folder.
   - For sRNA-seq: the default is to only perform *de novo* micro RNA identification (`--dn_mirna` argument in ShortStack). If you also want the known microRNAs, download the fasta file from [miRbase](https://www.mirbase.org), filter it for your species of interest, and add to the `srna_mapping_params` entry in the config file `--known_miRNAs <path/to/known_miRNA_file.fa>`.

### Running the Pipeline

1. To run the pipeline locally:
```bash
snakemake --use-conda --conda-frontend conda --cores 12
```

2. To run the pipeline on a HPC-slurm (using sbatch):
```bash
snakemake --profile profiles/slurm
```

If you do not want all the snakemake output (very talkative), instead of using `--quiet` I would recommmend redirecting it to a log and putting the run in the background:
```bash
snakemake --profile profiles/slurm > epigeneticbutton.log 2>&1 &
```

If you do not want all the snakemake output (very talkative), instead of using `--quiet` I would recommmend redirecting it to a log and putting the run in the background:
```bash
snakemake --profile profiles/slurm > epigeneticbutton.log 2>&1 &
```

3. Other option: To run the pipeline on a UGE cluster (using qsub):
```bash
mkdir hpclogs
snakemake --profile profiles/uge
```

*The commands for the clusters are specific to the CSHL environement. If using a profile, make sure these parameters are adapted to your cluster too or edit accordingly.*

4. Optional: to test the pipeline, consider generating a DAG first to make sure your samplefiles and parameters work:
```bash
snakemake --dag | dot -Tpng > dag.png
```

*Even if snakemake is launched on a cluster with a profile option, the run will output a lot on the terminal. It is recommended to launch the command from a screen, to start it from a script submitted to the cluster, or to put the command in the background (which will still output snakemake commands but allows further action).*\
*For full understanding of snakemake capabilities and options: https://snakemake.readthedocs.io/en/stable/*

## Sample file configuration

A template and more details can be found on the epicc-builder app:
https://epicc-builder.streamlit.app/
You can also use it to validate your entries.

### Common to all types of samples:
- Col2: *line*: Can be any information you want, such as `Col0` or `WT` to annotate and label samples
- Col3: *tissue*: Can be any information you want, such as `leaf` or `mutant` or `6h_stress` to annotate and label samples
The combination line x tissue will be the base for all comparisons (e.g `WT_leaf` vs `WT_roots` or `Col0_control` vs `Ler_stress`)
- Col5: *replicate*: Any value to match the different replicates (e.g Rep1, RepA, 1). All the different replicates are merged for samples with the same line, tissue and sample_type.
- Col6: *seq_id*: Unique identifier to identify the raw data. If the data is deposited in SRA, it can be an SRR number (e.g. SRR27821931) or a comma-delimited list of SRR numbers (e.g. SRR27821931,SRR27821932,SRR27821933) without spaces if multiple fastq files should be merged into 1 biological replicate. If the data is local, it must be a unique identifier of the file in this folder (e.g. `wt_k27`). This identifier should be shared by both 'R1' and 'R2' fastq files for paired-end data.
- Col7: *fastq_path*: Either `SRA` if raw data to be downloaded from SRA (the SRR number should be used as `seq-id`), or the path to the directory containing the fastq file (e.g. `/archive/fastq`), in which case the `seq_id` should be a unique identifier of the corresponding file in this folder (e.g. `/archive/fastq/raw.reads.wt_k27.fastq.gz`)
- Col8: *paired*: `PE` for paired-end data or `SE` for single-end data. PE samples should have two fastq files R1 and R2 at the location defined above, sharing the same identifier in Col6 (e.g. `/archive/fastq/raw.reads.wt_k27_R1.fastq.gz` and `/archive/fastq/raw.reads.wt_k27_R2.fastq.gz`)
- Col9: *ref_genome*: Name of the reference genome to use for mapping (e.g `tair10`). 
For each reference genome, a corresponding fasta, gff and gtf files are required. It can be a full path (including the extension) or relative to the main repo folder. These files can be gzipped. For example, if your sample file has `B73_NAM` as a reference genome (last column), there must be this entry in the config file:
```
B73_NAM:
	fasta_file: path/to/B73.fasta	# can be .fa(.gz) or .fasta(.gz)
	gff_file: B73.gff	# can be .gff*(.gz)
	gtf_file: B73.gtf	# can be .gtf(.gz)
```
Other files specific to each reference genome are optional.
The GTF file can be created from a GFF file with cufflinks `gffread -T <gff_file> -o <gtf_file>` and check that `transcript_id` and `gene_id` are correctly assigned in the 9th column. The GFF file should have `gene` and `exon` in the 3rd column. All files can be gzipped (.gz extension).

### Histone ChIP-seq
- Col1: *data_type*: `ChIP` or `ChIP_<id>` where `<id>` is an identifier to relate an IP sample to its corresponding input. Only necessary in case there are different inputs to be used for different IP samples that otherwise share the same `line` and `tissue` values.
For example: If you have H3K27meac IP samples which you want compared to an H3 sample, and H4K16ac to be compared to H4 samples. Both H3 and H4 samples should be labeled `Input` in sample_type, so to differentiate them, use `ChIP_H3` and `ChIP_H4` for their data_type and the ones of H3K27ac and H4K16ac, respectively. Example:\
`ChIP_H3	Col0	WT	H3K27ac	Rep1	wt_k27	./fastq/	PE	Tair10`\
`ChIP_H3	Col0	WT	IP	Rep1	wt_h3_ctrl	./fastq/	PE	Tair10`\
`ChIP_H4	Col0	WT	H3K27ac	Rep1	wt_h4k16	./fastq/	PE	Tair10`\
`ChIP_H4	Col0	WT	IP	Rep1	wt_h4_ctrl	./fastq/	PE	Tair10`
- Col4: *sample_type*: Either `Input` to be used as a control (even if it is actually H3 or IgG pull-down), or the histone mark IP (e.g. H3K9me2). If the mark is not already listed in the config file `chip_callpeaks: peaktype:`, add it to the desired category (either narrow or broad peaks).
- Option: Differential nucleosome sensitivity (DNS-seq) can be analyzed with `ChIP` data_type, using `MNase` for the light digest and `Input` for the heavy digest.

### Transcription factor ChIP-seq
- Col1: *data_type*: `TF_<tf_name>` where `<tf_name>` is the name of the transcirption factor (e.g. for `TB1` data, use `TF_TB1`). This name should be identical for the IP and its input, and for all replicates. Multiple TFs can be analyzed in parallel, each having its own set of IP and Input samples e.g. `TF_<name1>` and `TF_<name2>`.
- Col4: *sample_type*: Either `Input` or `IP`. This works for transcription factors with narrow peaks (default). Use `IPb` for broad peaks.

### RNA-seq
- Col1: *data_type*: `RNAseq`. No other options.
- Col4: *sample_type*: `RNAseq`. No other options.

### small RNA-seq
- Col1: *data_type*: `sRNA`. No other options.
- Col4: *sample_type*: `sRNA`. Can also be `smallRNA` or `shRNA`. Does not really matter but it is used in file names.

### Whole Genome Bisulfite Sequencing
- Col1: *data_type*: `mC`. No other options.
- Col4: *sample_type*: `mC`, `WGBS`, `Pico`, or `EMseq`. These labels help identify the chemical or enzymatic conversion method, but all are processed through the Bismark pipeline.

### Direct Methylation (Long-Read Sequencing)
- Col1: *data_type*: `mC`. No other options.
- Col4: *sample_type*: `dmC`. This identifies samples with native base modifications (Oxford Nanopore, PacBio) that have not undergone bisulfite conversion or other enzymatic treatments. However, this sample type can also be used with any upstream methylation analysis that produces modBAM or bedMethyl files.
- Col7: *fastq_path*: Path to input file or directory containing input files. Supports:
  - **modBAM**: BAM files with MM/ML methylation tags from basecalling (e.g., Dorado, Guppy)
  - **bedMethyl**: Pre-computed methylation calls in bedMethyl format (e.g., from modkit pileup)
- Col6: *seq_id*: Unique identifier used to locate files in the fastq_path directory. When fastq_path is a directory, files matching `*seq_id*.bam` (for modBAM) or `*seq_id*.bed*` (for bedMethyl) will be automatically detected. If both formats exist, bedMethyl is preferred as it's pre-computed. When fastq_path is a direct file path, seq_id serves as a sample identifier.
- Col8: *paired*: Only `SE` is currently supported for dmC samples (we assume long-read sequencing).
- Note: The pipeline automatically detects whether input is modBAM or bedMethyl format. modBAM files are aligned if necessary and processed through modkit pileup. Both formats are converted to a unified Bismark-compatible CX_report format for downstream analysis (bigwig generation, DMR calling) compatible with bisulfite samples.

## Configuration Options

More details can be found on the epicc-builder app or commented within the `config/config.yaml` file:
https://epicc-builder.streamlit.app/

### Main output options
- `full_analysis`: When `false`, only the mapping and the bigwigs will occur. When `true` (default), will also be performed: single-data analyses (e.g. peak calling for ChIP, differential expression for RNAseq, DMRs for mC) and combined analyses (e.g. Upset plots for ChIP/TF, heatmaps and metaplots on all genes).
- `te_analysis`: When `true`, small RNA differential expression will be performed (if such data is available), as well as heatmaps and metaplots of all the samples on TEs. The name and path to the TE file in bed format must be filled in the config file for the corresponding reference genome. The name of each TE (4th column of the bed file) must be unique. Default is `false`.
- `QC_option`: When `true`, runs fastQC on raw and trimmed fastq files. Default is `false`.

### Intermediate input formats
- `trimmed_fastqs`: When `false` (default), the analysis runs from raw, untrimmed fastq files and performs adapter trimming. If you already have trimmed fastqs, you can switch this config entry to `true` and no additional trimming will be performed (still compatible with nextflex_v3 deduplication and structural RNAs filtering for small RNAs).
- `aligned_bams`: When `true` you can directly provide alignment files for ChIP-seq data (either histone modifications or TF). A single SAM or BAM file must be present in the `fastq_path` folder matching the `seq_id` value in the metadata samplefile (same logic than when providing raw fastq file locally). No mapping stats plot will be available when providing bam files this way. Default is `false`.
- Note: These settings are applied to *all* samples in the analysis. If you have some samples to analyze from scratch and other already in an intermediate file: 
	- 1) run the pipeline once with the samples to run from scratch - potentially switching `full_analysis` to `false` for less output. 
	- 2) add the samples you already have intermediate files for in the samplefile and change the corresponding parameters in the config file. 
	- 3) run the pipeline normally again.
	These steps can be repeated if you have raw data, trimmed fastqs and bam files, first creating all the fastq files and then the bam files.

###  Intermediate Target Rules
- `map_only`: Only performs the alignement of all samples. It returns bam files, QC files and mapping metrics.
- `coverage_chip`: Creates bigwig files of coverage for all ChIP samples. The binsize is by default 1bp (can be updated in the config file `chip_tracks: binsize: 1`).

###  Plotting parameters
- `plot_allreps`: When `true`, all individual replicates are shown on heatmaps, metaplots and browsers (can be heavy). When `false` (default), one sample with all merged replicates is used for each sample.

### ChIP Mapping Parameters
- `default`: Standard mapping parameters
- `repeat`: Centromere-specific mapping (more sensitive)
- `repeatall`: Centromere mapping with relaxed MAPQ
- `all`: Relaxed mapping parameters

### DMRs parameters
- By default, DNA methylation data will be analyzed in all sequence contexts (CG, CHG and CHH, where H = A, T or C). The option for CG-only is under development.
- DMRs are called with the R package [DMRcaller](https://www.bioconductor.org/packages/release/bioc/html/DMRcaller.html) (DOI: 10.18129/B9.bioc.DMRcaller) for CG, CHG and CHH and the following (stringent) parameters:
	- CG: `method="noise-filter", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3`
	- CHG: `method="noise_filter", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.2, minGap=200, minSize=50, minReadsPerCytosine=3`
	- CHH: `method="bins", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.1, minGap=200, minSize=50, minReadsPerCytosine=3`
	These parameters were selected based on the most optimal results obtained by the authors [Catoni et al. 2018](https://academic.oup.com/nar/article/46/19/e114/5050634).
- A deeper analysis is available to try different parameters and methods to call the DMRs. Toggle the `use custom_script_dmrs` on the config file to use it. Feel free to edit it as well for different parameters.

##  Additional output options
Below is a list of *cool* outputs that can be generated once whole pipeline ran once. You'll find a basic structure for how to tell snakemake to generate them, feel free to replace the --cores 1 with the HPC profile you would rather use.

### **1. Plotting RNAseq expression levels on target genes**
Given a list of genes (and optional labels), it will plot the expression levels in all the different samples in the samplefile and analysis name defined. Genes uniquely differentially regulated in one sample versus one or more samples are color coded. It is based on a Rdata file created during the Differential Expression analysis.\
To run it, edit the config file with the target gene list file (`rnaseq_target_file`: 1 column list of genes ID that must match the gtf file of the reference genome used, optional second column for gene labels, additional columns can be present but will not be used) and a corresponding label (`rnaseq_target_file_label`: name which will be included in the name of the output pdf) and run the following command, replacing <analysis_name>, <ref_genome> and <rnaseq_target_file_label> with wanted values:
```bash 
snakemake --cores 1 results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.pdf
```
Note that the separators between the variables are two underscores next to each other `__`, except in `plot_expression`.\
An example where <analysis_name>="test_smk" and <ref_genome>="TAIR10", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:
```bash 
snakemake --cores 1 results/RNA/plots/plot_expression__test_smk__TAIR10__my_genes_of_interests.pdf --config rnaseq_target_file="data/target_genes.txt" rnaseq_target_file_label="my_genes_of_interests"
```
Output is a single pdf file named `results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.pdf` where each gene of the list is on an individual page.

### **2. Performing GO analysis on target genes**
Given a file containing a list of genes to do GO analysis on, and optionally a background file (default to all genes in the reference genome), it will perform Gene Ontology analysis.\
By default, GO is not performed since it requires manual input to build a database. To activate it, `GO` needs to be switched to `true` in the config file, and the files to make the GO database should be defined in the config file `gaf_file` and `gene_info_file` below the corresponding reference genome. See [Gene_Ontology.md](Help/Gene_Ontology.md) for more details on how to create the GO database.\
To run it, edit the config file with the target gene list file (`rnaseq_target_file`: 1 column list of genes ID that must match the gtf file of the reference genome used, optional second column for gene labels, additional columns can be present but will not be used) and a corresponding label (`rnaseq_target_file_label`: name which will be included in the name of the output files) and run the following command, replacing <analysis_name>, <ref_genome> and <rnaseq_target_file_label> with wanted values:
```bash 
snakemake --cores 1 results/RNA/GO/TopGO__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.done
```
Note that the separators between the variables are two underscores next to each other `__`.\
An example where <analysis_name>="test_smk" and <ref_genome>="ColCEN", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:
```bash 
snakemake --cores 1 results/RNA/GO/TopGO__test_smk__ColCEN__my_genes_of_interests.done --config rnaseq_target_file="data/target_genes.txt" rnaseq_target_file_label="my_genes_of_interests"
```
Output are two pdf files, one for the biological process terms `results/RNA/plots/topGO_<rnaseq_target_file_label>_BP_treemap.pdf` and one for the molecular function terms `results/RNA/plots/topGO_<rnaseq_target_file_label>_MF_treemap.pdf`. Corresponding tables listing the terms enriched for each gene of the `rnaseq_target_file` are also generated at `results/RNA/GO/topGO_<rnaseq_target_file_label>_<BP|MF>_<GOs|GIDs>.txt` for a focus on the GO terms or the GIDs, respectively.

### **3. Finding motifs on target regions**
Given a bed file containing different regions, it will perform a motifs analysis with meme.\
By default motifs analysis is only performed on the final selected TF peak files (`motifs: true` in the config file). Edit to `allrep: true` in the config file for motifs analysis to be performed on all replicates and pairwise idr peaks if available. A plant motifs database is used by default for tomtom. Download the appropriate file from JASPAR and replace its name in the config file `jaspar_db` and change the `motifs_ref_genome` to match the samples.\
To run the analysis:
```bash 
snakemake --cores 1 results/TF/chkpts/motifs__<motif_target_file_label>.done
```
Note that the separator is two underscores next to each other `__`.\
An example running the pipeline on a slurm hpc, for regions from <ref_genome>="ColCEN", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:
```bash 
snakemake --profile profiles/slurm results/TF/chkpts/motifs__my_regions_of_interests.done --config motifs_target_file="data/target_peaks.txt" motifs_target_file_label="my_regions_of_interests" motifs_ref_genome="ColCEN"
```
Output is the folder `results/TF/<motif_target_file_label>` containing a subdirectory called `meme` and potentially one called `tomtom` with all the results, as described in https://meme-suite.org/meme/index.html. \
When setting `motif_ref_genome`, it is safer to use a reference genome that has already been used in a run. Otherwise, it will be treated like the ref_genome of a sample, creating a fasta file in the genomes/<ref_genome> directory if a fasta file is found at ref_path.\
For the target file chosen `motif_target_file`, if the regions are over 500bp, only the middle 400bp will be used.

### **4. Performing sRNA differential analysis on regions**
Given a bed or gff file, it will perform the small RNA analysis with shortstack followed by differential analysis with edgeR, using all the samples from the sample file but limiting the mapping and counts to the loci in the target file. Edit `srna_target_file` and `srna_target_file_label` in the config file. To run the analysis: 
```bash 
snakemake --cores 1 results/sRNA/clusters/<analysis_name>__<ref_genome>__on_<srna_target_file_label>/Counts.txt
```
Note that the separators between the variables are two underscores next to each other `__` except between `on` and `<srna_target_file_label>` where it's only one `_`.\
An example running the pipeline on a slurm hpc, <analysis_name>="test_smk" and <ref_genome>="ColCEN", while setting the target file and its label "miRNAs" directly in the snakemake command:
```bash 
snakemake --profile profiles/slurm results/sRNA/clusters/test_smk__ColCEN__on_miRNAs/Counts.txt --config sRNA_target_file="data/miRNA.gff" sRNA_target_file_label="miRNAs"
```
Output is the results folder from Shortstack limited to this loci file, followed by the differential cluster analysis with edgeR.

If you only want the results of Shortstack and not the differential analysis, limit the run to the rule `analyze_all_srna_samples_on_target_file` instead, targeting: `results/sRNA/clusters/<analysis_name>__<ref_genome>__on_<srna_target_file_label>/Counts.txt`

The bed or gff file of regions **MUST HAVE** a header with a column called "Name" (the 4th column of a bed file or the 9th column of a gff3).

### **5. Plotting heatmap on regions**
Given a bed file, it will plot a heatmap using deeptools.
Edit `heatmap_target_file` and `heatmap_target_file_label` in the config file. To run the analysis: 
```bash 
snakemake --cores 1 results/combined/plots/Heatmap__<matrix_param>__<env>__<analysis_name>__<ref_genome>__<target_name>.pdf
```
- the <matrix_param> can be `regions` for scaled regions, `tss` for reference point on the TSS or `tes` for reference point on the TES.
- the <env> correspond to the data types to include. Since mC requires different parameters, it has to be done independently. If you have several different data types including mC, and want the order of the regions to be maintained in the mC heatmap based on the other samples, use:
```bash
snakemake --cores 1 results/combined/plots/Heatmap_sorted__<matrix_param>__mC__<analysis_name>__<ref_genome>__<target_name>.pdf
```
This will generate the heatmap for all the other samples first. If you want the regions sorted based on the mC samples only, use:
```bash
snakemake --cores 1 results/combined/plots/Heatmap__<matrix_param>__mC__<analysis_name>__<ref_genome>__<target_name>.pdf
```
To make a heatmap will all the samples (excluding mC), use <env>=`most`. If you want to include mC samples (will probbaly *not work*) use <env>=`all`.

An example running the pipeline on a slurm hpc, with <analysis_name>="test_smk", <ref_genome>="ColCEN", <matrix_param>="regions", on all samples but mC <env>="most", while setting the target file and its label "interesting_genes" directly in the snakemake command:
```bash 
snakemake --profile profiles/slurm results/combined/plots/Heatmap__regions__most__test_smk__ColCEN__interesting_genes.pdf --config heatmap_target_file="data/target_genes.bed" heatmap_target_file_label="interesting_genes"
```
Output is a pdf file, or two if sorted heatmap for mC samples was generated.\
By default, the heatmaps will be scaled by type (i.e. each ChIP mark, each TF, RNAseq, each sRNAseq size and each mC context on an appropriate scale based on the values in the heatmap). It can be changed to "default", where a single scale is used for the whole heatmap, or to "sample" where each sample is scaled individually. This can be changed in the config file `heatmaps_scales`.\
By default, the heatmaps are sorted based on "mean" of all samples accross all regions. This can be changed in the config file `heatmaps_sort_options` to "median" or to "none", keeping the regions in the order of the bedfile.\
If the given bedfile is stranded, the heatmap will be done by splitting the regions into plus and minus strand for proper stranded data (RNAseq and sRNAs) values. If this is not the wanted behavior, disable `stranded_heatmaps` in the config file.\
The color scheme of the heatmap is "seismic" for all samples and "Oranges" for mC. This can be changed manually in the config file `heatmaps_plot_params`.
The size of the scaled regions `middle` (-m in deeptools), the size of the surrounding regions `before` (-b in deeptools) and `after` (-a in deeptools) and the binsize `binsize` (-bs in deeptools) can be edited in the config file in `heatmaps` for each <matrix_params>.

### **6. Plotting metaplot profiles on regions**
Given a bed file, it will plot a metaplot profile using deeptools.
Edit `heatmap_target_file` and `heatmap_target_file_label` in the config file. To run the analysis: 
```bash 
snakemake --cores 1 results/combined/plots/Profile__<matrix_param>__<env>__<analysis_name>__<ref_genome>__<target_name>.pdf
```
Similar to heatmap above for the <matrix_param> options.\
Use <env>="all" to include all samples (mC and others).\
Output is two pdf files, where the samples are grouped by regions or not.\
By default, the heatmaps will be scaled by type (i.e. each ChIP mark, each TF, RNAseq, each sRNAseq size and each mC context on their appropriate scale based on the values in the heatmap). It can be changed to "default", where a single scale is used for the whole heatmap, or to "sample" where each sample is scaled individually. This can be changed in the config file `heatmaps_scales`.
By default, the profiles represent the "mean" accross all regions. This can be changed in the config file `profile_scale` to "median".
By default, the type of plots are "lines". See deeptools documentation for other options and edit `profiles_plot_params` in the config file.
The size of the scaled regions `middle` (-m in deeptools), the size of the surrounding regions `before` (-b in deeptools) and `after` (-a in deeptools) and the binsize `binsize` (-bs in deeptools) can be edited in the config file in `heatmaps` for each <matrix_params>.

### **7. Plotting browser screenshots on regions**
Given a region file, it will plot a browser screenshot using R packages.
Edit `browser_target_file` and `browser_target_file_label` in the config file. To run the analysis: 
```bash 
snakemake --cores 1 results/combined/plots/Browser_<target_name>__<env>__<analysis_name>__<ref_genome>.pdf
```
The target file is a bed-like file, with the following columns: Chr Start End ID Binsize Higlight_starts Higlight_widths\
Each region will be printed individually, and merged into a final PDF.\
Hightlights columns are optional, and correspond to regions of the browser that will be highlighted for this specific region (boxed). As many highlights can be used in a comma-separated lists, the first highlight will be in blue and all the others in red. For example, if the region to plot is chr1:1000-5000, using col6=3000,4000 col7=50,200 will make a blue box higlighting chr1:3000-3050 and a red one highlighting chr1:4000:4200.\
Use <env>="all" to include all samples, "most" for all data-types except mC, or any single environment for data type-specific browsers `[all, most, ChIP, TF, RNA, sRNA, mC]`.\
By default, no TE file is used. If you want to add TE annotations, supply a bed-file in the config file `browser_TE_file`.

### **8. Rerunning a specific analysis**
To rerun a specific analysis, force snakemake to recreate the target file, adding to the snakemake command: `<target_file> --force`
e.g `snakemake --cores 1 results/combined/plots/srna_sizes_stats_test_snakemake_sRNA.pdf --force`
If only the combined analysis is to be performed, and not everything else, delete all the chkpts files in `results/combined/chkpts/` as well as in the chkpt of each relevant environment `results/<env>/chkpts/<env>_analysis__<analysis_name>__<ref_genome>.done`.

## Output Structure

```
epigeneticbutton/
├── config/			# Location for the main config file and recommended location for sample files and target files
├── data/			# Location for test material and examples (e.g. zm_structural_RNAs.fa.gz)
├── Help/			# Location for help files (e.g. Structural_RNAs_Rfam.md)
├── profiles/
│	├── sge/		# Config file to run snakemake on a cluster managed by SGE
│	└── slurm/		# Config file to run snakemake on a cluster managed by SLURM
├── workflow/
│	├── envs/		# Conda environment file for depencies
│	├── rules/		# Snakemake files with data type analysis rules
│	├── scripts/		# R scripts for plots
│	└── snakefile		# main snakefile
├── genomes/			# Genome directories created upon run
│	└── {ref_genome}/	# Reference genome directories with sequence, annotation and indexes
└── results/			# Results directories created upon run
	├── combined/		# Combined analysis results
	│	├── bedfiles/	# Peak calling results
	│	├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
	│	├── logs/	# Log files
	│	├── matrix/	# Data matrices
	│	├── plots/	# Visualization plots
	│	└── reports/	# Analysis reports 
	└── <env>/	# Data type specific directories
		├── chkpts/	# Empty checkpoint files used for pipeline logic. Deleting them will trigger rerunning the corresponding analysis
		├── fastq/	# Processed FASTQ files
		├── logs/	# Log files
		├── mapped/	# Mapped reads (bam)
		├── plots/	# Data type specific plots
		├── reports/	# QC reports
		├── tracks/	# Track files (bigwigs)
		└── */		# data-specific directories (e.g. 'peaks' for ChIP, 'peaks' and 'motifs' for TF, 'DEG' for RNA, 'DMRs' and 'methylcall' for mC, 'clusters' for sRNA)
```

## Known potential issues

1. Relationship between IP and Input for ChIP-seq samples\
Whether a histone ChIP sample is to be compared to H3/H4 or to chromatin input, the sample it is compared to must be called 'Input'. It must also be sequenced either paired-end or single-end but the same than the IPs.

2. small RNA-seq libraries\
Different small RNAseq libraries have different chemistry and might need to be trimmed differently. For now, the code only works if all your samples were done using the same library preparation, either netflex v3 or not. If you have a mix of libraries, you should run the pipeline with each kind separately, and then rerun the analysis with all the samples you want to anlayze together.

3. idr/numpy version\
IDR relies on an older version of numpy to work (due to deprecated np.int) and needs to be loaded as a seperate environment. Not best practice, but more portable than patching idr (np.int=int).

4. Patched ComplexUpset\
Since ggplot2 version 4, the ComplexUpset version on CRAN is not compatible. A patch version exists which is installed from github and works fine for now.

5. Quality-Of-Service slurm configuration\
Due to the time limits on slurm at CSHL, a specific quality of service is used to allow potential long jobs to run for longer. This is likely specific to CSHL cluster. If you want to use slurm and do not have a quality of service setting called "slow_nice" then you can either delete the line `--qos={cluster.qos}` from the `profiles/slurm/config.yaml` file (which might lead to failed runs if you have a time limit), or replace the `qos: "slow_nice"` with another setting that allows longer time limit in the `config/config.yaml` file.

6. Help for local fasq files naming convention\
If using local fastq files for paired-end data, the two read files need to end with `*R1*.f(ast)q(.gz)` and `*R2*.f(ast)q(.gz)`, or `_1.f(ast)q(.gz)` and `_2.f(ast)q(.gz)`, i.e. extensions `fq` or `fastq` and gzipped `.gz` or not. The `<seq_id>` should be common between the two files but distinct from all other files in the same folder (i.e. only one file matching the `*<seq_id>*R1*.f(ast)q(.gz)` expression).

## Features under development
- RAMPAGE
- ATAC-seq

## FAQ

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0) license.

This means you are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material

Under the following terms:
- Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made
- NonCommercial — You may not use the material for commercial purposes without explicit permission
- ShareAlike — If you remix, transform, or build upon the material, you must distribute your contributions under the same license

For commercial use, please contact the author for permission.

See the [LICENSE](LICENSE) file for full details.

## Citation

If you use EpigeneticButton in your research, please cite:

```
Cahn, J., Regulski, M., Lynn, J. et al. MaizeCODE reveals bi-directionally expressed enhancers that harbor molecular signatures of maize domestication. Nat Commun 15, 10854 (2024). https://doi.org/10.1038/s41467-024-55195-w
```

## References 

This pipeline is only a combination of great tools developped by others. A non-exhaustive list of packages used are listed below. Please refer to them for more details.
- [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [ComplexUpset](https://krassowski.github.io/complex-upset/index.html)
- [Conda](https://anaconda.org/anaconda/conda)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [deepTools](https://deeptools.readthedocs.io/en/develop/index.html)
- [DMRcaller](https://bioconductor.org/packages/release/bioc/html/DMRcaller.html)
- [edgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [IDR](https://github.com/nboley/idr)
- [MACS2](https://pypi.org/project/MACS2/#description)
- [Python](https://www.python.org/)
- [R](https://www.r-project.org/)
- [Samtools](https://www.htslib.org/)
- [ShortStack](https://github.com/MikeAxtell/ShortStack)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [SRA-Toolkit](https://github.com/ncbi/sra-tools)
- [STAR](https://github.com/alexdobin/STAR)
- [The MEME suite](https://meme-suite.org/meme/doc/meme-chip.html?man_type=web)
- [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
- [UCSC-GenomeBrowser-kent](https://github.com/ucscGenomeBrowser/kent/)

## Contact

For questions or support, please open an issue in the GitHub repository.
