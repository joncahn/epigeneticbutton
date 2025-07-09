# Epigenetic Pipeline for Integrative Characterization (EPIC) 
# aka Epigenetic Button

A Snakemake-based pipeline for analyzing and integrating various types of (epi)genomics datasets, including histone and transcription factor ChIP-seq, RNA-seq, RAMPAGE, small RNA-seq, and methylC-seq.

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
  - RAMPAGE
  - small RNA-seq
  - MethylC-seq (mC)

- **Automated Analysis**:
  - Reference genome preparation
  - Sample-specific processing
  - Data type-specific analysis
  - Combined analysis across samples
  - Quality control and reporting

- **Flexible Configuration**:
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

2. Install snakemake:
```bash
conda install -c bioconda snakemake
```

## Usage

### Configuration

1. Prepare your sample metadata file (default to `samples.tsv`) with the required columns below (see Input requirements for more details specific to each data-type):
   - `data_type`: Type of data [RNAseq | ChIP_* | TF_* | mC | sRNA] (TF and RAMPAGE under development)
   - `line`: Sample line (e.g., B73)
   - `tissue`: Tissue type
   - `sample`: Sample identifier
   - `replicate`: Replicate ID
   - `seq_id`: Sequencing ID; use the corresponding SRR####### if downloading from SRA
   - `fastq_path`: Path to FASTQ files; if downloading from SRA, use "SRA" 
   - `paired`: [PE | SE]
   - `ref_genome`: Reference genome name

2. Update `config.yaml` with your paths and parameters:
   - Reference genome path
   - Sample file path
   - Analysis parameters / options
   - Species-specific parameters (e.g. genome size)
   - Resources allocation
   
3. If changing resource allocation for cluster submission, consider adjusting the `profiles/cluster.yaml` for job-specific resources, and the corresponding config file for your cluster scheduler (`profiles/sge/config.yaml` or `profiles/slurm/config.yaml` for SGE or SLURM, respectively). The default is to use 96 threads maximum in parallel. Keep in mind that units in the cluster file are in MB.

4. Create the `hpclogs` folder for cluster logs:
```bash
mkdir -p hpclogs
```

### Running the Pipeline

1. To run the pipeline locally:
```bash
snakemake --use-conda --conda-frontend conda --cores 12
```

2. To run the pipeline on a HPC-SGE (using qsub):
```bash
snakemake --profile profiles/sge
```

3. To run the pipeline on a HPC-slurm (using sbatch):
```bash
snakemake --profile profiles/slurm
```

4. Optional: consider prebuilding the environments to make sure no conflict arise (it should take ~5min):
```bash
snakemake --use-conda --conda-frontend conda --conda-create-envs-only --cores 1
```

5. Optional: to test the pipeline, consider generating a DAG first to make sure your sample files and parameters work:
```bash
snakemake --dag | dot -Tpng > dag.png
```
or to force all steps to be performed:
```bash
 snakemake --dag --forceall | dot -Tpng > dag.png
```

*For full understanding of snakemake capabilities and option: https://snakemake.readthedocs.io/en/stable/*

## Input requirements

### Transcription factor ChIP
- Set the sample_type to `IP` for transcription factors with narrow peaks (default); use `IPb` for broad peaks.

## Output Structure

```
epigeneticbutton/
├── config/ 		   # Location for the main config file and recommended location for sample files and target files
├── data/			   # Location for test material and examples (e.g. zm_structural_RNAs.fa.gz)
├── Help/			   # Location for help files (e.g. Help_structural_RNAs_database_with_Rfam)
├── profiles/
│	├── sge/		   # Config file to run snakemake on a cluster managed by SGE
│	├── slurm/		   # Config file to run snakemake on a cluster managed by SLURM
│	└──	cluster.yaml   # Config file with job-specific resources for cluster submission (used by sge and slurm)
├── workflow/
│	├── envs/		   # Conda environment file for depencies
│	├── rules/		   # Snakemake files with data type analysis rules
│	├── scripts/	   # R scripts for plots
│	└── snakefile	   # main snakefile
├── genomes/		   # Genome directories created upon run
│	└──	{ref_genome}/  # Genome-specific directories with sequence, annotation and indexes
└── results/           # Results directories created upon run
	├── combined/      # Combined analysis results
	│	├── logs/      # Log files
	│   ├── chkpts/    # Peak calling results
	│   ├── peaks/     # Peak calling results
	│   ├── DEG/       # Differential expression results
	│   ├── TSS/       # Transcription start site analysis
	│   ├── reports/   # Analysis reports
	│   ├── matrix/    # Data matrices
	│   └── plots/     # Visualization plots
	└── {data_type}/   # Data type specific directories
	    ├── fastq/     # Processed FASTQ files
	    ├── mapped/    # Mapped reads (bam)
	    ├── tracks/    # Track files (bigwigs)
	    ├── reports/   # QC reports
	    ├── */		   # data-specific directories (e.g. 'dmrs' for mC, 'clusters' for sRNA)
		└── plots/     # Data type specific plots
```

## Configuration Options

### Mapping Parameters
- `default`: Standard mapping parameters
- `colcen`: Centromere-specific mapping (more sensitive)
- `colcenall`: Centromere mapping with relaxed MAPQ
- `all`: Relaxed mapping parameters

###  Intermediate Target Rules
- `map_only`: Only performs the up to mapping of all samples. It returns bam files, QC files and mapping metrics.
- `coverage_chip`: Creates bigwig files of coverage for all ChIP samples. The binsize is by default 1bp (can be updated in config (chip_tracks: binsize: 1)

###  Additional output options
1. `rule plot_expression_levels`: Given a list of genes (and optional labels), it will plot the expression levels in all the different samples in the samplefile and analysis name defined. Genes uniquely differentially regulated in one sample versus one or more samples are color coded. It is based on a Rdata file created during the Differential Expression analysis (rule call_all_degs). To use it, edit the config file with the target gene list file (`rnaseq_target_file`: 1 column list of genes ID that must match the gtf file of the reference genome used, optional second column for gene labels, additional columns can be present but will not be used) and a corresponding label (`rnaseq_target_file_label`) and run the following command, replacing {analysis_name}, {ref_genome} and {target_label} with wanted values:
```bash 
snakemake --cores 1 results/RNA/plots/plot_expression__{analysis_name}__{ref_genome}__{rnaseq_target_file_label}.pdf
```
An example where {analysis_name}="test_smk" and {ref_genome}="TAIR10", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:
```bash 
snakemake --cores 1 results/RNA/plots/plot_expression__test_smk__TAIR10__my_genes_of_interests.pdf --config rnaseq_target_file="data/target_genes.txt" rnaseq_target_file_label="my_genes_of_interests"
```
Output is a single pdf file where each gene of the list is a page, named `results/RNA/plots/plot_expression__{analysis_name}__{ref_genome}__{rnaseq_target_file_label}.pdf`

2. Rerunning a specific analysis
To rerun a specific analysis, simply force snakemake to recreate the target file, adding to the snakemake command: `{target_file} --force`
e.g `snakemake --cores 1 results/combined/plots/srna_sizes_stats_test_snakemake_sRNA.pdf --force`

### DMRs parameters
- By default, DNA methylation data will be analyzed in all sequence contexts (CG, CHG and CHH, where H = A, T or C). The option for CG-only is under development.
- DMRs are called with the R package DMRcaller (DOI: 10.18129/B9.bioc.DMRcaller) for CG and CHH and the following (stringent) parameters: 
	CG: `method="noise-filter", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3`
	CHG: `method="noise_filter", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.2, minGap=200, minSize=50, minReadsPerCytosine=3`
	CHH: `method="bins", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.1, minGap=200, minSize=50, minReadsPerCytosine=3`
- Modify the script `scripts/R_call_DMRs.R` if other paramteres/contexts should be performed, or make a copy such as `scripts/R_call_DMRs_custom.R` and replace it in the `call_DMRs_pairwise` rule in the `mC.smk` file.

### Analysis Parameters
- `mark_of_interest`: Default histone mark (e.g., H3K27ac)
- `perform_analysis`: Enable/disable complete analysis
- `perform_combined`: Enable/disable combined analysis
- `perform_heatmaps`: Enable/disable heatmap generation
- `perform_te_analysis`: Enable/disable TE analysis

### Known potential issues

1. Relationship between IP and Input 
Whether a histone ChIP sample is to be compared to H3/H4 or to chromatin input, the sample it is compared to must be called 'Input'. It must also be sequenced either paired-end or single-end but the same than the IPs.

2. ShortStack version 
The 'epigenetic button' only works with ShortStack v4.0.x version. From v4.1, the developper created a new "condensed" bam format which breaks downstream analysis. New patches could be done in the future for v4.1 compatibility.

### Features under development
- Finishing ChIP-seq and RNA-seq
- Assignment of IP to Input based on suffix (e.g. ChIP_A)
- RAMPAGE
- small RNAseq
- Browser: create a hub/jbrowse session? invert minus stranded bigwigs.
- Plotting
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

## Contact

For questions or support, please open an issue in the GitHub repository.
