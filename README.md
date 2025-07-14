# Epigenetic Pipeline for Integrative Chromatin Characterization (EPICC) 
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
   - `data_type`: Type of data [RNAseq | ChIP_* | TF_* | mC | sRNA] (RAMPAGE under development)
   - `line`: Sample line (e.g., B73)
   - `tissue`: Tissue type
   - `sample_type`: Sample identifier
   - `replicate`: Replicate ID
   - `seq_id`: Sequencing ID; use the corresponding SRR####### if downloading from SRA
   - `fastq_path`: Path to FASTQ files; if downloading from SRA, use "SRA" 
   - `paired`: [PE | SE]
   - `ref_genome`: Reference genome name

2. Update `config.yaml` with your paths and parameters:
   - Reference genome path
   - Sample file path
   - Analysis parameters / options
   - Species-specific parameters
   - Resources allocation
   
3. If changing resource allocation for cluster submission, consider adjusting the `profiles/cluster.yaml` for job-specific resources, and the corresponding config file for your cluster scheduler (`profiles/sge/config.yaml` or `profiles/slurm/config.yaml` for SGE or SLURM, respectively). The default is to start 16 jobs maximum in parallel. Keep in mind that units in the cluster file are in MB.

4. Create the `hpclogs` folder for cluster logs if submitting to a hpc:
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

## Sample file configuration

### Common to all types of samples:
- Col2: *line*: Can be any information you want, such as `Col0` or `WT` to annotate and label samples
- Col3: *tissue*: Can be any information you want, such as `leaf` or `mutant` or `6h_stress` to annotate and label samples
The combination line x tissue will be the base for all comparisons (e.g `WT_leaf` vs `WT_roots` or `Col0_control` vs `Ler_stress`)
- Col5: *replicate*: Any value to match the different replicates (e.g Rep1, RepA, 1). All the different replicates are merged for 
- Col6: *seq_id*: Unique identifier to identify the raw data. Can be an SRR number (e.g. SRR27821931) if the data is deposited in SRA, or a unique identifier if the data is local (e.g. `wt_k27`).
- Col7: *fastq_path*: Either `SRA` if raw data to be downloaded from SRA (the SRR number should be used as `seq-id`), or the path to the directory containing the fastq file (e.g. `/archive/fastq`), in which case the `seq_id` should be a unique identifier of the corresponding file in this folder (e.g. `/archive/fastq/raw.reads.wt_k27_R1.fastq.gz`)
- Col8: *paired*: `PE` for paired-end data or `SE` for single-end data. PE samples should have two fastq files R1 and R2 at the location defined above.
- Col9: *ref_genome*: Name of the reference genome to use for mapping (e.g `tair10`). It should be the name of a directory found at the path defined in the config file `ref_path` which contains a single fasta file, a single gff file and a single gtf file. If mapping to multiple references, these directory should be organized in the same `ref_path`. For example, the following structure:
```
/home/
└── genomes/ 		# ref_path: "/home/genomes"
 	├── B73_NAM/	# ref_genome: "B73_NAM" (first ref genome)
	│	├──	B73.fasta	# can be .fa(.gz) or .fasta(.gz)
	│	├──	B73.gff		# can be .gff*(.gz)
	│	└──	B73.gtf		# can be .gtf(.gz)	
 	└── tair10/		# ref_genome: "tair10" (second ref genome)
	 	├──	Ath.fasta	# can be .fa(.gz) or .fasta(.gz)
	 	├──	Ath.gff		# can be .gff*(.gz)
	 	└──	Ath.gtf		# can be .gtf(.gz)
```
The GTF file can be created from a GFF file with cufflinks `gffread -T <gff_file> -o <gtf_file>` and check that `transcript_id` and `gene_id` are correctly assigned in the 9th column. The GFF file should have `gene` in the 3rd column. All files can be gzipped (.gz extension).

### Histone ChIP-seq
- Col1: *data_type*: `ChIP` or `ChIP_<id>` where `<id>` is an identifier to relate an IP sample to its corresponding input. Only necessary in case there are different inputs to be used for different IP samples that otherwise share the same `line` and `tissue` values.
For example: If you have H3K27meac IP samples which you want compared to an H3 sample, and H4K16ac to be compared to H4 samples. Both H3 and H4 samples should be labeled `Input` in sample_type, so to differentiate them, use `ChIP_H3` and `ChIP_H4` for their data_type and the ones of H3K27ac and H4K16ac, respectively. Example:
`ChIP_H3	Col0	WT	H3K27ac	Rep1	wt_k27	./fastq/	PE	Tair10`
`ChIP_H3	Col0	WT	IP	Rep1	wt_h3_ctrl	./fastq/	PE	Tair10`
`ChIP_H4	Col0	WT	H3K27ac	Rep1	wt_h4k16	./fastq/	PE	Tair10`
`ChIP_H4	Col0	WT	IP	Rep1	wt_h4_ctrl	./fastq/	PE	Tair10`
- Col4: *sample_type*: Either `Input` to be used as a control (even if it is actually H3 or IgG pull-down), or the histone mark IP (e.g. H3K9me2). If the mark is not already listed in the config file `chip_callpeaks: peaktype:`, add it to the desired category (either narrow or broad peaks).
- Option: Differential nucleosome sensitivity (DNS-seq) can be analyzed with `ChIP` data_type, using `MNase` for the light digest and `Input` for the heavy digest.

### Transcription factor ChIP
- Col1: *data_type*: `TF_<tf_name>` where `<tf_name>` is the name of the transcirption factor (e.g. for `TB1` data, use `TF_TB1`). This name should be identical for the IP and its input, and for all replicates. Multiple TFs can be analyzed in parallel, each having its own set of IP and Input samples e.g. `TF_<name1>` and `TF_<name2>`.
- Col4: *sample_type*: Either `Input` or `IP`. This works for transcription factors with narrow peaks (default). Use `IPb` for broad peaks.

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
snakemake --cores 1 results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.pdf
```
Note that the separator is two underscores next to each other `__`.
An example where <analysis_name>="test_smk" and <ref_genome>="TAIR10", while setting the target file and its label "my_genes_of_interests" directly in the snakemake command:
```bash 
snakemake --cores 1 results/RNA/plots/plot_expression__test_smk__TAIR10__my_genes_of_interests.pdf --config rnaseq_target_file="data/target_genes.txt" rnaseq_target_file_label="my_genes_of_interests"
```
Output is a single pdf file where each gene of the list is a page, named `results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<rnaseq_target_file_label>.pdf`

2. `rule perform_GO_on_target_file`: Given a file containing a list of genes to do GO analysis on, and a background file (default to all genes in the reference genome), it will perform Gene Ontology analysis. `GO` needs to be switched to `true` in the config file, and either the GO database need to be already made or the files to make it are defined in the config file `gaf_file` and `gene_info_file` below the corresponding reference genome. See `Help_Gene_Ontology` for more details on how to create the GO database. Output are two pdf files, one for the biological process terms `results/RNA/plots/topGO_<rnaseq_target_file_label>_BP_treemap.pdf` and one for the molecular function terms `results/RNA/plots/topGO_<rnaseq_target_file_label>_MF_treemap.pdf`. Corresponding tables listing all the terms enriched for each genes of the `rnaseq_target_file` are also generated at `results/RNA/GO/topGO_<rnaseq_target_file_label>_<BP|MF>_GOs.txt`.

3. Rerunning a specific analysis
To rerun a specific analysis, simply force snakemake to recreate the target file, adding to the snakemake command: `<target_file> --force`
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

## Known potential issues

1. Relationship between IP and Input 
Whether a histone ChIP sample is to be compared to H3/H4 or to chromatin input, the sample it is compared to must be called 'Input'. It must also be sequenced either paired-end or single-end but the same than the IPs.

2. ShortStack version 
The 'epigenetic button' only works with ShortStack v4.0.x version. From v4.1, the developper created a new "condensed" bam format which breaks downstream analysis. New patches could be done in the future for v4.1 compatibility.

## Features under development
- Finishing ChIP-seq and RNA-seq
- Assignment of IP to Input based on suffix (e.g. ChIP_A)
- RAMPAGE
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
