# EpigeneticButton

A Snakemake-based pipeline for analyzing and integrating various types of genomics datasets, including ChIP-seq, RNA-seq, RAMPAGE, shRNA, mC, and TF ChIP-seq data.

## Overview

EpigeneticButton is a comprehensive pipeline that processes and analyzes multiple types of genomics data. It provides an automated workflow for:
- Data preprocessing and quality control
- Read mapping and alignment
- Peak calling and differential expression analysis
- Data integration and visualization
- Transposable element analysis

## Features

- **Multiple Data Types Support**:
  - ChIP-seq
  - RNA-seq
  - RAMPAGE
  - shRNA
  - MethylC-seq (mC)
  - Transcription Factor ChIP-seq (TF)

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
```bash
cd epigeneticbutton
```

2. Install snakemake:
```bash
conda install -c bioconda snakemake
```

## Usage

### Configuration

1. Prepare your sample metadata file (`samples.tsv`) with the following columns:
   - `data_type`: Type of data [RNAseq | RAMPAGE | shRNA | ChIP_* | TF_* | mC]
   - `line`: Sample line (e.g., B73)
   - `tissue`: Tissue type
   - `sample`: Sample identifier
   - `replicate`: Replicate ID
   - `seq_id`: Sequencing ID
   - `fastq_path`: Path to FASTQ files
   - `paired`: [PE | SE]
   - `ref_genome`: Reference genome name

2. Update `config.yaml` with your paths and parameters:
   - Reference genome path
   - Sample file path
   - Analysis parameters
   - Resource allocation

### Running the Pipeline

1. To run the pipeline:
```bash
snakemake --use-conda --conda-frontend mamba --cores 4
```

2. To run the pipeline on a HPC using qsub:
```bash
qsub -V -cwd -pe threads 20 -l m_mem_free=8G -l tmp_free=16G -o pipeline.log -j y -N epigeneticbutton <<-'EOF'
#!/bin/bash
snakemake --use-conda --conda-frontend mamba --cores 20
EOF
```

This qsub command:
- Requests 20 threads
- Allocates 8GB of memory
- Allocates 16GB of temporary storage
- Outputs logs to pipeline.log
- Uses the conda environment
- Runs snakemake with 20 cores

3. Optional: for increased speed for solving environments, consider using mamba and/or prebuilding the environments:
```bash
conda install mamba -n base -c conda-forge # to install mamba, if not already
```
```bash
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only --cores 1
```
It is also recommended to set strict conda channel priorities:
```bash
conda config --set channel_priority strict
```

4. Optional: to test the pipeline, consider generating a DAG first to make sure your sample files and parameters work:
```bash
snakemake --dag | dot -Tpng > dag.png
```
or to force all steps to be performed:
```bash
 snakemake --dag --forceall | dot -Tpng > dag.png
```

### Output Structure

```
epigeneticbutton/
├── chkpts/              # Checkpoint files
├── combined/            # Combined analysis results
│   ├── peaks/          # Peak calling results
│   ├── DEG/            # Differential expression results
│   ├── TSS/            # Transcription start site analysis
│   ├── reports/        # Analysis reports
│   ├── matrix/         # Data matrices
│   └── plots/          # Visualization plots
├── {data_type}/        # Data type specific directories
│   ├── fastq/          # Processed FASTQ files
│   ├── mapped/         # Mapped reads
│   ├── tracks/         # Track files
│   ├── reports/        # QC reports
│   └── plots/          # Data type specific plots
└── logs/               # Log files
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

### Analysis Parameters
- `mark_of_interest`: Default histone mark (e.g., H3K27ac)
- `perform_analysis`: Enable/disable complete analysis
- `perform_combined`: Enable/disable combined analysis
- `perform_heatmaps`: Enable/disable heatmap generation
- `perform_te_analysis`: Enable/disable TE analysis

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
