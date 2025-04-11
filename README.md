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
cd epigeneticbutton
```

2. Install dependencies using conda:
```bash
conda env create -f environment.yaml
conda activate epigeneticbutton
```

For faster environment resolution, you can use mamba instead of conda:
```bash
# Install mamba in your base environment
conda install -n base -c conda-forge mamba

# Create the environment using mamba
mamba env create -f environment.yaml
mamba activate epigeneticbutton
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

```bash
snakemake --cores <N> --configfile config.yaml
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
