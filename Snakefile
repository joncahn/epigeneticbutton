import pandas as pd
import os
from snakemake.utils import min_version

min_version("6.0")

# Configuration
configfile: "config.yaml"

# Define sample metadata
samples = pd.read_csv(config["sample_file"], sep="\t", header=None,
                      names=["data_type", "line", "tissue", "sample", "replicate", 
                             "seq_id", "fastq_path", "paired", "ref_genome"])

# Define reference genomes
REF_GENOMES = set(samples["ref_genome"].unique())

# Define data types
DATA_TYPES = set(samples["data_type"].unique())

# Define label for the analysis
analysis_name = config["analysis_name"]

# Load the sample metadata and perform all operations in a single chain
analysis_samples = (
    samples
    .query("sample != 'Input'") # filter Input samples
    .assign(ref_dir=lambda df: df.apply(lambda row: os.path.join(config["ref_path"], row["ref_genome"]), axis=1)) # create a column with the path to reference genome
    [["data_type", "line", "tissue", "sample", "paired", "ref_dir"]] # select the necessary columns
    .drop_duplicates() # removes replicates
)

# Save the result to 'analysis_samplefile.txt'
analysis_samples.to_csv(f"{analysis_name}__analysis_samplefile.txt", sep="\t", index=False)

# Create dictionaries
refgenome_to_datatype = samples.groupby("ref_genome")["data_type"].unique().to_dict()
samples_to_replicates = samples.groupby("sample")["replicate"].unique().to_dict()
datatype_to_samples = analysis_samples.groupby("data_type")["sample"].unique().to_dict()

# Define output directories
DIRS = {
    "chkpts": "chkpts",
    "logs": "logs",
    "combined": {
        "peaks": "combined/peaks",
        "DEG": "combined/DEG",
        "TSS": "combined/TSS",
        "reports": "combined/reports",
        "matrix": "combined/matrix",
        "plots": "combined/plots",
        "chkpts": "combined/chkpts",
        "logs": "combined/logs"
    }
}

# Function to create directories
def create_directories(data_types, dirs):
    for data_type in data_types:
        for d in ["fastq", "mapped", "tracks", "reports", "logs", "chkpts", "plots"]:
            os.makedirs(f"{data_type}/{d}", exist_ok=True)
    
    for key, value in dirs.items():
        if isinstance(value, dict):
            for sub_value in value.values():
                os.makedirs(sub_value, exist_ok=True)
        else:
            os.makedirs(value, exist_ok=True)

# Call the function to create directories
create_directories(DATA_TYPES, DIRS)

# Rule all to specify final target
rule all:
	input:
		f"chkpts/combined_analysis__{analysis_name}.done"

# Rule to prepare reference genome for each data type
rule prepare_reference:
    input:
        refs = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
    output:
        chkpt = "chkpts/ref__{ref_genome}__{data_type}.done"
    params:
        ref_path = config["ref_path"],
        scripts_dir = config["scripts_dir"]
    log:
        "logs/prepare_ref__{ref_genome}__{data_type}.log"
    conda:
        "envs/reference.yaml"
    shell:
        """
        # Call the original environment preparation script
        qsub {params.scripts_dir}/MaizeCode_check_environment.sh \
            -p {params.ref_path} \
            -r {wildcards.ref_genome} \
            -d {wildcards.data_type} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to process samples based on data type
rule process_sample:
    input:
        ref_chkpt = "chkpts/ref__{ref_genome}__{data_type}.done"
    output:
        chkpt = "chkpts/sample__{data_type}__{sample}__{replicate}__{ref_genome}.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/process__{data_type}__{sample}__{replicate}__{ref_genome}.log"
    conda:
        "envs/{data_type}_sample.yaml"
    shell:
        """
        # Call the appropriate sample processing script based on data type
        qsub {params.scripts_dir}/MaizeCode_{wildcards.data_type}_sample.sh \
            -s {wildcards.sample} \
            -r {wildcards.replicate} \
            -d {wildcards.data_type} \
            -g {wildcards.ref_genome} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to perform data type specific analysis
rule analyze_sample:
    input:
        process_chkpt = lambda wildcards: expand(
            "chkpts/sample__{data_type}__{sample}__{replicate}__{ref_genome}.done",
            data_type = [wildcards.data_type],
            sample = datatype_to_samples[wildcards.data_type],
            replicate = [
                rep 
                for sample in datatype_to_samples[wildcards.data_type]
                for rep in samples_to_replicates[sample]
            ],
            ref_genome = [
                samples[samples["sample"] == sample]["ref_genome"].iloc[0]
                for sample in datatype_to_samples[wildcards.data_type]
            ]
        )
    output:
        chkpt = "chkpts/analysis__{data_type}__{analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"],
        analysis_samplefile = f"{analysis_name}__analysis_samplefile.txt"
    log:
        "logs/analysis__{data_type}__{analysis_name}.log"
    conda:
        "envs/{data_type}_analysis.yaml"
    shell:
        """
        # Call the appropriate analysis script based on data type
        qsub {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            -f {params.analysis_samplefile} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to perform combined analysis
rule combined_analysis:
    input:
        analysis_chkpt = expand("chkpts/analysis__{data_type}__{analysis_name}.done", 
            data_type = DATA_TYPES, 
            analysis_name = analysis_name)
    output:
        chkpt = f"chkpts/combined_analysis__{analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"],
        ref_path = config["ref_path"],
        analysis_samplefile = f"{analysis_name}__analysis_samplefile.txt"
    log:
        f"logs/combined_analysis__{analysis_name}.log"
    conda:
        "envs/combined.yaml"
    shell:
        """
        # Call the combined analysis script
        qsub {params.scripts_dir}/MaizeCode_analysis.sh \
            -f {params.analysis_samplefile} \
            -p {params.ref_path} > {log} 2>&1
        touch {output.chkpt}
        """ 
        
