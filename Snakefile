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

# Create dictionaries
refgenome_to_datatype = samples.groupby("ref_genome")["data_type"].unique().to_dict()
datatype_to_samples = samples.groupby("data_type")["sample"].unique().to_dict()
samples_to_replicates = samples.groupby("sample")["replicate"].unique().to_dict()

# Define label for the analysis
analysis_name = config["analysis_name"]

# Create the analysis sample file

import pandas as pd
import os

# Configuration (assuming config is a dictionary with the necessary entries)
config = {
    "sample_file": "path/to/sample_file.txt",
    "ref_path": "path/to/reference/genomes"
}

# Load the sample metadata and perform all operations in a single chain
analysis_samples = (
    samples
    .query("sample != 'Input'") # filter Input samples
    .assign(ref_dir=lambda df: df.apply(lambda row: os.path.join(config["ref_path"], row["ref_genome"]), axis=1)) # create a column with the path to reference genome
    [["data_type", "line", "tissue", "sample", "paired", "ref_dir"]] # select the necessary columns
    .drop_duplicates() # removes replicates
)

# Save the result to 'analysis_samplefile.txt'
analysis_samples.to_csv(f"{config["analysis_name"]}_analysis_samplefile.txt", sep="\t", index=False)

# Define output directories
DIRS = {
    "chkpts": "chkpts",
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
        os.makedirs(f"{data_type}/fastq", exist_ok=True)
        os.makedirs(f"{data_type}/mapped", exist_ok=True)
        os.makedirs(f"{data_type}/tracks", exist_ok=True)
        os.makedirs(f"{data_type}/reports", exist_ok=True)
        os.makedirs(f"{data_type}/logs", exist_ok=True)
        os.makedirs(f"{data_type}/chkpts", exist_ok=True)
        os.makedirs(f"{data_type}/plots", exist_ok=True)
    
    for key, value in dirs.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                os.makedirs(sub_value, exist_ok=True)
        else:
            os.makedirs(value, exist_ok=True)

# Call the function to create directories
create_directories(DATA_TYPES, DIRS)

# Rule all to specify final target
rule all:
	input:
		"chkpts/combined_analysis.done"

# Rule to prepare reference genome for each data type
rule prepare_reference:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
    output:
        chkpt = "chkpts/ref_{ref_genome}_{data_type}.done"
    params:
        ref_path = config["ref_path"]
    log:
        "logs/prepare_ref_{ref_genome}_{data_type}.log"
    conda:
        "envs/reference.yaml"
    shell:
        """
        # Call the original environment preparation script
        {config["scripts_dir"]}/MaizeCode_check_environment.sh \
            -p {params.ref_path} \
            -r {wildcards.ref_genome} \
            -d {wildcards.data_type} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to process samples based on data type
rule process_sample:
    input:
        ref_chkpt = expand("chkpts/ref_{ref_genome}_{data_type}.done", ref_genome = REF_GENOMES, data_type = lambda wildcards: refgenome_to_datatype[wildcards.ref_genome])
    output:
        chkpt = "chkpts/sample_{data_type}_{sample}_{replicate}.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/process_{data_type}_{sample}_{replicate}.log"
    conda:
        "envs/{data_type}_sample.yaml"
    shell:
        """
        # Call the appropriate sample processing script based on data type
        {params.scripts_dir}/MaizeCode_{wildcards.data_type}_sample.sh \
            -s {wildcards.sample} \
            -r {wildcards.replicate} \
            -d {wildcards.data_type} \
            -g {wildcards.ref_genome} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to perform data type specific analysis
rule analyze_sample:
    input:
        ref_chkpt = lambda wildcards: [ 
            f"chkpts/sample_{wildcards.data_type}_{sample}_{replicate}.done"
            for sample in datatype_to_samples.get(wildcards.data_type, [])
            for replicate in samples_to_replicates.get(sample, [])
        ]
    output:
        chkpt = "chkpts/analysis_{data_type}_{analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"],
        analysis_samplefile = f"{config["analysis_name"]}_analysis_samplefile.txt"
    log:
        "logs/analysis_{data_type}_{analysis_name}.log"
    conda:
        "envs/{data_type}_analysis.yaml"
    shell:
        """
        # Call the appropriate analysis script based on data type
        {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            -f {params.analysis_samplefile} > {log} 2>&1
        touch {output.chkpt}
        """

# # Rule to perform ChIP specific analysis
# rule analyze_ChIP:
    # input:
        # ref_chkpt = "chkpts/sample_ChIP_{sample}_{replicate}.done"
    # output:
        # chkpt = "chkpts/analysis_ChIP_{analysis_name}.done"
    # params:
        # scripts_dir = config["scripts_dir"]
    # log:
        # "logs/analysis_ChIP_{analysis_name}.log"
    # conda:
        # "envs/{data_type}_analysis.yaml"
    # shell:
        # """
        # # Call the appropriate analysis script based on data type
        # {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            # -f {input.sample_file} \
            # -p {config["ref_path"]} \
            # -r {wildcards.ref_genome} \
            # -d {wildcards.data_type} > {log} 2>&1
        # touch {output.chkpt}
        # """
		
# # Rule to perform RNA specific analysis
# rule analyze_RNA:
    # input:
        # ref_chkpt = "chkpts/sample_RNA_{sample}_{replicate}.done"
    # output:
        # chkpt = "chkpts/analysis_RNA_{analysis_name}.done"
    # params:
        # scripts_dir = config["scripts_dir"]
    # log:
        # "logs/analysis_RNA_{analysis_name}.log"
    # conda:
        # "envs/{data_type}_analysis.yaml"
    # shell:
        # """
        # # Call the appropriate analysis script based on data type
        # {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            # -f {input.sample_file} \
            # -p {config["ref_path"]} \
            # -r {wildcards.ref_genome} \
            # -d {wildcards.data_type} > {log} 2>&1
        # touch {output.chkpt}
        # """
		
# # Rule to perform TF specific analysis
# rule analyze_TF:
    # input:
        # ref_chkpt = "chkpts/sample_TF_{sample}_{replicate}.done"
    # output:
        # chkpt = "chkpts/analysis_TF_{analysis_name}.done"
    # params:
        # scripts_dir = config["scripts_dir"]
    # log:
        # "logs/analysis_TF_{analysis_name}.log"
    # conda:
        # "envs/{data_type}_analysis.yaml"
    # shell:
        # """
        # # Call the appropriate analysis script based on data type
        # {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            # -f {input.sample_file} \
            # -p {config["ref_path"]} \
            # -r {wildcards.ref_genome} \
            # -d {wildcards.data_type} > {log} 2>&1
        # touch {output.chkpt}
        # """

# # Rule to perform mC specific analysis
# rule analyze_mC:
    # input:
        # ref_chkpt = "chkpts/sample_mC_{sample}_{replicate}.done"
    # output:
        # chkpt = "chkpts/analysis_mC_{analysis_name}.done"
    # params:
        # scripts_dir = config["scripts_dir"]
    # log:
        # "logs/analysis_mC_{analysis_name}.log"
    # conda:
        # "envs/{data_type}_analysis.yaml"
    # shell:
        # """
        # # Call the appropriate analysis script based on data type
        # {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            # -f {input.sample_file} \
            # -p {config["ref_path"]} \
            # -r {wildcards.ref_genome} \
            # -d {wildcards.data_type} > {log} 2>&1
        # touch {output.chkpt}
        # """

		
# Rule to perform combined analysis
rule combined_analysis:
    input:
        ref_chkpt = expand("chkpts/analysis_{data_type}_{analysis_name}.done", data_type=DATA_TYPES, analysis_name=analysis_name)
    output:
        chkpt = "chkpts/combined_analysis.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/combined_analysis.log"
    conda:
        "envs/combined.yaml"
    shell:
        """
        # Call the combined analysis script
        {params.scripts_dir}/MaizeCode_analysis.sh \
            -f {input.sample_file} \
            -p {config["ref_path"]} > {log} 2>&1
        touch {output.chkpt}
        """ 
