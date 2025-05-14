import pandas as pd
import numpy as np
import os
import fnmatch
import sys
import re
from snakemake.utils import min_version

min_version("6.0")

# Configuration
configfile: "config.yaml"

# Define patterns and envs
env_patterns = [
    ("ChIP*", "ChIP"),
    ("RNAseq", "RNA"),
    ("RAMPAGE", "RNA"),
    ("shRNA", "shRNA"),
    ("mC", "mC"),
    ("TF_*", "TF"),
]

# Define sample metadata
samples = pd.read_csv(config["sample_file"], sep="\t", header=None,
                      names=["data_type", "line", "tissue", "sample_type", "replicate", 
                             "seq_id", "fastq_path", "paired", "ref_genome"])

# Create an "env" column based on the string patterns in data_type
conditions = [samples["data_type"].str.match(fnmatch.translate(pat)) for pat, _ in env_patterns]
choices = [env for _, env in env_patterns]
samples["env"] = np.select(conditions, choices, default="unknown")
    
# Check for unknown envs and exit if any
unknowns = samples.loc[samples["env"] == "unknown", "data_type"].unique()
if len(unknowns) > 0:
    unknown_list = "\n".join(f"  - {dt}" for dt in unknowns)
    raise ValueError(f"Type of data unknown for the following data types:\n"
                     f"{unknown_list}\n\n"
                     "Please check your sample file or update the env_patterns.")

# Create sample_name column
def create_sample_name(row):
    return f"{row['data_type']}__{row['line']}__{row['tissue']}__{row['sample_type']}__{row['replicate']}__{row['ref_genome']}"

samples["sample_name"] = samples.apply(create_sample_name, axis=1)

# Create a dictionary to store the information for each sample
sample_info_map = {
    (row["data_type"], row["line"], row["tissue"], row["sample_type"], row["replicate"], row["ref_genome"]): {
        "seq_id": row["seq_id"],
        "fastq_path": row["fastq_path"],
        "paired": row["paired"],
        "env": row["env"],
        "sample_name": row["sample_name"]
    }
    for _, row in samples.iterrows()
}

# Define reference genomes and the path to them
REF_GENOMES = set(samples["ref_genome"].unique())
REF_PATH = config["ref_path"]

# Define data types
DATA_TYPES = set(samples["data_type"].unique())

# Define the folde rin which the snakemake pipeline has been cloned
REPO_FOLDER = config["repo_folder"]

# Define label for the analysis
analysis_name = config["analysis_name"]

# Function to split the sample_name to recover its components
def parse_sample_name(sample_name):
    data_type, line, tissue, sample_type, rep, ref_genome = sample_name.split("__")
    parsed = {
        "data_type": data_type,
        "line": line,
        "tissue": tissue,
        "sample_type": sample_type,
        "replicate": rep,
        "ref_genome": ref_genome
    }    
    
    # To extract ChIP Input group from data_type
    if data_type.startswith("ChIP_"):
        chip_parts = data_type.split("_", 1)
        if len(chip_parts) == 2:
            parsed["group_label"] = chip_parts[1]   # e.g., "A"
    # To extract TF name from data_type
    if data_type.startswith("TF_"):
        tf_parts = data_type.split("_", 1)
        if len(tf_parts) == 2:
            parsed["tf_name"] = tf_parts[1]   # e.g., "TB1"

    return parsed

# Function to create a unique name for each sample based on the sample columns, and later based on wildcards
def sample_name(d):
    return f"{d['data_type']}__{d['line']}__{d['tissue']}__{d['sample_type']}__{d['replicate']}__{d['ref_genome']}"

# Function to access extra information form the samplefile using wildcards
def get_sample_info(wildcards, field):
    key = (wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.replicate, wildcards.ref_genome)
    return sample_info_map[key][field]

# Function to access extra information form the samplefile using the name
def get_sample_info_from_name(sname, samples, field):
    match = samples.loc[samples["sample_name"] == sname]
    if match.empty:
        raise ValueError(f"\nSample '{sname}' not found in {samples} table.")
    else:
        return match[field].iloc[0]

# Function to extract all samples of each data_type based on sample name
def get_sample_names_by_env(data_type):
    sample_names = samples.loc[samples['env'] == data_type, "sample_name"].tolist()
    return sample_names

# Map data types to environments
datatype_to_env = dict(zip(samples["data_type"], samples["env"]))

# Get unique list of environments
UNIQUE_ENVS = samples["env"].unique().tolist()

# Load the sample metadata and perform all operations in a single chain
analysis_samples = (
    samples
    .query("sample_type != 'Input'") # filter Input samples
    [["data_type", "line", "tissue", "sample_type", "ref_genome", "paired"]]
    .drop_duplicates()
)

# Save the result to 'analysis_samplefile.txt'
analysis_samples.to_csv(f"{analysis_name}__analysis_samplefile.txt", sep="\t", index=False)

# Function to create the unique name for each sample from analysis file
def sample_name_analysis(d):
    return f"{d['data_type']}__{d['line']}__{d['tissue']}__{d['sample_type']}__{d['ref_genome']}"

# To assign all replicates for each ChiP input
chip_input_to_replicates = (
    samples
    .query("sample_type == 'Input'")
    .groupby(["data_type", "line", "tissue", "ref_genome"])["replicate"]
    .apply(list)
    .to_dict()
)

# To assign all replicates to each analysis samples
analysis_to_replicates = (
    samples
    .groupby(["data_type", "line", "tissue", "sample_type", "ref_genome"])["replicate"]
    .apply(list)
    .to_dict()
)

# To assign all ChiP IP to their inputs
chip_samples_to_replicates = (
    analysis_samples
    .groupby(["data_type", "line", "tissue", "ref_genome"])["sample_type"]
    .apply(list)
    .to_dict()
)

# Define output directories
DIRS = {
    "genomes": "genomes",
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
def create_directories(unique_envs, dirs):
    for env in unique_envs:
        for d in ["fastq", "mapped", "tracks", "reports", "logs", "chkpts", "plots"]:
            os.makedirs(f"{env}/{d}", exist_ok=True)
        if env in ["ChIP", "TF"]:
            os.makedirs(f"{env}/peaks", exist_ok=True)
    
    for key, value in dirs.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                os.makedirs(sub_value, exist_ok=True)
        else:
            os.makedirs(value, exist_ok=True)

# Call the function to create directories
create_directories(UNIQUE_ENVS, DIRS)

# Include all rule files
include: "rules/environment_setup.smk"
include: "rules/sample_download.smk"
include: "rules/ChIPseq.smk"
include: "rules/RNAseq.smk"

# Rule all to specify final target
rule all:
	input:
		f"chkpts/combined_analysis__{analysis_name}.done"

# Rule to specify final target if only mapping is required
rule map_only:
    input:
        [
            f"{env}/chkpts/process__{sample_name}.done"
            for env in UNIQUE_ENVS
            for sample_name in get_sample_names_by_env(env)
        ]

rule coverage_chip:
    input: 
        [
            f"ChIP/tracks/coverage__{sample_name}.bw"
            for sample_name in get_sample_names_by_env("ChIP")
        ]

# Rule to perform combined analysis
rule combined_analysis:
    input:
        expand("ChIP/chkpts/ChIP_analysis__{ref_genome}.done", ref_genome=REF_GENOMES),
        expand("RNA/chkpts/RNA_analysis__{ref_genome}.done", ref_genome=REF_GENOMES),
        expand("chkpts/ref__{ref_genome}.done", ref_genome=REF_GENOMES)
    output:
        chkpt = f"chkpts/combined_analysis__{analysis_name}.done"
    params:
        region_file="all_genes.txt",
        scripts_dir = os.path.join(REPO_FOLDER,"scripts"),
        analysis_samplefile = f"{analysis_name}__analysis_samplefile.txt"
    log:
        f"logs/combined_analysis__{analysis_name}.log"
    conda:
        "envs/combined.yaml"
    shell:
        """
        touch {output.chkpt}
        """ 
