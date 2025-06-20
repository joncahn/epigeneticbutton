import pandas as pd
import numpy as np
import os
import fnmatch
import sys
import re
from snakemake.utils import min_version
from itertools import combinations

min_version("6.0")

# rule show_env:
    # shell:
        # "echo Conda version: $(conda --version) && "
        # "echo $(conda config --show solver) &&"
        # "echo Conda executable: $(which conda) && "
        # "echo Environment: $CONDA_DEFAULT_ENV && "
        # "conda info"

CONDA_ENV=os.path.join(REPO_FOLDER,"envs/epibutton.yaml")

# Configuration
configfile: "config.yaml"
analysis = config['full_analysis']

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

# Function to create a unique name for each sample based on the sample columns, and later based on wildcards
def sample_name_str(d, string):
    if string == 'sample':
        return f"{d['data_type']}__{d['line']}__{d['tissue']}__{d['sample_type']}__{d['replicate']}__{d['ref_genome']}"
    elif string == 'analysis':
        return f"{d['data_type']}__{d['line']}__{d['tissue']}__{d['sample_type']}__{d['ref_genome']}"

# Add a sample_name column to the sample file
samples["sample_name"] = samples.apply(lambda row: sample_name_str(row, 'sample'), axis=1)

# Function to extract extra information from data_Type (either the name of the TF, or the ChIP Input grouping if it exists)
def extract_extra_info(d):
    parts=d['data_type'].split("_", 1)
    if len(parts) == 2 and parts[0] in ["ChIP", "TF"]:
        return parts[1]
    else:
        return "N/A"   
    
# Add an extra_info column to the sample file
samples["extra_info"] = samples.apply(extract_extra_info, axis=1)

# Create a dictionary to store the information for each sample
sample_info_map = {
    (row["data_type"], row["line"], row["tissue"], row["sample_type"], row["replicate"], row["ref_genome"]): {
        "seq_id": row["seq_id"],
        "fastq_path": row["fastq_path"],
        "paired": row["paired"],
        "env": row["env"],
        "sample_name": row["sample_name"],
        "extra_info": row["extra_info"]
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

# Function to access extra information form the samplefile using the name
def get_sample_info_from_name(sname, samples, field):
    match = samples.loc[samples["sample_name"] == sname]
    if match.empty:
        raise ValueError(f"\nSample '{sname}' not found in {samples} table.")
    else:
        return match[field].iloc[0]

# Function to extract all samples of each env based on sample name
def get_sample_names_by_env(env, samples):
    sample_names = samples.loc[samples['env'] == env, "sample_name"].tolist()
    return sample_names

# Get unique list of environments
UNIQUE_ENVS = samples["env"].unique().tolist()

# Load the sample metadata and perform all operations in a single chain
analysis_samples = (
    samples
    .query("sample_type != 'Input'") # filter Input samples
    [["env", "data_type", "line", "tissue", "sample_type", "ref_genome", "paired"]]
    .drop_duplicates()
)
# Add a sample_name column to the analysis_sample file
analysis_samples["sample_name"] = analysis_samples.apply(lambda row: sample_name_str(row, 'analysis'), axis=1)

# Save the result to 'analysis_samplefile.txt'
analysis_samples.to_csv(f"{analysis_name}__analysis_samplefile.txt", sep="\t", index=False)

# To assign all replicates to each sample
analysis_to_replicates = (
    samples
    .groupby(["data_type", "line", "tissue", "sample_type", "ref_genome"])["replicate"]
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
        if env in ["mC"]:
            os.makedirs(f"{env}/methylcall", exist_ok=True)
            os.makedirs(f"{env}/DMRs", exist_ok=True)
        if env in ["RNA"]:
            os.makedirs(f"{env}/DEG", exist_ok=True)
    
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
include: "rules/mC.smk"
include: "rules/plotting_with_R.smk"

# Rule all to specify final target
rule all:
	input:
		f"chkpts/combined_analysis__{analysis_name}.done"

### Intermediate target rules
# Rule to specify final target if only mapping is required
rule map_only:
    input:
        expand("combined/plots/mapping_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=UNIQUE_ENVS)

# Rule to specify final target if only chip coverage is wanted
rule coverage_chip:
    input: 
        [
            f"ChIP/tracks/coverage__{sample_name}.bw"
            for sample_name in get_sample_names_by_env("ChIP", samples)
        ]

# Rule to perform combined analysis
rule combined_analysis:
    input:
        expand("ChIP/chkpts/ChIP_analysis__{ref_genome}.done", ref_genome=REF_GENOMES if "ChIP" in UNIQUE_ENVS else []),
        expand("TF/chkpts/ChIP_analysis__{ref_genome}.done", ref_genome=REF_GENOMES if "TF" in UNIQUE_ENVS else []),
        expand("RNA/chkpts/RNA_analysis__{ref_genome}.done", ref_genome=REF_GENOMES if "RNA" in UNIQUE_ENVS else []),
        expand("mC/chkpts/mC_analysis__{ref_genome}.done", ref_genome=REF_GENOMES if "mC" in UNIQUE_ENVS else []),
        expand("combined/plots/mapping_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","TF","RNA","mC"]]),
        expand("combined/plots/peak_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","TF"]])
    output:
        chkpt = f"chkpts/combined_analysis__{analysis_name}.done"
    params:
        region_file="all_genes.txt",
        scripts_dir = os.path.join(REPO_FOLDER,"scripts"),
        analysis_samplefile = f"{analysis_name}__analysis_samplefile.txt"
    log:
        f"logs/combined_analysis__{analysis_name}.log"
    shell:
        """
        touch {output.chkpt}
        """ 
