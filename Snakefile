import pandas as pd
import os
import fnmatch
import sys
from snakemake.utils import min_version

min_version("6.0")

# Configuration
configfile: "config.yaml"

# Define sample metadata
samples = pd.read_csv(config["sample_file"], sep="\t", header=None,
                      names=["data_type", "line", "tissue", "sample_type", "replicate", 
                             "seq_id", "fastq_path", "paired", "ref_genome"])

# Create a dictionary to store the information for each sample
sample_info_map = {
    (row["data_type"], row["line"], row["tissue"], row["sample_type"], row["replicate"], row["ref_genome"]): {
        "seq_id": row["seq_id"],
        "fastq_path": row["fastq_path"],
        "paired": row["paired"]
    }
    for _, row in samples.iterrows()
}

# Function to create a unique name for each sample based on the sample columns, and later based on wildcards
def sample_name(d):
    return f"{d['data_type']}__{d['line']}__{d['tissue']}__{d['sample_type']}__{d['replicate']}__{d['ref_genome']}"

# Function to split the sample_name to recover its components
def parse_sample_name(sample_name):
    data_type, line, tissue, sample_type, rep, ref_genome = sample_name.split("__")
    return {
        "data_type": data_type,
        "line": line,
        "tissue": tissue,
        "sample_type": sample_type,
        "replicate": rep,
        "ref_genome": ref_genome
    }

# Function to access extra information form the samplefile using wildcards
def get_sample_info(wildcards, field):
    key = (wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.replicate, wildcards.ref_genome)
    return sample_info_map[key][field]

# Function to access extra information form the samplefile using the name
def get_sample_info_from_name(sample_name, field):
    parts = sample_name.split("__")
    key = tuple(parts)
    return sample_info_map[key][field]

# Generate all sample output files required
all_sample_outputs = expand(
    "chkpts/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done",
    zip,
    data_type = samples["data_type"],
    line = samples["line"],
    tissue = samples["tissue"],
    sample_type = samples["sample_type"],
    replicate = samples["replicate"],
    ref_genome = samples["ref_genome"]
)

# Define reference genomes and the path to them
REF_GENOMES = set(samples["ref_genome"].unique())
REF_PATH = config["ref_path"]

# Define data types
DATA_TYPES = set(samples["data_type"].unique())

# Define the folde rin which the snakemake pipeline has been cloned
REPO_FOLDER = config["repo_folder"]

# Define label for the analysis
analysis_name = config["analysis_name"]

# Define patterns and envs
env_patterns = [
    ("ChIP*", "ChIP"),
    ("RNAseq", "RNA"),
    ("RAMPAGE", "RNA"),
    ("shRNA", "shRNA"),
    ("mC", "mC"),
    ("TF_*", "TF"),
]

# Function to determine environment
def get_env(data_type):
    for pattern, env in env_patterns:
        if fnmatch.fnmatch(data_type, pattern):
            return env
    return "unknown"

# Map data types to environments
datatype_to_env = {dt: get_env(dt) for dt in DATA_TYPES}
UNIQUE_ENVS = list(set(datatype_to_env.values()))

# Check for unknown envs and exit if any
unknowns = [dt for dt, env in datatype_to_env.items() if env == "unknown"]
if unknowns:
    print("Type of data unknown for the following data types:")
    for dt in unknowns:
        print(f"  - {dt}")
    print("\nPlease check your sample sheet or update the env_patterns.")
    sys.exit(1)

# Create dictionaries
refgenome_to_datatype = samples.groupby("ref_genome")["data_type"].unique().to_dict()
refgenome_to_env = {}
for ref, dtypes in refgenome_to_datatype.items():
    envs = {datatype_to_env.get(dt) for dt in dtypes if dt in datatype_to_env}
    if envs:
        refgenome_to_env[ref] = list(envs)

# Load the sample metadata and perform all operations in a single chain
analysis_samples = (
    samples
    .query("sample_type != 'Input'") # filter Input samples
    .assign(ref_dir=lambda df: df.apply(lambda row: os.path.join(REF_PATH, row["ref_genome"]), axis=1)) # create a column with the path to reference genome
    [["data_type", "line", "tissue", "sample_type", "paired", "ref_dir"]] # select the necessary columns
    .drop_duplicates() # removes replicates
)

# Save the result to 'analysis_samplefile.txt'
analysis_samples.to_csv(f"{analysis_name}__analysis_samplefile.txt", sep="\t", index=False)

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
    
    for key, value in dirs.items():
        if isinstance(value, dict):
            for sub_value in value.values():
                os.makedirs(sub_value, exist_ok=True)
        else:
            os.makedirs(value, exist_ok=True)

# Call the function to create directories
create_directories(UNIQUE_ENVS, DIRS)

# Include all rule files
include: "rules/environment_setup.smk"
include: "rules/ChIPseq.smk"

# Rule all to specify final target
rule all:
	input:
		f"chkpts/combined_analysis__{analysis_name}.done"

# # Rule to process samples based on data type
# rule process_rna_sample:
    # input:
        # ref_chkpt = lambda wildcards: f"chkpts/ref__{get_sample_info(wildcards, 'ref_genome')}.done"
    # output:
        # chkpt = "chkpts/process__RNA__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    # params:
        # scripts_dir = os.path.join(REPO_FOLDER,"scripts"),
        # ref_dir = lambda wildcards: os.path.join(REF_PATH, get_sample_info(wildcards, 'ref_genome')),
        # env = lambda wildcards: datatype_to_env[wildcards.data_type],
        # line = lambda wildcards: wildcards.line,
        # tissue = lambda wildcards: wildcards.tissue,
        # replicate = lambda wildcards: wildcards.replicate,
        # seq_id = lambda wildcards: get_sample_info(wildcards, 'seq_id'),
        # fastq_path = lambda wildcards: get_sample_info(wildcards, 'fastq_path'),
        # paired = lambda wildcards: get_sample_info(wildcards, 'paired'),
        # mapping_option = config["mapping_option"]
    # log:
        # "logs/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    # conda:
        # lambda wildcards: f"envs/{datatype_to_env[wildcards.data_type]}_sample.yaml"
    # shell:
        # """
        # cd {params.env}/
        # qsub ../{params.scripts_dir}/MaizeCode_{params.env}_sample.sh \
            # -x {wildcards.sample_type} \
            # -d {params.ref_dir} \
            # -l {params.line} \
            # -t {params.tissue} \
            # -m {wildcards.data_type} \
            # -r {params.replicate} \
            # -i {params.seq_id} \
            # -f {params.fastq_path} \
            # -p {params.paired} \
            # -s "download" \
            # -a {params.mapping_option} | tee {log}
        # cd ..
        # touch {output.chkpt}
        # """

# # Rule to prepare the file containing the path to regions bed files to use for analysis
# rule prepare_region_file:
    # input:
        # [
            # f"{env}/tracks/{ref}_all_genes.bed"
            # for ref, envs in refgenome_to_env.items()
            # for env in envs
        # ]
    # output:
        # region_file="all_genes.txt"
    # run:
        # with open(output.region_file, "w") as outfile:
            # for path in input:
                # if os.path.isfile(path) and os.path.getsize(path) > 0:
                    # outfile.write(f"{path}\n")        

# Rule to perform combined analysis
rule combined_analysis:
    input:
        expand("ChIP/chkpts/process__{sample_name}.done", sample_name=samples[samples["data_type"] == "ChIP"].apply(sample_name, axis=1)),
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
        # Call the combined analysis script
        qsub {params.scripts_dir}/MaizeCode_analysis.sh \
            -f {params.analysis_samplefile} \
            -r {params.region_file} | tee {log}
        touch {output.chkpt}
        """ 
        
