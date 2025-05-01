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
                      names=["data_type", "line", "tissue", "sample", "replicate", 
                             "seq_id", "fastq_path", "paired", "ref_genome"])

# Create a dictionary to store the information for each sample
sample_info_map = {
    (row.line, row.tissue, row.sample, row.replicate): {
        "seq_id": row.seq_id,
        "fastq_path": row.fastq_path,
        "paired": row.paired,
        "ref_genome": row.ref_genome,
        "data_type": row.data_type,
    }
    for _, row in samples.iterrows()
}

print(sample_info_map, file=sys.stderr)

# Function to access this information later on
def get_sample_info(wildcards, field):
    key = (wildcards.line, wildcards.tissue, wildcards.sample, wildcards.replicate)
    return sample_info_map[key][field]

# Generate all sample output files required
all_sample_outputs = expand(
    "chkpts/process__{data_type}__{line}__{tissue}__{sample}__{replicate}__{ref_genome}.done",
    zip,
    data_type = samples["data_type"],
    line = samples["line"],
    tissue = samples["tissue"],
    sample = samples["sample"],
    replicate = samples["replicate"],
    ref_genome = samples["ref_genome"]
)

# Define reference genomes
REF_GENOMES = set(samples["ref_genome"].unique())

# Define data types
DATA_TYPES = set(samples["data_type"].unique())

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

# Rule all to specify final target
rule all:
	input:
		f"chkpts/combined_analysis__{analysis_name}.done"

# Rule to prepare reference genome for each data type
rule prepare_reference:
    input:
        refs = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
    output:
        chkpt = "chkpts/ref__{ref_genome}__{env}.done"
    params:
        ref_path = config["ref_path"],
        scripts_dir = config["scripts_dir"]
    log:
        "logs/prepare_ref__{ref_genome}__{env}.log"
    conda:
        "envs/reference.yaml"
    shell:
        """
        # Call the original environment preparation script
        qsub {params.scripts_dir}/MaizeCode_check_environment.sh \
            -p {params.ref_path} \
            -r {wildcards.ref_genome} \
            -d {wildcards.env} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to process samples based on data type
rule process_sample:
    input:
        ref_chkpt = lambda wildcards: f"chkpts/ref__{get_sample_info(wildcards, 'ref_genome')}__{datatype_to_env[get_sample_info(wildcards, 'data_type')]}.done"
    output:
        chkpt = "chkpts/process__{data_type}__{line}__{tissue}__{sample}__{replicate}__{ref_genome}.done"
    params:
        scripts_dir = config["scripts_dir"],
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], get_sample_info(wildcards, 'ref_genome')),
        env = lambda wildcards: datatype_to_env[get_sample_info(wildcards, 'data_type')],
        line = lambda wildcards: wildcards.line,
        tissue = lambda wildcards: wildcards.tissue,
        replicate = lambda wildcards: wildcards.replicate,
        seq_id = lambda wildcards: get_sample_info(wildcards, 'seq_id'),
        fastq_path = lambda wildcards: get_sample_info(wildcards, 'fastq_path'),
        paired = lambda wildcards: get_sample_info(wildcards, 'paired'),
        mapping_option = config["mapping_option"]
    log:
        "logs/process__{data_type}__{line}__{tissue}__{sample}__{replicate}__{ref_genome}.log"
    conda:
        lambda wildcards: f"envs/{datatype_to_env[get_sample_info(wildcards, 'data_type')]}_sample.yaml"
    shell:
        """
        qsub {params.scripts_dir}/MaizeCode_{params.env}_sample.sh \
            -x {wildcards.sample} \
            -d {params.ref_dir} \
            -l {params.line} \
            -t {params.tissue} \
            -m {wildcards.data_type} \
            -r {params.replicate} \
            -i {params.seq_id} \
            -f {params.fastq_path} \
            -p {params.paired} \
            -s "download" \
            -a {params.mapping_option} > {log} 2>&1
        touch {output.chkpt}
        """

# # Rule to perform data type specific analysis
# rule analyze_sample:
    # input:
        # process_chkpt = lambda wildcards: [
            # f"chkpts/sample__{wildcards.data_type}__{sample}__{rep}__{samples[samples['sample'] == sample]['ref_genome'].iloc[0]}.done"
            # for sample in datatype_to_samples[wildcards.data_type]
            # for rep in samples_to_replicates[sample]
        # ]
    # output:
        # chkpt = "chkpts/analysis__{data_type}__{analysis_name}.done"
    # params:
        # scripts_dir = config["scripts_dir"],
        # analysis_samplefile = f"{analysis_name}__analysis_samplefile.txt",
        # env = lambda wildcards: datatype_to_env[wildcards.data_type]
    # log:
        # "logs/analysis__{data_type}__{analysis_name}.log"
    # conda:
        # lambda wildcards: f"envs/{datatype_to_env[wildcards.data_type]}_sample.yaml"
    # shell:
        # """
        # # Call the appropriate analysis script based on data type
        # qsub {params.scripts_dir}/MaizeCode_{params.env}_analysis.sh \
            # -f {params.analysis_samplefile} > {log} 2>&1
        # touch {output.chkpt}
        # """

# Rule to perform combined analysis
rule combined_analysis:
    input:
        sample_chkpt = all_sample_outputs
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
        
