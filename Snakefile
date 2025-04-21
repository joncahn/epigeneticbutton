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

print(DATA_TYPES)

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

# Create output directories
for data_type in DATA_TYPES:
    DIRS[data_type] = {
        "fastq": f"{data_type}/fastq",
        "mapped": f"{data_type}/mapped",
        "tracks": f"{data_type}/tracks",
        "reports": f"{data_type}/reports",
        "logs": f"{data_type}/logs",
        "chkpts": f"{data_type}/chkpts",
        "plots": f"{data_type}/plots"
    }

print(DIRS)

# Rule to create data type-specific directories
rule create_data_type_directories:
	output:
		expand("{data_type}/fastq", data_type=DATA_TYPES),
		expand("{data_type}/mapped", data_type=DATA_TYPES),
		expand("{data_type}/tracks", data_type=DATA_TYPES),
		expand("{data_type}/reports", data_type=DATA_TYPES),
		expand("{data_type}/logs", data_type=DATA_TYPES),
		expand("{data_type}/chkpts", data_type=DATA_TYPES),
		expand("{data_type}/plots", data_type=DATA_TYPES)
	shell:
		"""
		for data_type in {data_types}; do
			mkdir -p ${data_type}/fastq
			mkdir -p ${data_type}/mapped
			mkdir -p ${data_type}/tracks
			mkdir -p ${data_type}/reports
			mkdir -p ${data_type}/logs
			mkdir -p ${data_type}/chkpts
			mkdir -p ${data_type}/plots
		done
		""".format(data_types=" ".join(DATA_TYPES))

# Create all directories
rule create_directories:
    output:
        directory("chkpts"),
        directory("combined/peaks"),
        directory("combined/DEG"),
        directory("combined/TSS"),
        directory("combined/reports"),
        directory("combined/matrix"),
        directory("combined/plots"),
        directory("combined/chkpts"),
        directory("combined/logs")
    shell:
        """
        mkdir -p chkpts
        mkdir -p combined/peaks
        mkdir -p combined/DEG
        mkdir -p combined/TSS
        mkdir -p combined/reports
        mkdir -p combined/matrix
        mkdir -p combined/plots
        mkdir -p combined/chkpts
        mkdir -p combined/logs
        """

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
        ref_chkpt = "chkpts/ref_{ref_genome}_{data_type}.done"
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

# Rule to perform ChIP specific analysis
rule analyze_ChIP:
    input:
        ref_chkpt = "chkpts/sample_{data_type}_{sample}_{replicate}.done"
    output:
        chkpt = "chkpts/analysis_ChIP_{config.analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/analysis_ChIP_{config.analysis_name}.log"
    conda:
        "envs/{data_type}_analysis.yaml"
    shell:
        """
        # Call the appropriate analysis script based on data type
        {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            -f {input.sample_file} \
            -p {config["ref_path"]} \
            -r {wildcards.ref_genome} \
            -d {wildcards.data_type} > {log} 2>&1
        touch {output.chkpt}
        """
		
# Rule to perform RNA specific analysis
rule analyze_RNA:
    input:
        ref_chkpt = "chkpts/sample_{data_type}_{sample}_{replicate}.done"
    output:
        chkpt = "chkpts/analysis_RNA_{config.analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/analysis_RNA_{config.analysis_name}.log"
    conda:
        "envs/{data_type}_analysis.yaml"
    shell:
        """
        # Call the appropriate analysis script based on data type
        {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            -f {input.sample_file} \
            -p {config["ref_path"]} \
            -r {wildcards.ref_genome} \
            -d {wildcards.data_type} > {log} 2>&1
        touch {output.chkpt}
        """
		
# Rule to perform TF specific analysis
rule analyze_TF:
    input:
        ref_chkpt = "chkpts/sample_{data_type}_{sample}_{replicate}.done"
    output:
        chkpt = "chkpts/analysis_TF_{config.analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/analysis_TF_{config.analysis_name}.log"
    conda:
        "envs/{data_type}_analysis.yaml"
    shell:
        """
        # Call the appropriate analysis script based on data type
        {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            -f {input.sample_file} \
            -p {config["ref_path"]} \
            -r {wildcards.ref_genome} \
            -d {wildcards.data_type} > {log} 2>&1
        touch {output.chkpt}
        """

# Rule to perform mC specific analysis
rule analyze_mC:
    input:
        ref_chkpt = "chkpts/sample_{data_type}_{sample}_{replicate}.done"
    output:
        chkpt = "chkpts/analysis_mC_{config.analysis_name}.done"
    params:
        scripts_dir = config["scripts_dir"]
    log:
        "logs/analysis_mC_{config.analysis_name}.log"
    conda:
        "envs/{data_type}_analysis.yaml"
    shell:
        """
        # Call the appropriate analysis script based on data type
        {params.scripts_dir}/MaizeCode_{wildcards.data_type}_analysis.sh \
            -f {input.sample_file} \
            -p {config["ref_path"]} \
            -r {wildcards.ref_genome} \
            -d {wildcards.data_type} > {log} 2>&1
        touch {output.chkpt}
        """

		
# Rule to perform combined analysis
rule combined_analysis:
    input:
        ref_chkpt = "chkpts/analysis_{data_type}_{config.analysis_name}.done"
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
