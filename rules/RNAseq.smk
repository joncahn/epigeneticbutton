# function to access logs more easily
def return_log_rna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"RNA","logs",f"tmp_rna__{sample_name}__{step}__{paired}.log")    

def get_inputs_rna(wildcards):
    s = {k: getattr(wildcards, k) for k in ["data_type","line", "tissue", "sample_type", "replicate", "ref_genome"]}
    name = sample_name(s)
    paired = get_sample_info(wildcards, "paired")
    if paired == "PE":
        return f"RNA/chkpts/temp_pe__{name}.done"
    else:
        return f"RNA/chkpts/temp_se__{name}.done"
        
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/RNA_sample.yaml")

rule stat_file_rna:
    output:
        stat_file = "RNA/reports/summary_mapping_stats.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """

rule make_RNA_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa",
        gtf = "genomes/{ref_genome}/temp_{ref_genome}.gtf"
    output:
        indices = "genomes/{ref_genome}/STAR_index"
    log:
        os.path.join(REPO_FOLDER,"logs","STAR_index_{ref_genome}.log")
    threads: workflow.cores
    shell:
        """
        if [ ! -d {output.indices} ]; then
            printf "\nBuilding STAR index directory for {ref_genome}\n"
            mkdir {output.indices}
            STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.indices} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf}
        else
            printf "\nSTAR index already exists for {ref_genome}\n"
        fi
        """

rule STAR_map_pe:
    input:
        fastq1 = "RNA/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "RNA/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"combined/genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        prefix = "RNA/mapped/map_pe__{sample_name}_",
        touch = "RNA/chkpts/temp_pe__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        return_log_rna("{sample_name}", "mapping", "PE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nMapping {sample_name} to {ref_genome} with STAR version:\n"
        STAR --version
        STAR --runMode alignReads --genomeDir {input.indices} --readFilesIn {input.fastq1} {input.fastq2} --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix {output.prefix} --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 20 --quantMode GeneCounts
        touch {output.touch}
        """    

rule STAR_map_se:
    input:
        fastq0 = "RNA/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"combined/genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        prefix = "RNA/mapped/map_se__{sample_name}_",
        touch = "RNA/chkpts/temp_se__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']    
    log:
        return_log_rna("{sample_name}", "mapping", "SE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nMapping {sample_name} to {ref_genome} with STAR version:\n"
        STAR --version
        STAR --runMode alignReads --genomeDir {input.indices} --readFilesIn {input.fastq0} --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix {output.prefix} --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --quantMode GeneCounts
        touch {output.touch}
        """
        
rule check_pair_rna:
    input:
        lambda wildcards: get_inputs_rna(wildcards)
    output:
        touch = "RNA/chkpts/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    shell:
        """
        touch {output.touch}
        """
