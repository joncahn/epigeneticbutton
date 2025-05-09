# function to access logs more easily
def return_log_sample(data_type, sample_name, step, paired):
    return os.path.join(REPO_FOLDER,f"{data_type}","logs",f"tmp__{sample_name}__{step}__{paired}.log")

CONDA_ENV=os.path.join(REPO_FOLDER,"envs/download.yaml")
    
rule get_fastq_pe:
    output:
        fastq1 = "{data_type}/fastq/raw__{sample_name}__R1.fastq.gz",
        fastq2 = "{data_type}/fastq/raw__{sample_name}__R2.fastq.gz"
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type
    log:
        return_log_sample("{data_type}","{sample_name}", "downloading", "PE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {
        if [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "Using fasterq-dump for {params.sample_name} ({params.seq_id})\n"
            fasterq-dump -e {threads} --outdir "{params.data_type}/fastq" "{params.seq_id}"
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files\n"
            pigz -p {threads} "{params.data_type}/fastq/{params.seq_id}_1.fastq"
            mv "{params.data_type}/fastq/{params.seq_id}_1.fastq.gz" "{output.fastq1}"
            pigz -p {threads} "{params.data_type}/fastq/{params.seq_id}_2.fastq"
            mv "{params.data_type}/fastq/{params.seq_id}_2.fastq.gz" "{output.fastq2}"
        else
            printf "Copying PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R1*q.gz "{output.fastq1}"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R2*q.gz "{output.fastq2}"
        fi
        } 2>&1 | tee -a "{log}"
        """

        
rule get_fastq_se:
    output:
        fastq0 = "{data_type}/fastq/raw__{sample_name}__R0.fastq.gz"
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type
    log:
        return_log_sample("{data_type}","{sample_name}", "downloading", "SE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {
        if [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "Using fasterq-dump for {params.sample_name} ({params.seq_id})\n"
            fasterq-dump -e {threads} --outdir "{params.data_type}/fastq" "{params.seq_id}"
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files\n"
            pigz -p {threads} "{params.data_type}/fastq/{params.seq_id}.fastq"
            mv "{params.data_type}/fastq/{params.seq_id}.fastq.gz" "{output.fastq0}"
        else
            printf "\nCopying SE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/"{params.seq_id}"*q.gz "{output.fastq0}"
        fi
        } 2>&1 | tee -a "{log}"        
        """

rule process_fastq_pe:
    input:
        raw_fastq1 = "{data_type}/fastq/raw__{sample_name}__R1.fastq.gz",
        raw_fastq2 = "{data_type}/fastq/raw__{sample_name}__R2.fastq.gz"
    output:
        fastq1 = "{data_type}/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "{data_type}/fastq/trim__{sample_name}__R2.fastq.gz",
        metrics = "{data_type}/reports/trim_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        adapter1 = "AGATCGGAAGAGCACACGTCTGAAC",
        adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGA",
        trimming_quality = config['trimming_quality']
    log:
        return_log_sample("{data_type}","{sample_name}", "trimming", "PE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {
        #### printf "\nRunning fastQC for {params.sample_name} with fastqc version:\n"
        fastqc --version
        fastqc -o "{params.data_type}/reports/" "{input.raw_fastq1}"
        fastqc -o "{params.data_type}/reports/" {input.raw_fastq2}"
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a "{params.adapter1}" -A "{params.adapter2}" -o "{output.fastq1}" -p "{output.fastq2}" "{input.raw_fastq1}" "{input.raw_fastq2}" 2>&1 | tee "{output.metrics}"
		#### Removing untrimmed fastq
		rm -f "{input.raw_fastq1}" "{input.raw_fastq2}"
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for {params.sample_name}\n"
		fastqc -o "{params.data_type}/reports/" "{output.fastq1}"
		fastqc -o "{params.data_type}/reports/" "{output.fastq2}"
        } 2>&1 | tee -a "{log}"        
        """
        
rule process_fastq_se:
    input:
        raw_fastq = "{data_type}/fastq/raw__{sample_name}__R0.fastq.gz"
    output:
        fastq = "{data_type}/fastq/trim__{sample_name}__R0.fastq.gz",
        metrics = "{data_type}/reports/trim_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        adapter1 = "AGATCGGAAGAGCACACGTCTGAAC",
        trimming_quality = config['trimming_quality']
    log:
        return_log_sample("{data_type}","{sample_name}", "trimming", "SE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {
        ### QC of the raw reads with FastQC
        printf "\nRunning fastQC for {params.sample_name} with fastqc version:\n"
        fastqc --version
        fastqc -o "{params.data_type}/reports/" "{input.raw_fastq}"
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a "{params.adapter1}" -o "{output.fastq}" "{input.raw_fastq}" 2>&1 | tee "{output.metrics}"
		#### Removing untrimmed fastq
		rm -f "{input.raw_fastq}"
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for {params.sample_name}\n"
		fastqc -o "{params.data_type}/reports/" "{output.fastq}"
        } 2>&1 | tee -a "{log}"
        """
