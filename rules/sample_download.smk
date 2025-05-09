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
    run:
        from pathlib import Path
        import shutil, subprocess

        lf = open(log[0], "a")

        if params.fastq_path == "SRA":
            lf.write(f"\nUsing fasterq-dump for {params.sample_name} ({params.seq_id})\n")
            subprocess.run(["fasterq-dump", "-e", str(threads), "--outdir", f"{params.data_type}/fastq", params.seq_id], check=True, stdout=lf, stderr=lf)
            for r in [1, 2]:
                subprocess.run(["pigz", "-p", str(threads), f"{params.data_type}/fastq/{params.seq_id}_{r}.fastq"],check=True, stdout=lf, stderr=lf)
                shutil.move(f"{params.data_type}/fastq/{params.seq_id}_{r}.fastq.gz",output[f"fastq{r}"])
        else:
            lf.write(f"\nCopying PE fastq for {params.sample_name} from {params.fastq_path}\n")
            for r in [1, 2]:
                matches = list(Path(params.fastq_path).glob(f"*{params.seq_id}*R{r}*q.gz"))
                if not matches:
                    raise FileNotFoundError(f"Missing R{r} FASTQ for {params.seq_id}")
                shutil.copy(matches[0], output[f"fastq{r}"])
        lf.close()

        
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
        if [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "\nUsing fasterq-dump for {params.sample_name} ({params.seq_id})\n" >> {log} 2>&1
            fasterq-dump -e {threads} --outdir {params.data_type}/fastq {params.seq_id}
            printf "\n {params.sample_name} ({params.seq_id}) downloaded\nRenaming files..." >> {log} 2>&1
            pigz -p {threads} {params.data_type}/fastq/{params.seq_id}.fastq
            mv {params.data_type}/fastq/{params.seq_id}.fastq.gz {output.fastq0}
        else
            printf "\nCopying SE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n" >> {log} 2>&1
            cp {params.fastq_path}/${params.seq_id}*q.gz {output.fastq0}
        fi
        """

rule process_fastq_pe:
    input:
        raw_fastq1 = "{data_type}/fastq/raw__{sample_name}__R1.fastq.gz",
        raw_fastq2 = "{data_type}/fastq/raw__{sample_name}__R2.fastq.gz"
    output:
        fastq1 = "{data_type}/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "{data_type}/fastq/trim__{sample_name}__R2.fastq.gz",
        metrics_trim = "{data_type}/reports/trim_pe__{sample_name}.txt"
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
        #### printf "\nRunning fastQC for {params.sample_name} with fastqc version:\n"
        fastqc --version
        fastqc -o {params.data_type}/reports/ {input.raw_fastq1}
        fastqc -o {params.data_type}/reports/ {input.raw_fastq2}	
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a {params.adapter1} -A {params.adapter2} -o {output.fastq1} -p {output.fastq2} {input.raw_fastq1} {input.raw_fastq2} |& tee {output.metrics}
		#### Removing untrimmed fastq
		rm -f {input.raw_fastq1} {input.raw_fastq2}
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for {params.sample_name}\n"
		fastqc -o {params.data_type}/reports/ {output.fastq1}
		fastqc -o {params.data_type}/reports/ {output.fastq2}
        """
        
rule process_fastq_se:
    input:
        raw_fastq = "{data_type}/fastq/raw__{sample_name}__R0.fastq.gz"
    output:
        fastq = "{data_type}/fastq/trim__{sample_name}__R0.fastq.gz",
        metrics_trim = "{data_type}/reports/trim_se__{sample_name}.txt"
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
        ### QC of the raw reads with FastQC
        printf "\nRunning fastQC for {params.sample_name} with fastqc version:\n"
        fastqc --version
        fastqc -o {params.data_type}/reports/ {input.raw_fastq}
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a {params.adapter1} -o {output.fastq} {input.raw_fastq} |& tee {output.metrics}
		#### Removing untrimmed fastq
		rm -f {input.raw_fastq}
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for {params.sample_name}\n"
		fastqc -o {params.data_type}/reports/ {output.fastq}
        """
