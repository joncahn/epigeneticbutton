# function to access logs more easily
def return_log_sample(data_type, sample_name, step, paired):
    return os.path.join(REPO_FOLDER,f"{data_type}","logs",f"tmp__{sample_name}__{step}__{paired}.log")
    
rule get_fastq_pe:
    input:
        setup = "chkpts/directories_setup.done"
    output:
        fastq1 = temp("{data_type}/fastq/raw__{sample_name}__R1.fastq.gz"),
        fastq2 = temp("{data_type}/fastq/raw__{sample_name}__R2.fastq.gz")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "downloading", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["fastq_dump"]["threads"]
    resources:
        mem=config["resources"]["fastq_dump"]["mem"],
        tmp=config["resources"]["fastq_dump"]["tmp"]
    shell:
        """
        {{
        if [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "Using fasterq-dump for {params.sample_name} ({params.seq_id})\n"
            fasterq-dump -e {threads} --outdir "{params.data_type}/fastq" "{params.seq_id}"
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files\n"
            pigz -p {threads} "{params.data_type}/fastq/{params.seq_id}_1.fastq"
            mv "{params.data_type}/fastq/{params.seq_id}_1.fastq.gz" "{output.fastq1}"
            pigz -p {threads} "{params.data_type}/fastq/{params.seq_id}_2.fastq"
            mv "{params.data_type}/fastq/{params.seq_id}_2.fastq.gz" "{output.fastq2}"
        elif ls "{params.fastq_path}"/*"{params.seq_id}"*R1*q.gz 1> /dev/null 2>&1 && ls "{params.fastq_path}"/*"{params.seq_id}"*R2*q.gz 1> /dev/null 2>&1; then
            printf "Copying PE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R1*q.gz "{output.fastq1}"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R2*q.gz "{output.fastq2}"
        elif ls "{params.fastq_path}"/*"{params.seq_id}"*R1*q 1> /dev/null 2>&1 && ls "{params.fastq_path}"/*"{params.seq_id}"*R2*q 1> /dev/null 2>&1; then
            printf "Copying and gzipping PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*R1*q -c > "{output.fastq1}"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*R2*q -c > "{output.fastq2}"
        else
            printf "Error: No PE fastqs found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule get_fastq_se:
    input:
        setup = "chkpts/directories_setup.done"
    output:
        fastq0 = temp("{data_type}/fastq/raw__{sample_name}__R0.fastq.gz")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "downloading", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["fastq_dump"]["threads"]
    resources:
        mem=config["resources"]["fastq_dump"]["mem"],
        tmp=config["resources"]["fastq_dump"]["tmp"]
    shell:
        """
        {{
        if [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "Using fasterq-dump for {params.sample_name} ({params.seq_id})\n"
            fasterq-dump -e {threads} --outdir "{params.data_type}/fastq" "{params.seq_id}"
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files\n"
            pigz -p {threads} "{params.data_type}/fastq/{params.seq_id}.fastq"
            mv "{params.data_type}/fastq/{params.seq_id}.fastq.gz" "{output.fastq0}"
        elif ls "{params.fastq_path}"/*"{params.seq_id}"*q.gz 1> /dev/null 2>&1; then
            printf "\nCopying SE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*q.gz "{output.fastq0}"
        elif ls "{params.fastq_path}"/*"{params.seq_id}"*q 1> /dev/null 2>&1; then
            printf "\nCopying and gzipping SE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*q -c > "{output.fastq0}"          
        else
            printf "Error: No SE fastq found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
        fi
        }} 2>&1 | tee -a "{log}"        
        """

rule run_fastqc:
    input:
        fastq = "{data_type}/fastq/{step}__{sample_name}__{read}.fastq.gz"
    output:
        fastqc = "{data_type}/reports/{step}__{sample_name}__{read}_fastqc.html"
    params:
        data_type = lambda wildcards: wildcards.data_type,
        step = lambda wildcards: wildcards.step,
        sample_name = lambda wildcards: wildcards.sample_name,
        read = lambda wildcards: wildcards.read
    conda: CONDA_ENV
    threads: 1
    resources:
        mem=config["resources"]["fastqc"]["mem"],
        tmp=config["resources"]["fastqc"]["tmp"]
    shell:
        """
        fastqc -o "{params.data_type}/reports/" "{input.fastq}"
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
        trimming_quality = lambda wildcards: config['trimming_quality'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')]
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "trimming", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["process_fastq"]["threads"]
    resources:
        mem=config["resources"]["process_fastq"]["mem"],
        tmp=config["resources"]["process_fastq"]["tmp"]
    shell:
        """
        {{
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a "{params.adapter1}" -A "{params.adapter2}" -o "{output.fastq1}" -p "{output.fastq2}" "{input.raw_fastq1}" "{input.raw_fastq2}" 2>&1 | tee "{output.metrics}"
        }} 2>&1 | tee -a "{log}"        
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
        adapter1 = lambda wildcards: config['adapter1'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimming_quality = lambda wildcards: config['trimming_quality'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')]
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "trimming", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["process_fastq"]["threads"]
    resources:
        mem=config["resources"]["process_fastq"]["mem"],
        tmp=config["resources"]["process_fastq"]["tmp"]
    shell:
        """
        {{
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a "{params.adapter1}" -o "{output.fastq}" "{input.raw_fastq}" 2>&1 | tee "{output.metrics}"
        }} 2>&1 | tee -a "{log}"
        """
