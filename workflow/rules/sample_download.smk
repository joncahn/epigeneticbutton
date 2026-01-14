def return_log_sample(data_type, sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"results",data_type,"logs",f"tmp__{sample_name}__{step}__{paired}.log")
    
rule get_fastq_pe:
    output:
        fastq1 = temp("results/{data_type}/fastq/raw__{sample_name}__R1.fastq.gz"),
        fastq2 = temp("results/{data_type}/fastq/raw__{sample_name}__R2.fastq.gz")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        trimmed_fastqs = config['trimmed_fastqs'],
        exist_fastq1 = lambda wildcards: f"results/{wildcards.data_type}/fastq/trim__{wildcards.sample_name}__R1.fastq.gz",
        exist_fastq2 = lambda wildcards: f"results/{wildcards.data_type}/fastq/trim__{wildcards.sample_name}__R2.fastq.gz"
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "downloading", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["get_fastq_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["get_fastq_pe"]["mem_mb"],
        tmp_mb=config["resources"]["get_fastq_pe"]["tmp_mb"],
        qos=config["resources"]["get_fastq_pe"]["qos"]
    retries: 3
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" && -e "{params.exist_fastq1}" && -e "{params.exist_fastq2}" ]]; then
            printf "Fastqs already exist for PE {params.sample_name}\n"
            cp {params.exist_fastq1} {output.fastq1}
            cp {params.exist_fastq2} {output.fastq2}
        elif [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "Using fasterq-dump for PE {params.sample_name} ({params.seq_id})\n"
            numbers=$(echo "{params.seq_id}" | sed 's/,/ /g')
            fastq_files_r1=()
            fastq_files_r2=()
            for nb in ${{numbers}}; do
                fasterq-dump -e {threads} --outdir "results/{params.data_type}/fastq" "${{nb}}"
                fastq_files_r1+=("results/{params.data_type}/fastq/${{nb}}_1.fastq")
                fastq_files_r2+=("results/{params.data_type}/fastq/${{nb}}_2.fastq")
            done
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files\n"
            cat "${{fastq_files_r1[@]}}" > "results/{params.data_type}/fastq/raw__{params.sample_name}__R1.fastq"
            pigz -p {threads} "results/{params.data_type}/fastq/raw__{params.sample_name}__R1.fastq"
            rm -f "${{fastq_files_r1[@]}}"
            cat "${{fastq_files_r2[@]}}" > "results/{params.data_type}/fastq/raw__{params.sample_name}__R2.fastq"
            pigz -p {threads} "results/{params.data_type}/fastq/raw__{params.sample_name}__R2.fastq"
            rm -f "${{fastq_files_r2[@]}}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q.gz 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q.gz 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying PE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q.gz "{output.fastq1}"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q.gz "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying and gzipping PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q -c > "{output.fastq1}"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q -c > "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q.gz 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q.gz 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying PE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q.gz "{output.fastq1}"
            cp "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q.gz "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying and gzipping PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q -c > "{output.fastq1}"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q -c > "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*1*f*q 2>/dev/null | wc -l) -gt 1 ]]; then
            printf "Error: Too many fastqs found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\nThe seq_id used {params.seq_id} is likely not unique.\n"
        else
            printf "Error: No PE fastqs found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule get_fastq_se:
    output:
        fastq0 = temp("results/{data_type}/fastq/raw__{sample_name}__R0.fastq.gz")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        trimmed_fastqs = config['trimmed_fastqs'],
        exist_fastq0 = lambda wildcards: f"results/{wildcards.data_type}/fastq/raw__{wildcards.sample_name}__R0.fastq.gz"
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "downloading", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["get_fastq_se"]["threads"]
    resources:
        mem_mb=config["resources"]["get_fastq_se"]["mem_mb"],
        tmp_mb=config["resources"]["get_fastq_se"]["tmp_mb"],
        qos=config["resources"]["get_fastq_se"]["qos"]
    retries: 3
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" && -e "{params.exist_fastq0}" ]]; then
            printf "Fastq already existing for SE {params.sample_name}\n"
            cp {params.exist_fastq0} {output.fastq0}
        elif [[ "{params.fastq_path}" == "SRA" ]]; then
            printf "Using fasterq-dump for SE {params.sample_name} ({params.seq_id})\n"
            numbers=$(echo "{params.seq_id}" | sed 's/,/ /g')
            fastq_files=()
            for nb in ${{numbers}}; do
                fasterq-dump -e {threads} --outdir "results/{params.data_type}/fastq" "${{nb}}"
                fastq_files+=("results/{params.data_type}/fastq/${{nb}}.fastq")
            done
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files\n"
            cat "${{fastq_files[@]}}" > "results/{params.data_type}/fastq/raw__{params.sample_name}__R0.fastq"
            pigz -p {threads} "results/{params.data_type}/fastq/raw__{params.sample_name}__R0.fastq"
            rm -f "${{fastq_files[@]}}"
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
        fastq = "results/{data_type}/fastq/{step}__{sample_name}__{read}.fastq.gz"
    output:
        fastqc = "results/{data_type}/reports/{step}__{sample_name}__{read}_fastqc.html"
    params:
        data_type = lambda wildcards: wildcards.data_type,
        step = lambda wildcards: wildcards.step,
        sample_name = lambda wildcards: wildcards.sample_name,
        read = lambda wildcards: wildcards.read
    conda: CONDA_ENV
    threads: 1
    resources:
        mem_mb=config["resources"]["run_fastqc"]["mem_mb"],
        tmp_mb=config["resources"]["run_fastqc"]["tmp_mb"],
        qos=config["resources"]["run_fastqc"]["qos"]
    shell:
        """
        fastqc -o "results/{params.data_type}/reports/" "{input.fastq}"
        """

rule process_fastq_pe:
    input:
        raw_fastq1 = "results/{data_type}/fastq/raw__{sample_name}__R1.fastq.gz",
        raw_fastq2 = "results/{data_type}/fastq/raw__{sample_name}__R2.fastq.gz"
    output:
        fastq1 = "results/{data_type}/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "results/{data_type}/fastq/trim__{sample_name}__R2.fastq.gz",
        metrics = "results/{data_type}/reports/trim_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        adapter1 = lambda wildcards: config['adapter1'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        adapter2 = lambda wildcards: config['adapter2'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimming_quality = lambda wildcards: config['trimming_quality'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimmed_fastqs = config['trimmed_fastqs']
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "trimming", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["process_fastq_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["process_fastq_pe"]["mem_mb"],
        tmp_mb=config["resources"]["process_fastq_pe"]["tmp_mb"],
        qos=config["resources"]["process_fastq_pe"]["qos"]
    shell:
        """
        {{
		if [[ "{params.trimmed_fastqs}" == "True" ]]; then
            printf "\nFastq for {params.sample_name} is already trimmed\n"
            cp {input.raw_fastq1} {output.fastq1}
            cp {input.raw_fastq2} {output.fastq2}
            touch {output.metrics}
        else
            #### Trimming illumina adapters with Cutadapt
            printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
            cutadapt --version
            cutadapt -j {threads} {params.trimming_quality} -a "{params.adapter1}" -A "{params.adapter2}" -o "{output.fastq1}" -p "{output.fastq2}" "{input.raw_fastq1}" "{input.raw_fastq2}" 2>&1 | tee "{output.metrics}"
        fi
        }} 2>&1 | tee -a "{log}"        
        """
        
rule process_fastq_se:
    input:
        raw_fastq = "results/{data_type}/fastq/raw__{sample_name}__R0.fastq.gz"
    output:
        fastq = "results/{data_type}/fastq/trim__{sample_name}__R0.fastq.gz",
        metrics = "results/{data_type}/reports/trim_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        adapter1 = lambda wildcards: config['adapter1'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimming_quality = lambda wildcards: config['trimming_quality'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimmed_fastqs = config['trimmed_fastqs']
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "trimming", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["process_fastq_se"]["threads"]
    resources:
        mem_mb=config["resources"]["process_fastq_se"]["mem_mb"],
        tmp_mb=config["resources"]["process_fastq_se"]["tmp_mb"],
        qos=config["resources"]["process_fastq_se"]["qos"]
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" ]]; then
            printf "\nFastq for {params.sample_name} is already trimmed\n"
            cp {input.raw_fastq} {output.fastq}
        else
            printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
            cutadapt --version
            cutadapt -j {threads} {params.trimming_quality} -a "{params.adapter1}" -o "{output.fastq}" "{input.raw_fastq}" 2>&1 | tee "{output.metrics}"
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule get_available_bam:
    output: 
        bam = temp("results/{data_type}/mapped/copied__{sample_name}.bam")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        bam_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        aligned_bams = config['aligned_bams'],
        exist_bam = lambda wildcards: f"results/{wildcards.data_type}/mapped/final__{wildcards.sample_name}.bam"
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "copy_bam", "either"))
    conda: CONDA_ENV
    threads: config["resources"]["get_available_bam"]["threads"]
    resources:
        mem_mb=config["resources"]["get_available_bam"]["mem_mb"],
        tmp_mb=config["resources"]["get_available_bam"]["tmp_mb"],
        qos=config["resources"]["get_available_bam"]["qos"]
    shell:
        """
        {{
        if [[ "{params.aligned_bams}" == "True" && -e "{params.exist_bam}" ]]; then
            printf "\nFinal bam file already exists for {params.sample_name}\n"
            cp {params.exist_bam} {output.bam}
        elif ls "{params.bam_path}"/*"{params.seq_id}"*.bam 1> /dev/null 2>&1; then
            printf "\nCopying bam file for {params.sample_name} ({params.seq_id} in {params.bam_path})\n"
            samtools sort -@ {threads} -o "{output.bam}" "{params.bam_path}"/*"{params.seq_id}"*.bam
        elif ls "{params.bam_path}"/*"{params.seq_id}"*.sam 1> /dev/null 2>&1; then
            printf "\nCopying and gzipping sam file for {params.sample_name} ({params.seq_id} in {params.bam_path})\n"
            samtools sort -@ {threads} -b -o "{output.bam}" "{params.bam_path}"/*"{params.seq_id}"*.sam 
        else
            printf "Error: No bam or sam file found for {params.sample_name} ({params.seq_id} in {params.bam_path})\n"
        fi
        }} 2>&1 | tee -a "{log}"
        """