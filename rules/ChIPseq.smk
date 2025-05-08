# function to access logs more easily
def return_log_chip(sample_name, step):
    return os.path.join(REPO_FOLDER,"ChIP","logs",f"tmp__{sample_name}__{step}.log")

def get_fastq_inputs(wildcards):
    s = {k: getattr(wildcards, k) for k in ["data_type","line", "tissue", "sample_type", "replicate", "ref_genome"]}
    name = sample_name(s)
    paired = get_sample_info(wildcards, "paired")
    if paired == "PE":
        return [
            f"ChIP/fastq/trim__{name}__R1.fastq.gz",
            f"ChIP/fastq/trim__{name}__R2.fastq.gz"
        ]
    else:
        return f"ChIP/fastq/trim__{name}__R0.fastq.gz"

CONDA_ENV=os.path.join(REPO_FOLDER,"envs/chip.yaml")

rule stat_file:
    output:
        stat_file = "ChIP/reports/summary_mapping_stats.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """

rule make_ChIP_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa",
        gff = "genomes/{ref_genome}/temp_{ref_genome}.gff",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        indices = "combined/genomes/{ref_genome}"
    log:
        os.path.join(REPO_FOLDER,"logs","bowtie_index_{ref_genome}.log")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        ### There could be an issue between overlapping indices being built between ChIP and TF, to be resolved
        if ls {output.indices}/*.bt2* 1> /dev/null 2>&1; then
            printf "\nBowtie2 index already exists for {ref_genome}\n" >> {log} 2>&1
        else
            printf "\nBuilding Bowtie2 index for {ref_genome}\n" >> {log} 2>&1
            bowtie2-build --threads {threads} {input.fasta} {output.indices}
        fi
        """

rule get_fastq_pe:
    output:
        fastq1 = "ChIP/fastq/raw__{sample_name}__R1.fastq.gz",
        fastq2 = "ChIP/fastq/raw__{sample_name}__R2.fastq.gz"
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        return_log_chip("{sample_name}", "download_fastq")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        if [[ {params.fastq_path} == "SRA" ]]; then
            printf "\nUsing fasterq-dump for {params.sample_name} ({params.seq_id})\n" >> {log} 2>&1
            fasterq-dump -e {threads} --outdir ChIP/fastq {params.seq_id}
            printf "\n{params.sample_name} ({params.seq_id}) downloaded\nGzipping and renaming files..."
            pigz -p {threads} ChIP/fastq/{params.seq_id}_1.fastq
            mv ChIP/fastq/{params.seq_id}_1.fastq.gz {output.fastq1}
            pigz -p {threads} ChIP/fastq/{params.seq_id}_2.fastq
            mv ChIP/fastq/{params.seq_id}_2.fastq.gz {output.fastq2}
        else
            printf "\nCopying PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n" >> {log} 2>&1
            cp {params.fastq_path}/*{params.seq_id}*R1*q.gz {output.fastq1}
            cp {params.fastq_path}/*{params.seq_id}*R2*q.gz {output.fastq2}
        fi
        """
        
rule get_fastq_se:
    output:
        fastq0 = "ChIP/fastq/raw__{sample_name}__R0.fastq.gz"
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        return_log_chip("{sample_name}", "download_fastq")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        if [[ {params.fastq_path} == "SRA" ]]; then
            printf "\nUsing fasterq-dump for {params.sample_name} ({params.seq_id})\n" >> {log} 2>&1
            fasterq-dump -e {threads} --outdir ChIP/fastq {params.seq_id}
            printf "\n {params.sample_name} ({params.seq_id}) downloaded\nRenaming files..." >> {log} 2>&1
            pigz -p {threads} ChIP/fastq/{params.seq_id}.fastq
            mv ChIP/fastq/{params.seq_id}.fastq.gz {output.fastq0}
        else
            printf "\nCopying SE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n" >> {log} 2>&1
            cp {params.fastq_path}/${params.seq_id}*q.gz {output.fastq0}
        fi
        """

rule process_fastq_pe:
    input:
        raw_fastq1 = "ChIP/fastq/raw__{sample_name}__R1.fastq.gz",
        raw_fastq2 = "ChIP/fastq/raw__{sample_name}__R2.fastq.gz"
    output:
        fastq1 = "ChIP/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "ChIP/fastq/trim__{sample_name}__R2.fastq.gz"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        adapter1 = "AGATCGGAAGAGCACACGTCTGAAC",
        adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGA",
        trimming_quality = config['trimming_quality']
    log:
        return_log_chip("{sample_name}", "trimming")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        #### printf "\nRunning fastQC for {params.sample_name} with fastqc version:\n"
        fastqc --version
        fastqc -o ChIP/reports/ {input.raw_fastq1}
        fastqc -o ChIP/reports/ {input.raw_fastq2}	
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a {params.adapter1} -A {params.adapter2} -o {output.fastq1} -p {output.fastq2} {input.raw_fastq1} {input.raw_fastq2} |& tee {log}
		#### Removing untrimmed fastq
		rm -f {input.raw_fastq1} {input.raw_fastq2}
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for {params.sample_name}\n"
		fastqc -o ChIP/reports/ {output.fastq1}
		fastqc -o ChIP/reports/ {output.fastq2}
        """
        
rule process_fastq_se:
    input:
        raw_fastq = "ChIP/fastq/raw__{sample_name}__R0.fastq.gz"
    output:
        fastq = "ChIP/fastq/trim__{sample_name}__R0.fastq.gz"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        adapter1 = "AGATCGGAAGAGCACACGTCTGAAC",
        trimming_quality = config['trimming_quality']
    log:
        return_log_chip("{sample_name}", "trimming")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        ### QC of the raw reads with FastQC
        printf "\nRunning fastQC for {params.sample_name} with fastqc version:\n"
        fastqc --version
        fastqc -o ChIP/reports/ {input.raw_fastq}
		#### Trimming illumina adapters with Cutadapt
		printf "\nTrimming Illumina adapters for {params.sample_name} with cutadapt version:\n"
		cutadapt --version
		cutadapt -j {threads} {params.trimming_quality} -a {params.adapter1} -o {output.fastq} {input.raw_fastq} |& tee {log}
		#### Removing untrimmed fastq
		rm -f {input.raw_fastq}
		#### FastQC on trimmed data
		printf "\nRunning fastQC on trimmed files for {params.sample_name}\n"
		fastqc -o ChIP/reports/ {output.fastq}
        """

rule bowtie2_map_pe:
    input:
        fastq1 = "ChIP/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "ChIP/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = "combined/genomes/{ref_genome}"
    output:
        samfile = "ChIP/mapped/mapped__{sample_name}.bam",
        metrics = "ChIP/reports/bt2__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref = lambda wildcards: ref_genome,
        map_option = lambda wildcards: config['mapping_option'],
        mapping_params = lambda wildcards: config['mapping'][config['mapping_option']]['map_pe']    
    log:
        return_log_chip("{sample_name}", "mapping_on_{ref_genome}")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nMaping {params.sample_name} to {params.ref} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} --met-file {output.metrics} -x {input.indices} -1 {input.fastq1} -2 {input.fastq2} -S {output.sam} |& tee {log}
        """    
        
rule bowtie2_map_se:
    input:
        fastq = "ChIP/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = "combined/genomes/{ref_genome}"
    output:
        samfile = "ChIP/mapped/mapped__{sample_name}.bam",
        metrics = "ChIP/reports/bt2__{sample_name}.txt",
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref = lambda wildcards: ref_genome,
        map_option = lambda wildcards: config['mapping_option'],
        mapping_params = lambda wildcards: config['mapping'][config['mapping_option']]['map_se']    
    log:
        return_log_chip("{sample_name}", "mapping_on_{ref_genome}")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nMaping {params.sample_name} to {params.ref} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} --met-file {output.metrics} -x {input.indices} -U {input.fastq} -S {output.sam} |& tee {log}
        """        
        
rule check_pair:
    input:
        lambda wildcards: get_fastq_inputs(wildcards)
    output:
        touch = "ChIP/chkpts/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    shell:
        """
        touch {output.touch}
        """
     
# rule process_chip_sample:
    # input:
        # lambda wildcards: get_fastq_inputs(wildcards)
    # output:
        # chkpt = "ChIP/chkpts/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    # params:
        # scripts_dir = os.path.join(REPO_FOLDER,"scripts"),
        # ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref_genome),
        # data_type = lambda wildcards: wildcards.data_type,
        # line = lambda wildcards: wildcards.line,
        # tissue = lambda wildcards: wildcards.tissue,
        # sample_type = lambda wildcards: wildcards.sample_type,
        # replicate = lambda wildcards: wildcards.replicate,
        # seq_id = lambda wildcards: get_sample_info(wildcards, 'seq_id'),
        # fastq_path = lambda wildcards: get_sample_info(wildcards, 'fastq_path'),
        # paired = lambda wildcards: get_sample_info(wildcards, 'paired'),
        # mapping_option = config["mapping_option"],
        # log_path = lambda wildcards: return_log(f"ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}","process")
    # log:
        # "{params.log_path}"
    # conda:
        # CONDA_ENV
    # shell:
        # """
        # cd ChIP/
        # qsub ../{params.scripts_dir}/MaizeCode_ChIP_sample.sh \
            # -x "ChIP" \
            # -d {params.ref_dir} \
            # -l {params.line} \
            # -t {params.tissue} \
            # -m {params.data_type} \
            # -r {params.replicate} \
            # -i {params.seq_id} \
            # -f {params.fastq_path} \
            # -p {params.paired} \
            # -s "trim" \
            # -a {params.mapping_option} | tee {log}
        # cd ..
        # touch {output.chkpt}
        # """
