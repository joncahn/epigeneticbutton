# function to access logs more easily
def return_log(sample_name, step):
    return os.path.join(REPO_FOLDER,"ChIP","logs",f"tmp__{sample_name}__{step}.log")

def sample_name(s):
    return f"ChIP__{s['line']}__{s['tissue']}__{s['sample_type']}__{s['replicate']}__{s['ref_genome']}"

def get_fastq_inputs(wildcards):
    s = {k: getattr(wildcards, k) for k in ["data_type", "line", "tissue", "sample_type", "replicate", "ref_genome"]}
    name = sample_name(s)
    paired = get_sample_info(wildcards, "paired")
    if paired == "PE":
        return [
            f"ChIP/aligned/{name}__R1.fastq.gz",
            f"ChIP/aligned/{name}__R2.fastq.gz"
        ]
    else:
        return f"ChIP/aligned/{name}.fastq.gz"

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

rule download_fastq_pe:
    output:
        fastq1 = "ChIP/fastq/{sample_name}__R1.fastq.gz",
        fastq2 = "ChIP/fastq/{sample_name}__R2.fastq.gz"
    params:
        seq_id = lambda wildcards: get_sample_info(wildcards, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info(wildcards, "fastq_path"),
        paired = lambda wildcards: get_sample_info(wildcards, "paired"),
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        return_log("{sample_name}", "download_fastq")
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
        
rule download_fastq_se:
    output:
        fastq0 = "ChIP/fastq/{sample_name}.fastq.gz"
    params:
        seq_id = lambda wildcards: get_sample_info(wildcards, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info(wildcards, "fastq_path"),
        paired = lambda wildcards: get_sample_info(wildcards, "paired"),
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        return_log("{sample_name}", "download_fastq")
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
        
rule process_chip_sample:
    input:
        lambda wildcards: get_fastq_inputs(wildcards)
    output:
        chkpt = "ChIP/chkpts/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    params:
        scripts_dir = os.path.join(REPO_FOLDER,"scripts"),
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref_genome),
        data_type = lambda wildcards: wildcards.data_type,
        line = lambda wildcards: wildcards.line,
        tissue = lambda wildcards: wildcards.tissue,
        sample_type = lambda wildcards: wildcards.sample_type,
        replicate = lambda wildcards: wildcards.replicate,
        seq_id = lambda wildcards: get_sample_info(wildcards, 'seq_id'),
        fastq_path = lambda wildcards: get_sample_info(wildcards, 'fastq_path'),
        paired = lambda wildcards: get_sample_info(wildcards, 'paired'),
        mapping_option = config["mapping_option"]
    log:
        return_log("{sample_name}", "process")
    conda:
        CONDA_ENV
    shell:
        """
        cd ChIP/
        qsub ../{params.scripts_dir}/MaizeCode_ChIP_sample.sh \
            -x "ChIP" \
            -d {params.ref_dir} \
            -l {params.line} \
            -t {params.tissue} \
            -m {params.data_type} \
            -r {params.replicate} \
            -i {params.seq_id} \
            -f {params.fastq_path} \
            -p {params.paired} \
            -s "trim" \
            -a {params.mapping_option} | tee {log}
        cd ..
        touch {output.chkpt}
        """
