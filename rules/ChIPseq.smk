# function to access logs more easily
def return_log(line, tissue, sample_type, rep, ref, step):
    return os.path.join(REPO_FOLDER,"ChIP","logs",f"tmp_chip_{step}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log")

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
        fasta = "genomes/{ref}/temp_{ref}.fa",
        gff = "genomes/{ref}/temp_{ref}.gff",
        chrom_sizes = "genomes/{ref}/chrom.sizes"
      output:
        indices = "combined/genomes/{ref}"
      log:
        os.path.join(REPO_FOLDER,"logs",f"bowtie_index_{ref}.log")
      threads: workflow.cores
      shell:
        """
        ### There could be an issue between overlapping indices being built between ChIP and TF, to be resolved
        if ls {output.indices}/*.bt2* 1> /dev/null 2>&1; then
            printf "\nBowtie2 index already exists for {wildcards.ref}\n" >> {log} 2>&1
        else
            printf "\nBuilding Bowtie2 index for {wildcards.ref}\n" >> {log} 2>&1
            bowtie2-build --threads {threads} {input.fasta} {output.indices}
        fi
        """

rule download_fastq:
    output:
        touch = "ChIP/chkpts/ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.rep}__{wildcards.ref_genome}.done"
        lambda wildcards: [
            fastq1 = f"ChIP/fastq/ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.rep}__{wildcards.ref_genome}__R1.fastq.gz",
            fastq2 = f"ChIP/fastq/ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.rep}__{wildcards.ref_genome}__R2.fastq.gz"
        ] if get_sample_info(wildcards, "paired") == "PE" else
            fastq0 = f"ChIP/fastq/ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.rep}__{wildcards.ref_genome}.fastq.gz"
    params:
        data_type = lambda wildcards: wildcards.data_type,
        line = lambda wildcards: wildcards.line,
        tissue = lambda wildcards: wildcards.tissue,
        sample_type = lambda wildcards: wildcards.sample_type,
        replicate = lambda wildcards: wildcards.replicate,
        ref_genome = lambda wildcards: wildcards.ref_genome,
        seq_id = lambda wildcards: get_sample_info(wildcards, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info(wildcards, "fastq_path"),
        paired = lambda wildcards: get_sample_info(wildcards, "paired")
        sample_name = lambda wildcards: f"ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.rep}__{wildcards.ref_genome}"
    log:
        return_log("{line}", "{tissue}", "{sample_type}", "{rep}", "{ref}", "download_fastq")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        if [[ {params.paired} == "PE" ]]; then
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
        else
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
        fi
        touch {output.touch}
        """
        
rule process_chip_sample:
    input:
        touch = "ChIP/chkpts/ChIP__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.rep}__{wildcards.ref_genome}.done"
    output:
        chkpt = "ChIP/chkpts/process__ChIP__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    params:
        scripts_dir = os.path.join(REPO_FOLDER,"scripts"),
        ref_dir = lambda wildcards: os.path.join(REF_PATH, get_sample_info(wildcards, 'ref_genome')),
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
        return_log("{line}", "{tissue}", "{sample_type}", "{rep}", "{ref}", "process")
    conda:
        CONDA_ENV
    shell:
        """
        cd ChIP/
        qsub ../{params.scripts_dir}/MaizeCode_{params.env}_sample.sh \
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
