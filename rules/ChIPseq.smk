# function to access logs more easily
def return_log_chip(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"ChIP","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def get_inputs_chip(wildcards):
    s = {k: getattr(wildcards, k) for k in ["data_type","line", "tissue", "sample_type", "replicate", "ref_genome"]}
    name = sample_name(s)
    paired = get_sample_info(wildcards, "paired")
    if paired == "PE":
        return f"ChIP/logs/process_pe_sample__{name}.log"
    else:
        return f"ChIP/logs/process_se_sample__{name}.log"
        
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/chip.yaml")

rule stat_file_chip:
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
        indices = "genomes/{ref_genome}"
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

rule bowtie2_map_pe:
    input:
        fastq1 = "ChIP/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "ChIP/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"combined/genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}"
    output:
        samfile = "ChIP/mapped/mapped_pe__{sample_name}.sam",
        metrics = "ChIP/reports/bt2_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: config['mapping_option'],
        mapping_params = lambda wildcards: config['mapping'][config['mapping_option']]['map_pe']    
    log:
        return_log_chip("{sample_name}", "mapping", "PE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nMaping {params.sample_name} to {params.ref_genome} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} --met-file {log} -x {input.indices} -1 {input.fastq1} -2 {input.fastq2} -S {output.sam} |& tee {output.metrics}
        """    
        
rule bowtie2_map_se:
    input:
        fastq = "ChIP/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"combined/genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}"
    output:
        samfile = "ChIP/mapped/mapped_se__{sample_name}.sam",
        metrics = "ChIP/reports/bt2_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: config['mapping_option'],
        mapping_params = lambda wildcards: config['mapping'][config['mapping_option']]['map_se']    
    log:
        return_log_chip("{sample_name}", "mapping", "SE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nMaping {params.sample_name} to {params.ref} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} --met-file {log} -x {input.indices} -U {input.fastq} -S {output.sam} |& tee {output.metrics}
        """

rule filter_results_pe:
    input:
        samfile = "ChIP/mapped/mapped_pe__{sample_name}.sam"
    output:
        bamfile = "ChIP/mapped/{sample_name}.bam",
        metrics_dup = "ChIP/reports/markduppe__{sample_name}.txt",
        metrics_flag = "ChIP/reports/flagstatpe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        map_option = lambda wildcards: config['mapping_option'],
        filtering_params = lambda wildcards: config['mapping'][config['mapping_option']]['filter']    
    log:
        return_log_chip("{sample_name}", "filtering", "PE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o ChIP/mapped/temp1_{sample_name}.bam {input.samfile}
        rm -f {input.samfile}
        samtools fixmate -@ {threads} -m ChIP/mapped/temp1_{sample_name}.bam ChIP/mapped/temp2_{sample_name}.bam
        samtools sort -@ {threads} -o ChIP/mapped/temp3_{sample_name}.bam ChIP/mapped/temp2_{sample_name}.bam
        samtools markdup -r -s -f {output.metrics_dup} -@ {threads} ChIP/mapped/temp3_{sample_name}.bam {output.bamfile}
        samtools index -@ {threads} {output.bamfile}
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} {output.bamfile} > {output.metrics_flag}
        rm -f ChIP/mapped/temp*_{sample_name}.bam
        """

rule filter_results_se:
    input:
        samfile = "ChIP/mapped/mapped_se__{sample_name}.sam"
    output:
        bamfile = "ChIP/mapped/{sample_name}.bam",
        metrics_dup = "ChIP/reports/markdupse__{sample_name}.txt",
        metrics_flag = "ChIP/reports/flagstatse__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        map_option = lambda wildcards: config['mapping_option'],
        filtering_params = lambda wildcards: config['mapping'][config['mapping_option']]['filter']    
    log:
        return_log_chip("{sample_name}", "filtering", "SE")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o ChIP/mapped/temp1_{sample_name}.bam {input.samfile}
        rm -f {input.samfile}
        samtools sort -@ {threads} -o ChIP/mapped/temp2_{sample_name}.bam ChIP/mapped/temp1_{sample_name}.bam
        samtools markdup -r -s -f {output.metrics_dup} -@ {threads} ChIP/mapped/temp2_{sample_name}.bam {output.bamfile}
        samtools index -@ {threads} {output.bamfile}
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} {output.bamfile} > {output.metrics_flag}
        rm -f ChIP/mapped/temp*_{sample_name}.bam
        """

rule make_statistics_file_pe:
    input:
        stat_file = "ChIP/reports/summary_mapping_stats.txt",
        metrics_trim = "ChIP/reports/trim_pe__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        metrics_map = "ChIP/reports/bt2_pe__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        logs = lambda wildcards: [ return_log_chip(sample_name(wildcards), step, get_sample_info(wildcards, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "ChIP/logs/process_pe_sample__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" {input.metrics_trim} | awk '{print $NF}' | sed 's/,//g')
        filt=$(grep "reads" {input.metrics_map} | awk '{print $1}')
        multi=$(grep "aligned concordantly >1 times" {input.metrics_map} | awk '{print $1}')
        single=$(grep "aligned concordantly exactly 1 time" {input.metrics_map} | awk '{print $1}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={line} -v t={tissue} -v m={sample_type} -v r={rep} -v g={ref_genome} -v a=${tot} -v b=${filt} -v c=${allmap} -v d=${single} 'BEGIN {print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' >> ChIP/reports/summary_mapping_stats.txt
        cat {input.logs} > {output.log}
        rm -f {input.logs}
        """

rule make_statistics_file_se:
    input:
        stat_file = "ChIP/reports/summary_mapping_stats.txt",
        metrics_trim = "ChIP/reports/trim_se__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        metrics_map = "ChIP/reports/bt2_se__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        logs = lambda wildcards: [ return_log_chip(sample_name(wildcards), step, get_sample_info(wildcards, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "ChIP/logs/process_se_sample__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total reads processed:" {input.metrics_trim} | awk '{print $NF}' | sed 's/,//g')
        filt=$(grep "reads" {input.metrics_map} | awk '{print $1}')
        multi=$(grep "aligned >1 times" {input.metrics_map} | awk '{print $1}')
        single=$(grep "aligned exactly 1 time" {input.metrics_map} | awk '{print $1}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={line} -v t={tissue} -v m={sample_type} -v r={rep} -v g={ref_genome} -v a=${tot} -v b=${filt} -v c=${allmap} -v d=${single} 'BEGIN {print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}' >> ChIP/reports/summary_mapping_stats.txt
        cat {input.logs} > {output.log}
        rm -f {input.logs}
        """
        
rule check_pair_chip:
    input:
        lambda wildcards: get_inputs_chip(wildcards)
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
