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
        
def assign_mapping_paired(wildcards, rulename, outputfile):
    sname = sample_name(wildcards)
    paired = get_sample_info_from_name(sname,'paired')
    if paired == "PE":
        rule_obj = getattr(rules, f"{rulename}_pe")
    else:
        rule_obj = getattr(rules, f"{rulename}_pe")
    return getattr(rule_obj.output, outputfile).format(sample_name=sname)
        
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/chip.yaml")

rule stat_file_chip:
    output:
        stat_file = f"ChIP/reports/summary_mapping_stats_{analysis_name}.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """

rule make_bt2_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa",
        gff = "genomes/{ref_genome}/temp_{ref_genome}.gff",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        indices = directory("genomes/{ref_genome}/bt2_index")
    log:
        os.path.join(REPO_FOLDER,"logs","bowtie_index_{ref_genome}.log")
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nBuilding Bowtie2 index for {wildcards.ref_genome}\n"
        mkdir genomes/{wildcards.ref_genome}/bt2_index
        bowtie2-build --threads {threads} "{input.fasta}" "{output.indices}/{wildcards.ref_genome}"
        }} 2>&1 | tee -a "{log}"
        """

rule bowtie2_map_pe:
    input:
        fastq1 = "ChIP/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "ChIP/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/bt2_index"
    output:
        samfile = temp("ChIP/mapped/mapped_pe__{sample_name}.sam"),
        metrics = "ChIP/reports/bt2_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: config['chip_mapping_option'],
        mapping_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['map_pe']    
    log:
        temp(return_log_chip("{sample_name}", "mapping", "PE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nMaping {params.sample_name} to {params.ref_genome} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} -x "{input.indices}/{params.ref_genome}" -1 "{input.fastq1}" -2 "{input.fastq2}" -S "{output.samfile}" 2>&1 | tee "{output.metrics}"
        }} 2>&1 | tee -a "{log}"
        """    
        
rule bowtie2_map_se:
    input:
        fastq = "ChIP/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/bt2_index"
    output:
        samfile = temp("ChIP/mapped/mapped_se__{sample_name}.sam"),
        metrics = "ChIP/reports/bt2_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: config['chip_mapping_option'],
        mapping_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['map_se']    
    log:
        temp(return_log_chip("{sample_name}", "mapping", "SE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nMaping {params.sample_name} to {params.ref_genome} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} -x "{input.indices}/{params.ref_genome}" -U "{input.fastq}" -S "{output.samfile}" 2>&1 | tee "{output.metrics}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_chip_pe:
    input:
        samfile = "ChIP/mapped/mapped_pe__{sample_name}.sam"
    output:
        bamfile = "ChIP/mapped/{sample_name}.bam",
        metrics_dup = "ChIP/reports/markdup_pe__{sample_name}.txt",
        metrics_flag = "ChIP/reports/flagstat_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        map_option = lambda wildcards: config['chip_mapping_option'],
        filtering_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['filter']    
    log:
        temp(return_log_chip("{sample_name}", "filtering", "PE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o "ChIP/mapped/temp1_{sample_name}.bam" "{input.samfile}"
        rm -f "{input.samfile}"
        samtools fixmate -@ {threads} -m "ChIP/mapped/temp1_{sample_name}.bam" "ChIP/mapped/temp2_{sample_name}.bam"
        samtools sort -@ {threads} -o "ChIP/mapped/temp3_{sample_name}.bam" "ChIP/mapped/temp2_{sample_name}.bam"
        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} "ChIP/mapped/temp3_{sample_name}.bam" "{output.bamfile}"
        samtools index -@ {threads} "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f ChIP/mapped/temp*"_{sample_name}.bam"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_chip_se:
    input:
        samfile = "ChIP/mapped/mapped_se__{sample_name}.sam"
    output:
        bamfile = "ChIP/mapped/{sample_name}.bam",
        metrics_dup = "ChIP/reports/markdup_se__{sample_name}.txt",
        metrics_flag = "ChIP/reports/flagstat_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        map_option = lambda wildcards: config['chip_mapping_option'],
        filtering_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['filter']    
    log:
        temp(return_log_chip("{sample_name}", "filtering", "SE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o "ChIP/mapped/temp1_{sample_name}.bam" "{input.samfile}"
        rm -f "{input.samfile}"
        samtools sort -@ {threads} -o "ChIP/mapped/temp2_{sample_name}.bam" "ChIP/mapped/temp1_{sample_name}.bam"
        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} "ChIP/mapped/temp2_{sample_name}.bam" "{output.bamfile}"
        samtools index -@ {threads} "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f ChIP/mapped/temp*"_{sample_name}.bam"
        }} 2>&1 | tee -a "{log}"
        """

rule make_chip_stats_pe:
    input:
        stat_file = f"ChIP/reports/summary_mapping_stats_{analysis_name}.txt",
        metrics_trim = "ChIP/reports/trim_pe__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        metrics_map = "ChIP/reports/bt2_pe__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        logs = lambda wildcards: [ return_log_chip(sample_name(wildcards), step, get_sample_info(wildcards, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "ChIP/logs/process_pe_sample__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned concordantly >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned concordantly exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={wildcards.line} -v t={wildcards.tissue} -v m={wildcards.sample_type} -v r={wildcards.replicate} -v g={wildcards.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{input.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule make_chip_stats_se:
    input:
        stat_file = f"ChIP/reports/summary_mapping_stats_{analysis_name}.txt",
        metrics_trim = "ChIP/reports/trim_se__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        metrics_map = "ChIP/reports/bt2_se__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        logs = lambda wildcards: [ return_log_chip(sample_name(wildcards), step, get_sample_info(wildcards, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "ChIP/logs/process_se_sample__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total reads processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={wildcards.line} -v t={wildcards.tissue} -v m={wildcards.sample_type} -v r={wildcards.replicate} -v g={wildcards.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{input.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """
        
rule check_pair_chip:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "make_chip_stats", "log")
    output:
        touch = "ChIP/chkpts/process__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    shell:
        """
        touch {output.touch}
        """

rule map_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "mapping", "bamfile")
    output:
        "ChIP/mapped/{sample_name}.bam"
    
rule make_coverage_chip:
    input: 
        bamfile = "ChIP/mapped/{sample_name}.bam"
    output:
        bigwigcov = "ChIP/tracks/coverage_{sample_name}.bw"
    params:
        binsize = config['chip_tracks']['binsize']
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        bamCoverage -b {input.bamfile} -o {output.bigwigcov} -bs {params.binsize} -p {threads}
        """
            
rule merging replicates:
    input:
        bamfiles = expand("ChIP/mapped/{sample_name}.bam", 
            sample_name = lambda wildcards: analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []))
    output:
        mergefile = "ChIP/mapped/merged_{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}.bam"
    params:
        sample_name = "{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}"
    log:
        temp(return_log_chip("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging", "merged"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sample_name}\n"
		samtools merge -@ {threads} ChIP/mapped/temp_{params.sample_name}.bam {input.bamfiles}
		samtools sort -@ {threads} -o {output.mergefile} ChIP/mapped/temp_{params.sample_name}.bam
		rm -f ChIP/mapped/temp_{params.sample_name}.bam
		samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """