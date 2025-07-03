# CONDA_ENV_sRNA=os.path.join(REPO_FOLDER,"envs/epibutton_srna.yaml")

# function to access logs more easily
def return_log_smallrna(sample_name, step, size):
    return os.path.join(REPO_FOLDER,"results","sRNA","logs",f"tmp__{sample_name}__{step}__{size}.log")

def define_input_file_for_shortstack(sample_name):
    paired = get_sample_info_from_name(sample_name, samples, 'paired')
    if paired == "SE":
        return "filtered__{sample_name}__R0" if config['structural_rna_depletion'] else "trim__{sample_name}__R0"
    
def define_final_srna_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    analysis_name = config['analysis_name']
    srna_min = config['srna_min_size']
    srna_max = config['srna_max_size']
    map_files = []
    bigwig_files = []
    qc_files = []
    deg_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'sRNA') & (samples['ref_genome'] == ref_genome) ].copy()
    
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')        
        map_files.append(f"results/sRNA/reports/sizes_stats__{sname}.txt")
        qc_files.append(f"results/sRNA/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file
        qc_files.append(f"results/sRNA/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
        
        for size in range(srna_min, srna_max + 1):
            bigwig_files.append(f"results/sRNA/tracks/{sname}__{size}nt__plus.bw")
            bigwig_files.append(f"results/sRNA/tracks/{sname}__{size}nt__minus.bw")
        
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'sRNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            for size in range(srna_min, srna_max + 1):
                bigwig_files.append(f"results/sRNA/tracks/{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}__{size}nt__plus.bw")
                bigwig_files.append(f"results/sRNA/tracks/{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}__{size}nt__minus.bw")
            
    results = map_files
    
    if qc_option == "all":
        results += qc_files
        
    if analysis:
        results += bigwig_files

    return results

rule filter_structural_rna:
    input:
        fastq = "results/sRNA/fastq/trim__{sample_name}__R0.fastq.gz",
        fasta = config['structural_rna_fafile']
    output:
        filtered_fastq = temp("results/sRNA/fastq/filtered__{sample_name}__R0.fastq"),
        gzipped_fastq = "results/sRNA/fastq/filtered__{sample_name}__R0.fastq.gz"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_smallrna("{sample_name}", "filter_structural_rna", "all"))
    conda: CONDA_ENV
    threads: config["resources"]["filter_structural_rna"]["threads"]
    resources:
        mem=config["resources"]["filter_structural_rna"]["mem"],
        tmp=config["resources"]["filter_structural_rna"]["tmp"]
    shell:
        """
        {{
        if [[ ! -d genomes/structural_RNAs/{params.ref_genome}/bt2_index ]]; then
            if [[ {input.fasta} =~ \.gz$ ]]; then
                printf "Generating bowtie2 index for structural RNAs of {params.ref_genome} using gzipped file: {input.fasta}\n"
                mkdir -p genomes/structural_RNAs/{params.ref_genome}/bt2_index
                gunzip {input.fasta} -c > "genomes/structural_RNAs/{params.ref_genome}/temp.fa"
                bowtie2-build --threads {threads} "genomes/structural_RNAs/{params.ref_genome}/temp.fa" "genomes/structural_RNAs/{params.ref_genome}/bt2_index"
                rm -f "genomes/structural_RNAs/{params.ref_genome}/temp.fa"
            else
                printf "Generating bowtie2 index for structural RNAs of {params.ref_genome} using file: {input.fasta}\n"
                mkdir -p genomes/structural_RNAs/{params.ref_genome}/bt2_index
                bowtie2-build --threads {threads} "{input.fasta}" "genomes/structural_RNAs/{params.ref_genome}/bt2_index"
            fi
        fi
        bowtie2 --very-sensitive -p {threads} -x genomes/structural_RNAs/{params.ref_genome}/bt2_index -U {input.fastq} | samtools view -@ {threads} -f 0x4 | samtools fastq -@ {threads} > {output.filtered_fastq}
        pigz -p {threads} {output.filtered_fastq} -c > {output.gzipped_fastq}
        }} 2>&1 | tee -a "{log}"
        """

rule dispatch_srna_fastq:
    input:
        fastq = lambda wildcards: f"results/sRNA/fastq/{define_input_file_for_shortstack(wildcards.sample_name)}.fastq.gz"
    output:
        fastq_file = temp("results/sRNA/mapped/clean__{sample_name}.fastq.gz")
    conda: CONDA_ENV
    localrule: True
    shell:
        """
        cp {input.fastq} {output.fastq_file}
        """

rule shortstack_map:
    input:
        fastq = "results/sRNA/mapped/clean__{sample_name}.fastq.gz",
        fasta = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa"
    output:
        count_file = "results/sRNA/mapped/{sample_name}/Results.txt",
        bam_file = temp("results/sRNA/mapped/{sample_name}/clean__{sample_name}.bam"),
        bai_file = temp("results/sRNA/mapped/{sample_name}/clean__{sample_name}.bam.bai")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        srna_params = config['srna_mapping_params']
    log:
        temp(return_log_smallrna("{sample_name}", "mapping_shortstack", "all"))
    conda: CONDA_ENV
    threads: config["resources"]["shortstack_map"]["threads"]
    resources:
        mem=config["resources"]["shortstack_map"]["mem"],
        tmp=config["resources"]["shortstack_map"]["tmp"]
    shell:
        """
        {{
        rm -rf sRNA/mapped/{params.sample_name}
        printf "\nMapping {params.sample_name} to {params.ref_genome} with Shortstack version:\n"
        ShortStack --version
        ShortStack --readfile {input.fastq} --genomefile {input.fasta} --threads {threads} {params.srna_params} --outdir results/sRNA/mapped/{params.sample_name}
        samtools index -@ {threads} {output.bam_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule make_srna_size_stats:
    input:
        bamfile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}.bam",
        baifile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}.bam.bai"
    output:
        report = "results/sRNA/reports/sizes_stats__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        temp(return_log_smallrna("{sample_name}", "make_srna_stats", "all"))
    conda: CONDA_ENV
    threads: config["resources"]["make_srna_size_stats"]["threads"]
    resources:
        mem=config["resources"]["make_srna_size_stats"]["mem"],
        tmp=config["resources"]["make_srna_size_stats"]["tmp"]
    shell:
        """
        {{
        printf "\nGetting stats for {params.sample_name}\n"
        printf "Sample\tType\tSize\tCount\n" > {output.report}
        zcat results/sRNA/fastq/trim__{params.sample_name}__R0.fastq.gz | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"trimmed",$2,$1}}' >> "{output.report}"
        printf "\nGetting filtered stats for {params.sample_name}\n"
        if [[ -s results/sRNA/fastq/filtered__{params.sample_name}__R0.fastq.gz ]]; then
            zcat results/sRNA/fastq/filtered__{params.sample_name}__R0.fastq.gz | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"filtered",$2,$1}}' >> "{output.report}"
        fi
        samtools view {input.bamfile} | awk '$2==0 || $2==16 {{print length($10)}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"mapped",$2,$1}}' >> "{output.report}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_size_srna_sample:
    input:
        bamfile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}.bam",
        baifile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}.bam.bai"
    output:
        filtered_file = "results/sRNA/mapped/sized__{size}nt__{sample_name}.bam"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        size = lambda wildcards: wildcards.size
    log:
        temp(return_log_smallrna("{sample_name}", "filter_size_srna", "{size}"))
    conda: CONDA_ENV
    threads: config["resources"]["filter_size_srna_sample"]["threads"]
    resources:
        mem=config["resources"]["filter_size_srna_sample"]["mem"],
        tmp=config["resources"]["filter_size_srna_sample"]["tmp"]
    shell:
        """
        {{
        printf "Filtering only {params.size} nucleotides sRNAs for {params.sample_name}\n"
        samtools view -h {input.bamfile} | awk -v n={params.size} '(length($10) == n) || $1 ~ /^@/' | samtools view -bS - > {output.filtered_file}
        samtools index -@ {threads} {output.filtered_file}
        }} 2>&1 | tee -a "{log}"
        """

rule merging_srna_replicates:
    input:
        bamfiles = lambda wildcards: [ f"results/sRNA/mapped/sized__{wildcards.size}nt__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.bam" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        tempfile = temp("results/sRNA/mapped/temp__{size}nt__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"),
        mergefile = "results/sRNA/mapped/merged__{size}nt__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis'),
        size = lambda wildcards: wildcards.size
    log:
        temp(return_log_smallrna("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging_srna_reps", "{size}"))
    conda: CONDA_ENV
    threads: config["resources"]["merging_srna_replicates"]["threads"]
    resources:
        mem=config["resources"]["merging_srna_replicates"]["mem"],
        tmp=config["resources"]["merging_srna_replicates"]["tmp"]
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname} {params.size}\n"
		samtools merge -@ {threads} {output.tempfile} {input.bamfiles}
		samtools sort -@ {threads} -o {output.mergefile} {output.tempfile}
		samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """

rule make_srna_stranded_bigwigs:
    input: 
        bamfile = lambda wildcards: f"results/sRNA/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'sized'}__{wildcards.size}nt__{wildcards.sample_name}.bam"
    output:
        bw_plus = "results/sRNA/tracks/{sample_name}__{size}nt__plus.bw",
        bw_minus = "results/sRNA/tracks/{sample_name}__{size}nt__minus.bw"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        size = lambda wildcards: wildcards.size,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_smallrna("{sample_name}", "making_bigiwig", "{size}"))
    conda: CONDA_ENV
    threads: config["resources"]["make_srna_stranded_bigwigs"]["threads"]
    resources:
        mem=config["resources"]["make_srna_stranded_bigwigs"]["mem"],
        tmp=config["resources"]["make_srna_stranded_bigwigs"]["tmp"]
    shell:
        """
        {{
        printf "Getting stranded coverage for {params.sample_name} {params.size}nt\n"
		bamCoverage --filterRNAstrand reverse -bs 1 -p {threads} --normalizeUsing CPM -b {input.bamfile} -o {output.bw_plus}
		bamCoverage --filterRNAstrand forward -bs 1 -p {threads} --normalizeUsing CPM -b {input.bamfile} -o {output.bw_minus}
        }} 2>&1 | tee -a "{log}"
        """

rule all_srna:
    input:
        final = lambda wildcards: define_final_srna_output(wildcards.ref_genome)
    output:
        touch = "results/sRNA/chkpts/sRNA_analysis__{analysis_name}__{ref_genome}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """