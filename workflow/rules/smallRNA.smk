CONDA_ENV_SRNA=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_srna.yaml")

def return_log_smallrna(sample_name, step, size):
    return os.path.join(REPO_FOLDER,"results","sRNA","logs",f"tmp__{sample_name}__{step}__{size}.log")

def define_input_file_for_structural(sample_name):
    paired = get_sample_info_from_name(sample_name, samples, 'paired')
    nextflex_v3 = config['nextflex_v3_deduplication']
    if paired == "SE":
        return "deduplicated__{sample_name}__R0" if nextflex_v3 else "trim__{sample_name}__R0"

def get_bt1_indices(wildcards):
    ref_genome = parse_sample_name(wildcards.sample_name)['ref_genome']
    genomesize = float(config[config[ref_genome]['species']]['genomesize'])
    if genomesize > 4e9:
        return multiext(f"genomes/{ref_genome}/{ref_genome}.fa", ".1.ebwtl", ".2.ebwtl",".3.ebwtl",".4.ebwtl",".rev.1.ebwtl",".rev.2.ebwtl")
    else: 
        return multiext(f"genomes/{ref_genome}/{ref_genome}.fa", ".1.ebwt", ".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt")
        
def define_input_file_for_shortstack(sample_name):
    paired = get_sample_info_from_name(sample_name, samples, 'paired')
    rna_depletion = config['structural_rna_depletion']
    netflex_v3 = config['nextflex_v3_deduplication']
    if paired == "SE":
        if rna_depletion:
            return "filtered__{sample_name}__R0" 
        elif nextflex_v3:
            return "deduplicated__{sample_name}__R0"
        else:
            return "trim__{sample_name}__R0"
    elif paired == "PE":
        raise ValueError(
            "Paired-end small RNA is not yet a feature."
            f"{wildcards.sample_name} is set to be PE."
        )

def define_input_for_grouped_analysis(ref_genome):
    bamfiles = []
    filtered_rep_samples = samples[ (samples['env'] == 'sRNA') & (samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        bamfiles.append(f"results/sRNA/mapped/{sname}/clean__{sname}_condensed.bam")
    
    return bamfiles
    
def define_srna_target_file(wildcards):
    tname = config['srna_target_file_label']
    if wildcards.target_name == "new_clusters":
        return []
    elif wildcards.target_name == "all_genes":
        return f"results/combined/bedfiles/{wildcards.ref_genome}__all_genes.bed"
    elif wildcards.target_name == "protein_coding_genes":
        return f"results/combined/bedfiles/{wildcards.ref_genome}__protein_coding_genes.bed"
    elif wildcards.target_name == "all_TEs":
        return f"genomes/{wildcards.ref_genome}/{wildcards.ref_genome}__TE_file.bed"
    elif wildcards.target_name == tname:
        return config['srna_target_file']
    else:
        raise ValueError(   
            f"{wildcards.target_name} does not match possible files." 
            "It can be 'all_genes', 'protein_coding_genes' or the value of "
            "'srna_target_file_name' in the config file"
        )

def define_final_srna_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    te_analysis = config['te_analysis']
    analysis_name = config['analysis_name']
    srna_min = config['srna_min_size']
    srna_max = config['srna_max_size']
    trimmed_fastqs = config['trimmed_fastqs']
    rna_depletion = config['structural_rna_depletion']
    nextflex_v3 = config['nextflex_v3_deduplication']
    map_files = []
    bigwig_files = []
    qc_files = []
    deg_files = []
    analysis_files = []
    te_analysis_files = []
    
    filtered_rep_samples = samples[ (samples['env'] == 'sRNA') & (samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        qc_files.append(f"results/sRNA/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed and potentially filtered fastq files
        map_files.append(f"results/sRNA/reports/sizes_stats__{sname}.txt")
        if not trimmed_fastqs:
            qc_files.append(f"results/sRNA/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw fastq file
            if nextflex_v3:
                qc_files.append(f"results/sRNA/reports/deduplicated__{sname}__R0_fastqc.html") # fastqc of deduplicated (after nextflex v3) fastq file
            if rna_depletion:
                qc_files.append(f"results/sRNA/reports/filtered__{sname}__R0_fastqc.html") # fastqc of structural RNA depleted fastq file
        
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
    
    if len(filtered_analysis_samples) >= 2 and any(len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2 for _, row in filtered_analysis_samples.iterrows()):
        analysis_files.append(f"results/sRNA/chkpts/calling_differential_sRNA_clusters__{analysis_name}__{ref_genome}__on_new_clusters.done")
        analysis_files.append(f"results/sRNA/chkpts/calling_differential_sRNA_clusters__{analysis_name}__{ref_genome}__on_all_genes.done")
        te_analysis_files.append(f"results/sRNA/chkpts/calling_differential_sRNA_clusters__{analysis_name}__{ref_genome}__on_all_TEs.done")
    elif len(filtered_analysis_samples) >= 1:
        analysis_files.append(f"results/sRNA/clusters/{analysis_name}__{ref_genome}__on_new_clusters/Counts.txt")
        analysis_files.append(f"results/sRNA/clusters/{analysis_name}__{ref_genome}__on_all_genes/Counts.txt")
        te_analysis_files.append(f"results/sRNA/clusters/{analysis_name}__{ref_genome}__on_all_TEs/Counts.txt")
    
    results = map_files + bigwig_files

    if qc_option == "all":
        results += qc_files
        
    if analysis:
        results += analysis_files
   
    if te_analysis:
        results += te_analysis_files

    return results

rule deduplicate_srna_nextflexv3:
    input:
        fastq = "results/sRNA/fastq/trim__{sample_name}__R0.fastq.gz"
    output:
        collapse_folder = temp(directory("results/sRNA/fastq/collapsed_{sample_name}")),
        collapsed_fastq = temp("results/sRNA/fastq/collapsed_{sample_name}/trim__{sample_name}__R0_trimmed.fastq"),
        deduplicated_fastq = temp("results/sRNA/fastq/deduplicated__{sample_name}__R0.fastq"),
        gzipped_fastq = temp("results/sRNA/fastq/deduplicated__{sample_name}__R0.fastq.gz"),
        report = "results/sRNA/reports/deduplicated_sizes_stats__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_smallrna("{sample_name}", "deduplicate_srna_nextflexv3", "all"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["deduplicate_srna_nextflexv3"]["threads"]
    resources:
        mem_mb=config["resources"]["deduplicate_srna_nextflexv3"]["mem_mb"],
        tmp_mb=config["resources"]["deduplicate_srna_nextflexv3"]["tmp_mb"],
        qos=config["resources"]["deduplicate_srna_nextflexv3"]["qos"]
    shell:
        """
        {{
        #### 1) Collapse the PCR-duplicated reads
        seqcluster collapse -f {input.fastq} -o {output.collapse_folder}

        ### 2) Trimming the read-specific UMIs (first and last 4bp, and adding a minimum 15bp)
        seqtk trimfq -b 4 -e 4 {output.collapsed_fastq} | seqtk seq -L 15 > {output.deduplicated_fastq}
        cat {output.deduplicated_fastq} | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"deduplicated",$2,$1}}' >> "{output.report}"
        
        pigz -p {threads} -c {output.deduplicated_fastq} > {output.gzipped_fastq}
        }} 2>&1 | tee -a "{log}"
        """

rule make_bt2_indices_for_structural_RNAs:
    input:
        fasta = lambda wildcards: config[wildcards.ref_genome]['structural_rna_fafile']
    output:
        temp_fasta = temp("genomes/structural_RNAs/{ref_genome}/temp.fa"),
        indices = directory("genomes/structural_RNAs/{ref_genome}_bt2_index")
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","structural_RNA_bt2_index_{ref_genome}.log"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["make_bt2_indices_for_structural_RNAs"]["threads"]
    resources:
        mem_mb=config["resources"]["make_bt2_indices_for_structural_RNAs"]["mem_mb"],
        tmp_mb=config["resources"]["make_bt2_indices_for_structural_RNAs"]["tmp_mb"],
        qos=config["resources"]["make_bt2_indices_for_structural_RNAs"]["qos"]
    shell:
        """
        {{
        if [[ {input.fasta} == *.gz ]]; then
            printf "Fasta of structural RNAs for {wildcards.ref_genome} is gzipped file: {input.fasta}\n"
            gunzip {input.fasta} -c > {output.temp_fasta}
        else
            printf "Fasta of structural RNAs for {wildcards.ref_genome} is unzipped file: {input.fasta}\n"
            cp {input.fasta} {output.temp_fasta}
        fi
        printf "\nBuilding Bowtie2 index for {wildcards.ref_genome}\n"
        mkdir -p genomes/structural_RNAs/{wildcards.ref_genome}_bt2_index
        bowtie2-build --threads {threads} "{output.temp_fasta}" "{output.indices}/{wildcards.ref_genome}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_structural_rna:
    input:
        fastq = lambda wildcards: f"results/sRNA/fastq/{define_input_file_for_structural(wildcards.sample_name)}.fastq.gz",
        indices = lambda wildcards: f"genomes/structural_RNAs/{parse_sample_name(wildcards.sample_name)['ref_genome']}_bt2_index"
    output:
        filtered_fastq = temp("results/sRNA/fastq/filtered__{sample_name}__R0.fastq"),
        gzipped_fastq = temp("results/sRNA/fastq/filtered__{sample_name}__R0.fastq.gz"),
        report = "results/sRNA/reports/filtered_sizes_stats__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_smallrna("{sample_name}", "filter_structural_rna", "all"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["filter_structural_rna"]["threads"]
    resources:
        mem_mb=config["resources"]["filter_structural_rna"]["mem_mb"],
        tmp_mb=config["resources"]["filter_structural_rna"]["tmp_mb"],
        qos=config["resources"]["filter_structural_rna"]["qos"]
    shell:
        """
        {{
        bowtie2 --very-sensitive -p {threads} -x "{input.indices}/{params.ref_genome}" -U {input.fastq} | samtools view -@ {threads} -f 0x4 | samtools fastq -@ {threads} > {output.filtered_fastq}
        cat {output.filtered_fastq} | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"filtered",$2,$1}}' > "{output.report}"
        pigz -p {threads} {output.filtered_fastq} -c > {output.gzipped_fastq}
        }} 2>&1 | tee -a "{log}"
        """

rule dispatch_srna_fastq:
    input:
        fastq = lambda wildcards: f"results/sRNA/fastq/{define_input_file_for_shortstack(wildcards.sample_name)}.fastq.gz"
    output:
        fastq_file = "results/sRNA/fastq/clean__{sample_name}.fastq.gz"
    conda: CONDA_ENV_SRNA
    localrule: True
    shell:
        """
        cp {input.fastq} {output.fastq_file}
        """

rule make_bowtie1_indices:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    output:
        indices = multiext("genomes/{ref_genome}/{ref_genome}.fa", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","bt1_index_{ref_genome}.log"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["make_bowtie1_indices"]["threads"]
    resources:
        mem_mb=config["resources"]["make_bowtie1_indices"]["mem_mb"],
        tmp_mb=config["resources"]["make_bowtie1_indices"]["tmp_mb"],
        qos=config["resources"]["make_bowtie1_indices"]["qos"]
    shell:
        """
        {{
        printf "\nMaking Bowtie1 indices for {wildcards.ref_genome}\n"
        bowtie-build {input.fasta} {input.fasta}
        }} 2>&1 | tee -a "{log}"
        """

rule make_bowtie1_indices_large:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    output:
        indices = multiext("genomes/{ref_genome}/{ref_genome}.fa", ".1.ebwtl", ".2.ebwtl", ".3.ebwtl", ".4.ebwtl", ".rev.1.ebwtl", ".rev.2.ebwtl")
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","large_bt1_index_{ref_genome}.log"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["make_bowtie1_indices"]["threads"]
    resources:
        mem_mb=config["resources"]["make_bowtie1_indices_large"]["mem_mb"],
        tmp_mb=config["resources"]["make_bowtie1_indices_large"]["tmp_mb"],
        qos=config["resources"]["make_bowtie1_indices_large"]["qos"]
    shell:
        """
        {{
        printf "\nMaking large Bowtie1 indices for {wildcards.ref_genome}\n"
        bowtie-build {input.fasta} {input.fasta}
        }} 2>&1 | tee -a "{log}"
        """
        
rule shortstack_map:
    input:
        fastq = "results/sRNA/fastq/clean__{sample_name}.fastq.gz",
        fasta = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa",
        indices = get_bt1_indices
    output:
        count_file = "results/sRNA/mapped/{sample_name}/Results.txt",
        bam_file = "results/sRNA/mapped/{sample_name}/clean__{sample_name}_condensed.bam",
        bai_file = "results/sRNA/mapped/{sample_name}/clean__{sample_name}_condensed.bam.csi",
        touch_file = "results/sRNA/chkpts/map_sRNA__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        srna_params = config['srna_mapping_params']
    log:
        temp(return_log_smallrna("{sample_name}", "mapping_shortstack", "all"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["shortstack_map"]["threads"]
    resources:
        mem_mb=config["resources"]["shortstack_map"]["mem_mb"],
        tmp_mb=config["resources"]["shortstack_map"]["tmp_mb"],
        qos=config["resources"]["shortstack_map"]["qos"]
    shell:
        """
        {{
        rm -rf results/sRNA/mapped/{params.sample_name}
        printf "\nMapping {params.sample_name} to {params.ref_genome} with Shortstack version:\n"
        ShortStack --version
        ShortStack --readfile {input.fastq} --genomefile {input.fasta} --threads {threads} {params.srna_params} --outdir results/sRNA/mapped/{params.sample_name}        
        touch {output.touch_file}
        }} 2>&1 | tee -a "{log}"
        """

rule make_cluster_bedfiles:
    input:
        count_file = "results/sRNA/mapped/{sample_name}/Results.txt"
    output:
        cluster_bedfile = "results/sRNA/mapped/{sample_name}/clusters.bed"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        srna_min = config['srna_min_size'],
        srna_max = config['srna_max_size']
    log:
        temp(return_log_smallrna("{sample_name}", "make_cluster_bedfiles", "all"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["make_cluster_bedfiles"]["threads"]
    resources:
        mem_mb=config["resources"]["make_cluster_bedfiles"]["mem_mb"],
        tmp_mb=config["resources"]["make_cluster_bedfiles"]["tmp_mb"],
        qos=config["resources"]["make_cluster_bedfiles"]["qos"]
    shell:
        """
        {{
        ## To create a bedfile of clusters for Upset plots
        awk -v OFS="\t" -v m={params.srna_min} -v n={params.srna_max} 'NR==1 {{for (i=1; i<=NF; i++) {{if ($i=="DicerCall") dicer_col=i; if ($i=="MIRNA") mirna_col=i}}}} NR>1 {{if ($mirna_col=="Y") t="MIRNA"; else if ($dicer_col>=m && $dicer_col<=n) t=$dicer_col"nt"; else t="Others"; print $3, $4-1, $5, t}}' {input.count_file} > {output.cluster_bedfile}
        }} 2>&1 | tee -a "{log}"
        """
        
rule make_srna_size_stats:
    input:
        bamfile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}_condensed.bam",
        baifile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}_condensed.bam.csi"
    output:
        report = "results/sRNA/reports/sizes_stats__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        temp(return_log_smallrna("{sample_name}", "make_srna_stats", "all"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["make_srna_size_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["make_srna_size_stats"]["mem_mb"],
        tmp_mb=config["resources"]["make_srna_size_stats"]["tmp_mb"],
        qos=config["resources"]["make_srna_size_stats"]["qos"]
    shell:
        """
        {{
        printf "\nGetting stats for {params.sample_name}\n"
        printf "Sample\tType\tSize\tCount\n" > {output.report}
        printf "\nGetting filtered stats for {params.sample_name}\n"
        if [[ -s results/sRNA/reports/deduplicated_sizes_stats__{params.sample_name}.txt ]]; then
            cat results/sRNA/reports/deduplicated_sizes_stats__{params.sample_name}.txt >> "{output.report}"
            zcat results/sRNA/fastq/trim__{params.sample_name}__R0.fastq.gz | awk '{{if(NR%4==2) {{a=length($1)-8; if (a>=15) print a}}}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"trimmed",$2,$1}}' >> "{output.report}"
        else
            zcat results/sRNA/fastq/trim__{params.sample_name}__R0.fastq.gz | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c | awk -v OFS="\t" -v n={params.sample_name} '{{print n,"trimmed",$2,$1}}' >> "{output.report}"
        fi
        if [[ -s results/sRNA/reports/filtered_sizes_stats__{params.sample_name}.txt ]]; then
            cat results/sRNA/reports/filtered_sizes_stats__{params.sample_name}.txt >> "{output.report}"
        fi
        samtools view {input.bamfile} | awk '$2==0 || $2==16 {{print length($10)"_"$1}}' | sort -u | awk -F'_' '{{print $1,$NF}}' | sort -n | awk -v OFS="\t" -v n={params.sample_name} '{{sum[$1]+=$2}} END {{for (i in sum) print n,"mapped",i,sum[i]}}' >> "{output.report}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_size_srna_sample:
    input:
        bamfile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}_condensed.bam",
        baifile = "results/sRNA/mapped/{sample_name}/clean__{sample_name}_condensed.bam.csi"
    output:
        filtered_file = temp("results/sRNA/mapped/sized__{size}nt__{sample_name}.bam"),
        index_file = temp("results/sRNA/mapped/sized__{size}nt__{sample_name}.bam.bai")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        size = lambda wildcards: wildcards.size
    log:
        temp(return_log_smallrna("{sample_name}", "filter_size_srna", "{size}"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["filter_size_srna_sample"]["threads"]
    resources:
        mem_mb=config["resources"]["filter_size_srna_sample"]["mem_mb"],
        tmp_mb=config["resources"]["filter_size_srna_sample"]["tmp_mb"],
        qos=config["resources"]["filter_size_srna_sample"]["qos"]
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
        mergefile = temp("results/sRNA/mapped/merged__{size}nt__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"),
        indexfile = temp("results/sRNA/mapped/merged__{size}nt__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam.bai")
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis'),
        size = lambda wildcards: wildcards.size
    log:
        temp(return_log_smallrna("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging_srna_reps", "{size}"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["merging_srna_replicates"]["threads"]
    resources:
        mem_mb=config["resources"]["merging_srna_replicates"]["mem_mb"],
        tmp_mb=config["resources"]["merging_srna_replicates"]["tmp_mb"],
        qos=config["resources"]["merging_srna_replicates"]["qos"]
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
        bamfile = lambda wildcards: f"results/sRNA/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'sized'}__{wildcards.size}nt__{wildcards.sample_name}.bam",
        baifile = lambda wildcards: f"results/sRNA/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'sized'}__{wildcards.size}nt__{wildcards.sample_name}.bam.bai",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        temp_folder = temp(directory("results/sRNA/tracks/{sample_name}__{size}")),
        temp_minus = temp("results/sRNA/tracks/{sample_name}__{size}/{sample_name}__{size}nt__minus.bg"),
        temp_minus_rev = temp("results/sRNA/tracks/{sample_name}__{size}/{sample_name}__{size}nt__minus_rev.bg"),
        temp_minus_sort = temp("results/sRNA/tracks/{sample_name}__{size}/{sample_name}__{size}nt__minus_sort.bg"),
        bw_plus = "results/sRNA/tracks/{sample_name}__{size}nt__plus.bw",
        bw_minus = "results/sRNA/tracks/{sample_name}__{size}nt__minus.bw"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        size = lambda wildcards: wildcards.size,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_smallrna("{sample_name}", "making_bigiwig", "{size}"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["make_srna_stranded_bigwigs"]["threads"]
    resources:
        mem_mb=config["resources"]["make_srna_stranded_bigwigs"]["mem_mb"],
        tmp_mb=config["resources"]["make_srna_stranded_bigwigs"]["tmp_mb"],
        qos=config["resources"]["make_srna_stranded_bigwigs"]["qos"]
    shell:
        """
        {{
        printf "Getting stranded coverage for {params.sample_name} {params.size}nt\n"
        input_bamfile="{input.bamfile}"
        basename=${{input_bamfile%.bam}}
        ShortTracks --mode simple --stranded --bamfile {input.bamfile}
        mv ${{basename}}_p.bw {output.bw_plus}
        printf "Inverting minus strand (back to positive values)\n"
        bigWigToBedGraph ${{basename}}_m.bw {output.temp_minus}
        awk -v OFS="\t" '{{print $1,$2,$3,-$4}}' {output.temp_minus} > {output.temp_minus_rev}
        bedSort {output.temp_minus_rev} {output.temp_minus_sort}
        bedGraphToBigWig {output.temp_minus_sort} {input.chrom_sizes} {output.bw_minus}
        rm -f ${{basename}}_m*
        rm -f ${{basename}}_p*
        }} 2>&1 | tee -a "{log}"
        """

rule analyze_all_srna_samples_on_target_file:
    input:
        bamfiles = lambda wildcards: define_input_for_grouped_analysis(wildcards.ref_genome),
        fasta = lambda wildcards: f"genomes/{wildcards.ref_genome}/{wildcards.ref_genome}.fa",
        target_file = lambda wildcards: define_srna_target_file(wildcards)
    output:
        count_file = "results/sRNA/clusters/{analysis_name}__{ref_genome}__on_{target_name}/Counts.txt"
    wildcard_constraints:
        ref_genome = r"[^/]+"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        srna_min = config['srna_min_size'],
        srna_max = config['srna_max_size']
    log:
        temp(return_log_smallrna("{ref_genome}", "{analysis_name}_analysis", "{target_name}"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["analyze_all_srna_samples_on_target_file"]["threads"]
    resources:
        mem_mb=config["resources"]["analyze_all_srna_samples_on_target_file"]["mem_mb"],
        tmp_mb=config["resources"]["analyze_all_srna_samples_on_target_file"]["tmp_mb"],
        qos=config["resources"]["analyze_all_srna_samples_on_target_file"]["qos"]
    shell:
        """
        {{
        rm -rf results/sRNA/clusters/{params.analysis_name}__{params.ref_genome}__on_{params.target_name}
        if [[ "{params.target_name}" == "new_clusters" ]]; then
            printf "\nAnalyszing all samples from {params.analysis_name} on {params.ref_genome} with Shortstack version:\n"
            ShortStack --version
            ShortStack --bamfile {input.bamfiles} --genomefile {input.fasta} --threads {threads} --dicermin {params.srna_min} --dicermax {params.srna_max} --outdir results/sRNA/clusters/{params.analysis_name}__{params.ref_genome}__on_{params.target_name}
        else
            printf "\nAnalyszing all samples from {params.analysis_name} on {params.ref_genome} limited to {params.target_name} with Shortstack version:\n"
            ShortStack --version
            ShortStack --bamfile {input.bamfiles} --genomefile {input.fasta} --threads {threads} --dicermin {params.srna_min} --dicermax {params.srna_max} --locifile {input.target_file} --outdir results/sRNA/clusters/{params.analysis_name}__{params.ref_genome}__on_{params.target_name}
        fi
        if [[ ! -e {output.count_file} ]]; then
            touch {output.count_file}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule prep_files_for_differential_srna_clusters:
    input: 
        count_file = "results/sRNA/clusters/{analysis_name}__{ref_genome}__on_{target_name}/Counts.txt"        
    output:
        srna_samples = "results/sRNA/clusters/{analysis_name}__{ref_genome}__on_{target_name}/samples_for_edgeR.txt",
        srna_counts = "results/sRNA/clusters/{analysis_name}__{ref_genome}__on_{target_name}/counts_for_edgeR.txt"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name
    log:
        temp(return_log_smallrna("{ref_genome}", "{analysis_name}_prep", "{target_name}"))
    threads: config["resources"]["prep_files_for_differential_srna_clusters"]["threads"]
    resources:
        mem_mb=config["resources"]["prep_files_for_differential_srna_clusters"]["mem_mb"],
        tmp_mb=config["resources"]["prep_files_for_differential_srna_clusters"]["tmp_mb"],
        qos=config["resources"]["prep_files_for_differential_srna_clusters"]["qos"]
    run:
        filtered_samples = samples[ (samples['data_type'] == 'sRNA') & (samples['ref_genome'] == params.ref_genome) ].copy()
        filtered_samples['Sample'] = filtered_samples['line'] + "__" + filtered_samples['tissue']
        filtered_samples['Replicate'] = filtered_samples['Sample'] + "__" + filtered_samples['replicate'].astype(str)
        
        sRNA_samples = filtered_samples[['Replicate','Sample']].drop_duplicates()    
        sRNA_samples = sRNA_samples.sort_values(by=['Sample', 'Replicate'],ascending=[True, True]).reset_index(drop=True)
        sRNA_samples['Color'] = pd.factorize(sRNA_samples['Sample'])[0] + 1

        sRNA_samples.to_csv(output.srna_samples, sep="\t", index=False)
        
        column_order = ['Name']
        for _, row in sRNA_samples.iterrows():
            ROW = filtered_samples.loc[filtered_samples["Replicate"] == row["Replicate"]].iloc[0]
            sname = sample_name_str(ROW, 'sample')
            column_order.append(sname)
            
        temp = pd.read_csv(input.count_file, sep="\t", header=0)
        temp = temp.rename(columns=lambda x: x[7:] if x.startswith("clean__") else x)
        temp = temp.rename(columns=lambda x: x[:-10] if x.endswith("_condensed") else x)
        sRNA_counts = temp[column_order]
        sRNA_counts.to_csv(output.srna_counts, sep="\t", index=False)

rule call_all_differential_srna_clusters:
    input:
        srna_samples = "results/sRNA/clusters/{analysis_name}__{ref_genome}__on_{target_name}/samples_for_edgeR.txt",
        srna_counts = "results/sRNA/clusters/{analysis_name}__{ref_genome}__on_{target_name}/counts_for_edgeR.txt"
    output:
        mdsplot = "results/combined/plots/MDS_sRNA_{analysis_name}_{ref_genome}__on_{target_name}_d12.pdf",
        touch = "results/sRNA/chkpts/calling_differential_sRNA_clusters__{analysis_name}__{ref_genome}__on_{target_name}.done"
    params:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_call_srna_DEGs.R"),
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        region_file = lambda wildcards: f"results/sRNA/clusters/{analysis_name}__{wildcards.ref_genome}__on_{wildcards.target_name}/Results.txt" if define_srna_target_file(wildcards) == [] else define_srna_target_file(wildcards)
    log:
        temp(return_log_smallrna("{ref_genome}", "{analysis_name}_cluster_degs", "{target_name}"))
    conda: CONDA_ENV_SRNA
    threads: config["resources"]["call_all_DEGs"]["threads"]
    resources:
        mem_mb=config["resources"]["call_all_DEGs"]["mem_mb"],
        tmp_mb=config["resources"]["call_all_DEGs"]["tmp_mb"],
        qos=config["resources"]["call_all_DEGs"]["qos"]
    shell:
        """
        printf "running edgeR for all samples in {params.ref_genome}\n"
        Rscript "{params.script}" "{input.srna_counts}" "{input.srna_samples}" "{params.analysis_name}" "{params.ref_genome}" "{params.target_name}" "{params.region_file}"
        touch {output.touch}
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
