# function to access logs more easily
def return_log_chip(env, sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"results",env,"logs",f"tmp__{sample_name}__{step}__{paired}.log")

def assign_mapping_paired(wildcards, rulename, outputfile):
    sname = wildcards.sample_name
    env = get_sample_info_from_name(sname, samples, 'env')
    paired = get_sample_info_from_name(sname, samples, 'paired')
    if paired == "PE":
        rule_obj = getattr(rules, f"{rulename}_pe")
    elif paired == "SE":
        rule_obj = getattr(rules, f"{rulename}_se")
        
    return getattr(rule_obj.output, outputfile).format(sample_name=sname, env=env)

def assign_chip_input(wildcards):
    inputname = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__Input__{wildcards.replicate}__{wildcards.ref_genome}"
    if wildcards.file_type in ['merged', 'pseudo1', 'pseudo2']:
        return inputname
    elif inputname in samples['sample_name']:
        return inputname
    else:
        ipname = sample_name_str(wildcards, 'sample')
        ippaired = get_sample_info_from_name(ipname, samples, 'paired')
        alts = []
        for rep in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, "Input", wildcards.ref_genome), []):
            alt_inputname = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__Input__{rep}__{wildcards.ref_genome}"
            alts.append(f"{alt_inputname}")
            if get_sample_info_from_name(alt_inputname, samples, 'paired') == ippaired:
                return alt_inputname
            else:
                alts.append(f"{ippaired}")
        
        raise ValueError(f"\nSample '{ipname}' does not have corresponding Input among:\n{alts}")
                
def get_peaktype(sample_type, peaktype_config):
    if sample_type == "IP":
        return "narrow"
    for pattern, peaktype in peaktype_config.items():
        if re.search(pattern, sample_type):
            return peaktype
    raise ValueError(f"\nNo peaktype found for sample_type '{sample_type}")

def assign_peak_files_for_idr(wildcards):
    sname = sample_name_str(wildcards, 'analysis')
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    peaktype = get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"])
    replicates = analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), [])
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    if paired == "PE":
        return [ f"results/{env}/peaks/peaks_pe__final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}_peaks.{peaktype}Peak"
                for replicate in replicates ]
    else:
        return [ f"results/{env}/peaks/peaks_se__final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}_peaks.{peaktype}Peak"
                for replicate in replicates ]

def input_peak_files_for_best_peaks(wildcards):
    peaktype = get_peaktype(wildcards.sample_type, config["chip_callpeaks"]['peaktype'])
    sname = sample_name_str(wildcards, 'analysis')
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')

    if len(analysis_to_replicates[(wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome)]) >= 2:
        if paired == "PE":
            result = [ f"results/{env}/peaks/peaks_pe__merged__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__merged__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_pe__pseudo1__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__merged__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_pe__pseudo2__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__merged__{wildcards.ref_genome}_peaks.{peaktype}Peak" ]
        else:
            result = [ f"results/{env}/peaks/peaks_pe__merged__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__merged__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_pe__pseudo1__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__merged__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_pe__pseudo2__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__merged__{wildcards.ref_genome}_peaks.{peaktype}Peak" ]
    else:
        one_rep = analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), [])[0]
        if paired == "PE":
            result = [ f"results/{env}/peaks/peaks_se__final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{one_rep}__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_se__pseudo1__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{one_rep}__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_se__pseudo2__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{one_rep}__{wildcards.ref_genome}_peaks.{peaktype}Peak" ]
        else:
            result = [ f"results/{env}/peaks/peaks_se__final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{one_rep}__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_se__pseudo1__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{one_rep}__{wildcards.ref_genome}_peaks.{peaktype}Peak",
                       f"results/{env}/peaks/peaks_se__pseudo2__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{one_rep}__{wildcards.ref_genome}_peaks.{peaktype}Peak" ]

    return result

def get_replicate_name(wildcards, pos):
    sname = wildcards.sample_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    prefix = "peaks_pe" if paired == "PE" else "peaks_se"
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    data_type = get_sample_info_from_name(sname, analysis_samples, 'data_type')
    line = get_sample_info_from_name(sname, analysis_samples, 'line')
    tissue = get_sample_info_from_name(sname, analysis_samples, 'tissue')
    sample_type = get_sample_info_from_name(sname, analysis_samples, 'sample_type')
    peaktype = get_peaktype(sample_type, config["chip_callpeaks"]["peaktype"])
    ref_genome = get_sample_info_from_name(sname, analysis_samples, 'ref_genome')
    rep_list = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])
    
    if pos >= len(rep_list): 
        return "missingrep"
    else:
        return f"results/{env}/peaks/{prefix}__final__{data_type}__{line}__{tissue}__{sample_type}__{rep_list[pos]}__{ref_genome}_peaks.{peaktype}Peak"

def get_replicate_pairs(wildcards):
    sname = sample_name_str(wildcards, 'analysis')
    reps = analysis_to_replicates.get((wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), [])
    pairs = []
    for i in range(0,len(reps)):
        for j in range(i+1, len(reps)):
            rep_i = reps[i]
            rep_j = reps[j]
            pairs.append(f"{rep_i} {rep_j}")
    return pairs

def define_logs_final_input(wildcards):
    log_files = []
    sname = wildcards.sample_name
    data_type = get_sample_info_from_name(sname, analysis_samples, 'data_type')
    line = get_sample_info_from_name(sname, analysis_samples, 'line')
    tissue = get_sample_info_from_name(sname, analysis_samples, 'tissue')
    sample_type = get_sample_info_from_name(sname, analysis_samples, 'sample_type')
    ref_genome = get_sample_info_from_name(sname, analysis_samples, 'ref_genome')
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    peaktype = get_peaktype(sample_type, config["chip_callpeaks"]['peaktype'])
    for rep in analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), []):
        namerep = f"{data_type}__{line}__{tissue}__{sample_type}__{rep}__{ref_genome}"
        log_files.append(return_log_chip(env, namerep, f"final__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, namerep, f"making_bigwig_final", ""))
        log_files.append(return_log_chip(env, namerep, f"making_fingerprint_final", ""))
    
    if len(analysis_to_replicates[(data_type, line, tissue, sample_type, ref_genome)]) >= 2:
        log_files.append(return_log_chip(env, sname, "IDR", ""))
        log_files.append(return_log_chip(env, sname, "merging_reps", ""))
        
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}", f"merged__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}", f"pseudo1__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}", f"pseudo2__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}", "splitting_pseudreps", ""))
    else:
        one_rep = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])[0]
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__{one_rep}__{ref_genome}", "splitting_pseudreps", ""))
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__{one_rep}__{ref_genome}", f"pseudo1__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__{one_rep}__{ref_genome}", f"pseudo2__{peaktype}peak_calling", paired))
        
    log_files.append(return_log_chip(env, f"{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "selecting_best_peaks", ""))
    
    return log_files

def define_final_chip_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    map_files = []
    stat_files = []
    qc_files = []
    peak_files = []
    bigwig_files = []
    filtered_rep_samples = samples[ ((samples['env'] == 'ChIP') | (samples['env'] == 'TF')) & (samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        env = get_sample_info_from_name(sname, samples, 'env')
        if paired == "PE":
            map_files.append(f"results/{env}/logs/process_chip_pe_sample__{sname}.log") # mapping stats for each paired-end replicate
            qc_files.append(f"results/{env}/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
            qc_files.append(f"results/{env}/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            qc_files.append(f"results/{env}/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
            qc_files.append(f"results/{env}/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
        else:
            map_files.append(f"results/{env}/logs/process_chip_pse_sample__{sname}.log") # mapping stats for each single-end replicate
            qc_files.append(f"results/{env}/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file
            qc_files.append(f"results/{env}/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
            
    filtered_rep_samples_no_input = filtered_rep_samples[ (filtered_rep_samples['sample_type'] != "Input") ].copy()
    for _, row in filtered_rep_samples_no_input.iterrows():
        peaktype = get_peaktype(row.sample_type, config["chip_callpeaks"]['peaktype'])
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        env = get_sample_info_from_name(sname, samples, 'env')
        bigwig_files.append(f"results/{env}/tracks/FC__final__{sname}.bw") # bigwig log2FC enrichment vs input for each replicate
        stat_files.append(f"results/{env}/plots/Fingerprint__final__{sname}.png") # fingerprint plots for each replicate and its input
        if paired == "PE":
            peak_files.append(f"results/{env}/peaks/peaks_pe__final__{sname}_peaks.{peaktype}Peak") # peak file for each paired-end replicate
        else:
            peak_files.append(f"results/{env}/peaks/peaks_se__final__{sname}_peaks.{peaktype}Peak") # peak file for each single-end replicate
            
    filtered_analysis_samples = analysis_samples[ ((analysis_samples['env'] == 'ChIP') | (analysis_samples['env'] == 'TF')) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        env = get_sample_info_from_name(spname, analysis_samples, 'env')
        peak_files.append(f"results/{env}/peaks/selected_peaks__{spname}.bed") # best peak file for each analysis sample
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            stat_files.append(f"results/{env}/chkpts/idr__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{row.ref_genome}.done") # idr analyses between each pair of replicates
            bigwig_files.append(f"results/{env}/tracks/FC__merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.bw") # bigiwig log2FC for merged replicates vs merged inputs
        
    if qc_option == "all":
        results = map_files + qc_files
    else:
        results = map_files 
        
    if analysis:
        results += bigwig_files + peak_files + stat_files
    
    return results
        
rule make_bt2_indices:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gff = "genomes/{ref_genome}/{ref_genome}.gff",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        indices = directory("genomes/{ref_genome}/bt2_index")
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","bowtie_index_{ref_genome}.log"))
    conda: CONDA_ENV
    threads: config["resources"]["make_bt2_indices"]["threads"]
    resources:
        mem=config["resources"]["make_bt2_indices"]["mem"],
        tmp=config["resources"]["make_bt2_indices"]["tmp"]
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
        fastq1 = "results/{env}/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "results/{env}/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/bt2_index"
    output:
        samfile = temp("results/{env}/mapped/mapped_pe__{sample_name}.sam"),
        metrics = "results/{env}/reports/bt2_pe__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: config['chip_mapping_option'],
        mapping_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['map_pe']    
    log:
        temp(return_log_chip("{env}","{sample_name}", "mappingBT2", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["bowtie2_map_pe"]["threads"]
    resources:
        mem=config["resources"]["bowtie2_map_pe"]["mem"],
        tmp=config["resources"]["bowtie2_map_pe"]["tmp"]
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
        fastq = "results/{env}/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/bt2_index"
    output:
        samfile = temp("results/{env}/mapped/mapped_se__{sample_name}.sam"),
        metrics = "results/{env}/reports/bt2_se__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: config['chip_mapping_option'],
        mapping_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['map_se']    
    log:
        temp(return_log_chip("{env}","{sample_name}", "mappingBT2", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["bowtie2_map_se"]["threads"]
    resources:
        mem=config["resources"]["bowtie2_map_se"]["mem"],
        tmp=config["resources"]["bowtie2_map_se"]["tmp"]
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
        samfile = "results/{env}/mapped/mapped_pe__{sample_name}.sam"
    output:
        bamfile = temp("results/{env}/mapped/mapped_pe__{sample_name}.bam"),
        metrics_dup = "results/{env}/reports/markdup_pe__{sample_name}.txt",
        metrics_flag = "results/{env}/reports/flagstat_pe__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        map_option = lambda wildcards: config['chip_mapping_option'],
        filtering_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['filter']    
    log:
        temp(return_log_chip("{env}","{sample_name}", "filteringChIP", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["filter_chip_pe"]["threads"]
    resources:
        mem=config["resources"]["filter_chip_pe"]["mem"],
        tmp=config["resources"]["filter_chip_pe"]["tmp"]
    shell:
        """
        {{
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {params.sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o "results/{params.env}/mapped/temp1_{params.sample_name}.bam" "{input.samfile}"
        samtools fixmate -@ {threads} -m "results/{params.env}/mapped/temp1_{params.sample_name}.bam" "results/{params.env}/mapped/temp2_{params.sample_name}.bam"
        samtools sort -@ {threads} -o "results/{params.env}/mapped/temp3_{params.sample_name}.bam" "results/{params.env}/mapped/temp2_{params.sample_name}.bam"
        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} "results/{params.env}/mapped/temp3_{params.sample_name}.bam" "{output.bamfile}"
        samtools index -@ {threads} "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f results/{params.env}/mapped/temp*"_{params.sample_name}.bam"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_chip_se:
    input:
        samfile = "results/{env}/mapped/mapped_se__{sample_name}.sam"
    output:
        bamfile = temp("results/{env}/mapped/mapped_se__{sample_name}.bam"),
        metrics_dup = "results/{env}/reports/markdup_se__{sample_name}.txt",
        metrics_flag = "results/{env}/reports/flagstat_se__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        map_option = lambda wildcards: config['chip_mapping_option'],
        filtering_params = lambda wildcards: config['chip_mapping'][config['chip_mapping_option']]['filter']    
    log:
        temp(return_log_chip("{env}","{sample_name}", "filteringChIP", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["filter_chip_se"]["threads"]
    resources:
        mem=config["resources"]["filter_chip_se"]["mem"],
        tmp=config["resources"]["filter_chip_se"]["tmp"]
    shell:
        """
        {{
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {params.sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o "results/{params.env}/mapped/temp1_{params.sample_name}.bam" "{input.samfile}"
        samtools sort -@ {threads} -o "results/{params.env}/mapped/temp2_{params.sample_name}.bam" "results/{params.env}/mapped/temp1_{params.sample_name}.bam"
        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} "results/{params.env}/mapped/temp2_{params.sample_name}.bam" "{output.bamfile}"
        samtools index -@ {threads} "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f results/{params.env}/mapped/temp*"_{params.sample_name}.bam"
        }} 2>&1 | tee -a "{log}"
        """

rule make_chip_stats_pe:
    input:
        metrics_trim = "results/{env}/reports/trim_pe__{sample_name}.txt",
        metrics_map = "results/{env}/reports/bt2_pe__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.env, wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingBT2", "filteringChIP"] ]
    output:
        stat_file = "results/{env}/reports/summary_ChIP_PE_mapping_stats_{sample_name}.txt",
        log = "results/{env}/logs/process_chip_pe_sample__{sample_name}.log"        
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        line = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'line'),
        tissue = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'tissue'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'sample_type'),
        replicate = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'replicate'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'ref_genome')
    threads: config["resources"]["make_chip_stats_pe"]["threads"]
    resources:
        mem=config["resources"]["make_chip_stats_pe"]["mem"],
        tmp=config["resources"]["make_chip_stats_pe"]["tmp"]
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned concordantly >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned concordantly exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule make_chip_stats_se:
    input:
        metrics_trim = "results/{env}/reports/trim_se__{sample_name}.txt",
        metrics_map = "results/{env}/reports/bt2_se__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.env, wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingBT2", "filteringChIP"] ]
    output:
        stat_file = "results/{env}/reports/summary_ChIP_SE_mapping_stats_{sample_name}.txt",
        log = "results/{env}/logs/process_chip_se_sample__{sample_name}.log"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        line = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'line'),
        tissue = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'tissue'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'sample_type'),
        replicate = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'replicate'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'ref_genome')
    threads: config["resources"]["make_chip_stats_se"]["threads"]
    resources:
        mem=config["resources"]["make_chip_stats_se"]["mem"],
        tmp=config["resources"]["make_chip_stats_se"]["tmp"]
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total reads processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        """

rule pe_or_se_chip_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "filter_chip", "bamfile")
    output:
        bam = "results/{env}/mapped/final__{sample_name}.bam",
        touch = temp("results/{env}/chkpts/map_chip__{sample_name}.done")
    wildcard_constraints:
        env = "ChIP|TF"
    localrule: True
    shell:
        """
        mv {input} {output.bam}
        mv {input}.bai {output.bam}.bai
        touch {output.touch} 
        """
    
rule make_coverage_chip:
    input: 
        bamfile = "results/{env}/mapped/final__{sample_name}.bam"
    output:
        bigwigcov = "results/{env}/tracks/coverage__{sample_name}.bw"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        binsize = config['chip_tracks']['binsize']
    conda: CONDA_ENV
    threads: config["resources"]["make_coverage_chip"]["threads"]
    resources:
        mem=config["resources"]["make_coverage_chip"]["mem"],
        tmp=config["resources"]["make_coverage_chip"]["tmp"]
    shell:
        """
        bamCoverage -b {input.bamfile} -o {output.bigwigcov} -bs {params.binsize} -p {threads}
        """

rule make_bigwig_chip:
    input: 
        ipfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        bigwigfile = "results/{env}/tracks/FC__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bw"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        ipname = lambda wildcards: f"{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{wildcards.file_type}__{assign_chip_input(wildcards)}",
        binsize = config['chip_tracks']['binsize'],
        params = config['chip_tracks']['params']
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "making_bigwig_{file_type}", ""))
    conda: CONDA_ENV
    threads: config["resources"]["make_bigwig_chip"]["threads"]
    resources:
        mem=config["resources"]["make_bigwig_chip"]["mem"],
        tmp=config["resources"]["make_bigwig_chip"]["tmp"]
    shell:
        """
        {{
        printf "\nMaking bigwig files for {params.ipname} (vs {params.inputname}) with deeptools version:\n"
        deeptools --version
        bamCompare -b1 {input.ipfile} -b2 {input.inputfile} -o {output.bigwigfile} -p {threads} --binSize {params.binsize} {params.params}
        }} 2>&1 | tee -a "{log}"
        """

rule make_fingerprint_plot:
    input: 
        ipfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        pngplot = "results/{env}/plots/Fingerprint__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.png"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        ipname = lambda wildcards: f"{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{wildcards.file_type}__{assign_chip_input(wildcards)}"
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "making_fingerprint_{file_type}", ""))
    conda: CONDA_ENV
    threads: config["resources"]["make_fingerprint_plot"]["threads"]
    resources:
        mem=config["resources"]["make_fingerprint_plot"]["mem"],
        tmp=config["resources"]["make_fingerprint_plot"]["tmp"]
    shell:
        """
        {{
        printf "\nPlotting fingerprint for {params.ipname} (vs {params.inputname}) with deeptools version:\n"
        deeptools --version
        plotFingerprint -b {input.ipfile} {input.inputfile} -o {output.pngplot} -p {threads} -l {params.ipname} {params.inputname}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_macs2_pe:
    input:
        ipfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        peakfile = "results/{env}/peaks/peaks_pe__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}_peaks.{peaktype}Peak"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        ipname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{assign_chip_input(wildcards)}",
        peaktype = lambda wildcards: get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"]),
        filetype = lambda wildcards: {wildcards.file_type},
        env = lambda wildcards: {wildcards.env},
        params = config["chip_callpeaks"]['params'],
        genomesize = config[config['species']]['genomesize']
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "{file_type}__{peaktype}peak_calling", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["calling_peaks_macs2_pe"]["threads"]
    resources:
        mem=config["resources"]["calling_peaks_macs2_pe"]["mem"],
        tmp=config["resources"]["calling_peaks_macs2_pe"]["tmp"]
    shell:
        """
        {{
        if [[ "{params.peaktype}" == "broad" ]]; then
            add="--broad"
        else
            add=""
        fi        
        printf "\nCalling {params.peaktype} peaks for paired-end {params.ipname} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAMPE -g {params.genomesize} {params.params} -n peaks_pe__{params.filetype}__{params.ipname} --outdir results/{params.env}/peaks/ ${{add}}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_macs2_se:
    input:
        ipfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        peakfile = "results/{env}/peaks/peaks_se__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}_peaks.{peaktype}Peak"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        ipname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{assign_chip_input(wildcards)}",
        peaktype = lambda wildcards: get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"]),
        filetype = lambda wildcards: {wildcards.file_type},
        env = lambda wildcards: {wildcards.env},
        params = config["chip_callpeaks"]['params'],
        genomesize = config[config['species']]['genomesize']
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "{file_type}__{peaktype}peak_calling", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["calling_peaks_macs2_se"]["threads"]
    resources:
        mem=config["resources"]["calling_peaks_macs2_se"]["mem"],
        tmp=config["resources"]["calling_peaks_macs2_se"]["tmp"]
    shell:
        """
        {{
        if [[ "{params.peaktype}" == "broad" ]]; then
            add="--broad"
        else
            add=""
        fi
        printf "\nCalling {params.peaktype} peaks for single-end {params.ipname} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAM -g {params.genomesize} {params.params} -n peaks_se__{params.filetype}__{params.ipname} --outdir results/{params.env}/peaks/ ${{add}}
        }} 2>&1 | tee -a "{log}"
        """
        
rule idr_analysis_replicates:
    input:
        peak_file = lambda wildcards: assign_peak_files_for_idr(wildcards)
    output:
        touch = "results/{env}/chkpts/idr__{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}.done"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis'),
        peaktype = lambda wildcards: get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"]),
        paired = lambda wildcards: get_sample_info_from_name(sample_name_str(wildcards, 'analysis'), analysis_samples, 'paired'),
        data_type = lambda wildcards: wildcards.data_type,
        line = lambda wildcards: wildcards.line,
        tissue = lambda wildcards: wildcards.tissue,
        sample_type = lambda wildcards: wildcards.sample_type,
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        replicate_pairs = lambda wildcards: get_replicate_pairs(wildcards)
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "IDR", ""))
    conda: CONDA_ENV
    threads: config["resources"]["idr_analysis_replicates"]["threads"]
    resources:
        mem=config["resources"]["idr_analysis_replicates"]["mem"],
        tmp=config["resources"]["idr_analysis_replicates"]["tmp"]
    shell:
        """
        {{
        printf "\nLooping over each unique pair of biological replicates for {params.sname} to perform IDR with:\n"
        idr --version
		if [[ "{params.paired}" == "PE" ]]; then
            pre="pe"
        else
            pre="se"
        fi
        for pair in {params.replicate_pairs}; do
            rep1=$(echo ${{pair}} | cut -d":" -f1)
            rep2=$(echo ${{pair}} | cut -d":" -f2)
            file1="results/{params.env}/peaks/peaks_${{pre}}__final__{params.data_type}__{params.line}__{params.tissue}__{params.sample_type}__${{rep1}}__{params.ref_genome}_peaks.{params.peaktype}Peak"
            file2="results/{params.env}/peaks/peaks_${{pre}}__final__{params.data_type}__{params.line}__{params.tissue}__{params.sample_type}__${{rep2}}__{params.ref_genome}_peaks.{params.peaktype}Peak"
            outfile="results/{params.env}/peaks/idr_${{pre}}__{params.data_type}__{params.line}__{params.tissue}__{params.sample_type}__${{rep1}}_vs_${{rep2}}__{params.ref_genome}_peaks.{params.peaktype}Peak"
            printf "\nPerforming IDR for ${{rep1}} vs ${{rep2}}\n"
            idr --input-file-type {params.peaktype}Peak --output-file-type {params.peaktype}Peak --samples ${{file1}} ${{file2}} -o ${{outfile}} -l results/{params.env}/reports/idr_{params.sname}.log --plot || true
            ## I think "|| true" is to avoid potential pipeline breaking errors if no positive peaks were found
            mv "${{outfile}}.png" results/{params.env}/plots/
        done
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule merging_chip_replicates:
    input:
        bamfiles = lambda wildcards: [ f"results/{wildcards.env}/mapped/final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.bam" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        temp_merge = temp("results/{env}/mapped/temp_merged__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"),
        mergefile = "results/{env}/mapped/merged__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis'),
        env = lambda wildcards: wildcards.env
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging_reps", ""))
    conda: CONDA_ENV
    threads: config["resources"]["merging_chip_replicates"]["threads"]
    resources:
        mem=config["resources"]["merging_chip_replicates"]["mem"],
        tmp=config["resources"]["merging_chip_replicates"]["tmp"]
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
		samtools merge -@ {threads} {output.temp_merge} {input.bamfiles}
		samtools sort -@ {threads} -o {output.mergefile} {output.temp_merge}
		samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """
        
rule making_pseudo_replicates:
    input:
        bamfile = lambda wildcards: f"results/{wildcards.env}/mapped/{'merged' if wildcards.replicate == 'merged' else 'final'}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam"
    output:
        temp_pseudo1 = temp("results/{env}/mapped/temp_pseudo1__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bam"),
        temp_pseudo2 = temp("results/{env}/mapped/temp_pseudo2__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bam"),   
        pseudo1 = temp("results/{env}/mapped/pseudo1__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bam"),
        pseudo2 = temp("results/{env}/mapped/pseudo2__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bam")
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis'),
        env = lambda wildcards: wildcards.env
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "splitting_pseudreps", ""))
    conda: CONDA_ENV
    threads: config["resources"]["making_pseudo_replicates"]["threads"]
    resources:
        mem=config["resources"]["making_pseudo_replicates"]["mem"],
        tmp=config["resources"]["making_pseudo_replicates"]["tmp"]
    shell:
        """
        {{
        printf "\nSplitting {params.sname} in two pseudo-replicates\n"
        samtools view -b -h -s 1.5 -@ {threads} -U {output.temp_pseudo2} -o {output.temp_pseudo1} {input.bamfile}
		samtools sort -@ {threads} -o {output.pseudo1} {output.temp_pseudo1}
		samtools sort -@ {threads} -o {output.pseudo2} {output.temp_pseudo2}
        }} 2>&1 | tee -a "{log}"
        """

rule best_peaks_pseudoreps:
    input:
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        peakfiles = input_peak_files_for_best_peaks
    output:
        bestpeaks = "results/{env}/peaks/selected_peaks__{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}.bed",
        stats_pseudoreps = temp("results/{env}/reports/stats_pseudoreps__{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}.txt")
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis'),
        env = lambda wildcards: wildcards.env,
        peaktype = lambda wildcards: get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"])
    log:
        temp(return_log_chip("{env}","{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "selecting_best_peaks", ""))
    conda: CONDA_ENV
    threads: config["resources"]["best_peaks_pseudoreps"]["threads"]
    resources:
        mem=config["resources"]["best_peaks_pseudoreps"]["mem"],
        tmp=config["resources"]["best_peaks_pseudoreps"]["tmp"]
    shell:
        """
        {{
        printf "\nIntersecting merged peaks ({input.peakfiles[0]}), and both pseudo replicates to select the best peaks for {params.sname}\n"
        awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[0]} | sort -k1,1 -k2,2n -u > "results/{params.env}/peaks/temp_{params.sname}_merged.bed"
		awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[1]} | sort -k1,1 -k2,2n -u > "results/{params.env}/peaks/temp_{params.sname}_pseudo1.bed"
		awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[2]} | sort -k1,1 -k2,2n -u > "results/{params.env}/peaks/temp_{params.sname}_pseudo2.bed"
		bedtools intersect -a results/{params.env}/peaks/temp_{params.sname}_pseudo1.bed -b results/{params.env}/peaks/temp_{params.sname}_pseudo2.bed > "results/{params.env}/peaks/temp_{params.sname}_pseudos.bed"
		bedtools intersect -a results/{params.env}/peaks/temp_{params.sname}_merged.bed -b results/{params.env}/peaks/temp_{params.sname}_pseudo1.bed -u > "results/{params.env}/peaks/temp_{params.sname}_selected.bed"
		bedtools intersect -a {input.peakfiles[0]} -b results/{params.env}/peaks/temp_{params.sname}_selected.bed -u > "results/{params.env}/peaks/selected_peaks_{params.sname}.{params.peaktype}Peak"
        printf "\nGetting best quality peaks peaks\n"
        ## Note: If broadpeak, an additional "summit" column will be added for potential downstream processes, which only represent the middle of the peak, not its summit.
        sort -k1,1 -k2,2n -k5nr results/{params.env}/peaks/selected_peaks_{params.sname}.{params.peaktype}Peak | awk -v OFS="\t" '{{print $1";"$2";"$3,$4,$5,$6,$7,$8,$9,$10}}' | awk 'BEGIN {{a=0}} {{b=$1; if (b!=a) print $0; a=$1}}' | awk -F"[;\t]" -v OFS="\t" -v t={params.peaktype} '{{if (t=="broad") $10=int(($3-$2)/2); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' | bedtools sort -g {input.chrom_sizes} > "{output.bestpeaks}"
        printf "\nExtracting peak stats for {params.sname}\n"
        merged=$(wc -l results/{params.env}/peaks/temp_{params.sname}_merged.bed | cut -d" " -f1)
		pseudos=$(awk '{{print $1,$2,$3}}' results/{params.env}/peaks/temp_{params.sname}_pseudos.bed | sort -k1,1 -k2,2n -u | wc -l)
		selected=$(cat results/{params.env}/peaks/temp_{params.sname}_selected.bed | sort -k1,1 -k2,2n -u | wc -l)
		printf "Merged=${{merged}}\nPseudos=${{pseudos}}\nSelected=${{selected}}\n" > "results/{params.env}/reports/stats_pseudoreps__{params.sname}.txt"
        rm -f "results/{params.env}/peaks/temp_{params.sname}"*
        }} 2>&1 | tee -a "{log}"
        """    

rule make_peak_stats:
    input:
        logs = lambda wildcards: define_logs_final_input(wildcards),
        stats_pseudoreps = lambda wildcards: f"results/{wildcards.env}/reports/stats_pseudoreps__{wildcards.sample_name}.txt"
        ## maybe a better solution is to append a stat file with wc -l as they are generated, or to create a new stat file for each file, then accessible in bash by regex on the samplename
    output:
        stat_file = "results/{env}/reports/summary_ChIP_peak_stats_{sample_name}.txt",
        log = "results/{env}/logs/called_peaks__{sample_name}.log"
    wildcard_constraints:
        env = "ChIP|TF"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        paired = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'paired'),
        peaktype = lambda wildcards: get_peaktype(get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'sample_type'), config["chip_callpeaks"]["peaktype"]),
        line = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'line'),
        tissue = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'tissue'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'sample_type'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome'),
        rep1 = lambda wildcards: get_replicate_name(wildcards, 0),
        rep2 = lambda wildcards: get_replicate_name(wildcards, 1)
    threads: config["resources"]["make_peak_stats"]["threads"]
    resources:
        mem=config["resources"]["make_peak_stats"]["mem"],
        tmp=config["resources"]["make_peak_stats"]["tmp"]
    shell:
        """
        nrep1=$(awk '{{print $1,$2,$3}}' {params.rep1} | sort -k1,1 -k2,2n -u | wc -l)
        if [[ "{params.rep2}" == "missingrep" ]]; then
            nrep2=0
        else
            nrep2=$(awk '{{print $1,$2,$3}}' {params.rep2} | sort -k1,1 -k2,2n -u | wc -l)
        fi
        merged=$(grep "Merged" {input.stats_pseudoreps} | cut -d"=" -f2)
        pseudos=$(grep "Pseudos" {input.stats_pseudoreps} | cut -d"=" -f2)
        selected=$(grep "Selected" {input.stats_pseudoreps} | cut -d"=" -f2)
        ## Need to add idr if present; limited to 2 reps for now: Estimate max number of reps first, then filling all samples with 0? Easier in R?
        printf "Line\tTissue\tMark\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tCommon_peaks\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.ref_genome} -v a=${{nrep1}} -v b=${{nrep2}} -v c=${{merged}} -v d=${{pseudos}} -v e=${{selected}} 'BEGIN {{if (c==0) {{x=a}} else {{x=c}}; print l,t,m,r,a,b,c,d,e" ("e/x*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        """

rule all_chip:
    input:
        final = lambda wildcards: define_final_chip_output(wildcards.ref_genome)
    output:
        touch = "results/{env}/chkpts/ChIP_analysis__{analysis_name}__{ref_genome}.done"
    wildcard_constraints:
        env = "ChIP|TF"
    localrule: True
    shell:
        """
        touch {output.touch}
        """
       