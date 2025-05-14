# function to access logs more easily
def return_log_chip(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"ChIP","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def assign_mapping_paired(wildcards, rulename, outputfile):
    sname = wildcards.sample_name
    paired = get_sample_info_from_name(sname, samples, 'paired')
    if paired == "PE":
        rule_obj = getattr(rules, f"{rulename}_pe")
    elif paired == "SE":
        rule_obj = getattr(rules, f"{rulename}_se")
        
    return getattr(rule_obj.output, outputfile).format(sample_name=sname)

def assign_chip_input(wildcards):
    inputname = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__Input__{wildcards.replicate}__{wildcards.ref_genome}"
    if wildcards.file_type in ['merged', 'pseudo1', 'pseudo2']:
        return inputname
    elif inputname in samples['sample_name']:
        return inputname
    else:
        ipname = sample_name(wildcards, 'sample')
        ippaired = get_sample_info_from_name(ipname, samples, 'paired')
        alts = []
        for rep in chip_input_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.ref_genome), []):
            alt_inputname = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__Input__{rep}__{wildcards.ref_genome}"
            alts.append(f"{alt_inputname}")
            if get_sample_info_from_name(alt_inputname, samples, 'paired') == ippaired:
                return alt_inputname
            else:
                alts.append(f"{ippaired}")
        
        raise ValueError(f"\nSample '{ipname}' does not have corresponding Input among:\n{alts}")
                
def get_peaktype(sample_type, peaktype_config):
    for pattern, peaktype in peaktype_config.items():
        if re.search(pattern, sample_type):
            return peaktype
    raise ValueError(f"\nNo peaktype found for sample_type '{sample_type}")

def define_final_chip_output(ref_genome):
    peak_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'ChIP') & (samples['ref_genome'] == ref_genome) & (samples['sample_type'] != "Input") ]
    
    for _, row in filtered_rep_samples.iterrows():
        peaktype = get_peaktype(row.sample_type, config["chip_callpeaks"]['peaktype'])
        sname = sample_name(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        if paired == "PE":
            peak_files.append(f"ChIP/peaks/peaks_pe__final__{sname}.{peaktype}Peak")
        else:
            peak_files.append(f"ChIP/peaks/peaks_se__final__{sname}.{peaktype}Peak")
        
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            if paired == "PE":
                peak_files.append(f"ChIP/peaks/peaks_pe__merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.{peaktype}Peak")
            else:
                peak_files.append(f"ChIP/peaks/peaks_se__merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.{peaktype}Peak")
   
    return peak_files
        
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
        bamfile = temp("ChIP/mapped/mapped_pe__{sample_name}.bam"),
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
        bamfile = temp("ChIP/mapped/mapped_se__{sample_name}.bam"),
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
        metrics_trim = "ChIP/reports/trim_pe__{sample_name}.txt",
        metrics_map = "ChIP/reports/bt2_pe__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "ChIP/logs/process_pe_sample__{sample_name}.log"
    params:
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned concordantly >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned concordantly exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{input.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule make_chip_stats_se:
    input:
        stat_file = f"ChIP/reports/summary_mapping_stats_{analysis_name}.txt",
        metrics_trim = "ChIP/reports/trim_se__{sample_name}.txt",
        metrics_map = "ChIP/reports/bt2_se__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "ChIP/logs/process_se_sample__{sample_name}.log"
    params:
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total reads processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{input.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule pe_or_se_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "filter_chip", "bamfile")
    output:
        bam = "ChIP/mapped/final__{sample_name}.bam",
        touch = "ChIP/chkpts/map__{sample_name}.done"
    shell:
        """
        mv {input} {output.bam}
        touch {output.touch} 
        """
    
rule make_coverage_chip:
    input: 
        bamfile = "ChIP/mapped/final__{sample_name}.bam"
    output:
        bigwigcov = "ChIP/tracks/coverage__{sample_name}.bw"
    params:
        binsize = config['chip_tracks']['binsize']
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        bamCoverage -b {input.bamfile} -o {output.bigwigcov} -bs {params.binsize} -p {threads}
        """

rule make_bigwig_chip:
    input: 
        ipfile = lambda wildcards: f"ChIP/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"ChIP/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        bigwigfile = "ChIP/tracks/FC__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bw"
    params:
        ipname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{assign_chip_input(wildcards)}",
        binsize = config['chip_tracks']['binsize'],
        params = config['chip_tracks']['params']
    log:
        temp(return_log_chip("{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "making_bigwig", "either"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nCalling {params.peaktype} peaks for {params.ipname} (vs {params.inputname}) using macs2 version:\n"
        bamCompare -b1 {ipfile} -b2 {inputfile} -o {output.bigwigfile} -p {threads} --binSize {params.binsize} {params.params}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_macs2_pe:
    input:
        ipfile = lambda wildcards: f"ChIP/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"ChIP/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        peakfile = "ChIP/peaks/peaks_pe__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.{peaktype}Peak"
    params:
        ipname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{assign_chip_input(wildcards)}",
        peaktype = lambda wildcards: get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"]),
        params = config["chip_callpeaks"]['params'],
        genomesize = config["chip_callpeaks"]['genomesize']
    log:
        temp(return_log_chip("{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "{file_type}__{peaktype}peak_calling", "PE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nCalling {params.peaktype} peaks for {params.ipname} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t ${ipfile} -c {inputfile} -f BAMPE -g {params.genomesize} {params.params} -n {param.ipname} --outdir ChIP/peaks/ --{params.peaktype}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_macs2_se:
    input:
        ipfile = lambda wildcards: f"ChIP/mapped/{wildcards.file_type}__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}.bam",
        inputfile = lambda wildcards: f"ChIP/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        peakfile = "ChIP/peaks/peaks_se__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.{peaktype}Peak"
    params:
        ipname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        inputname = lambda wildcards: f"{assign_chip_input(wildcards)}",
        peaktype = lambda wildcards: get_peaktype(wildcards.sample_type, config["chip_callpeaks"]["peaktype"]),
        params = config["chip_callpeaks"]['params'],
        genomesize = config["chip_callpeaks"]['genomesize']
    log:
        temp(return_log_chip("{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "{file_type}__{peaktype}peak_calling", "SE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nCalling {params.peaktype} peaks for {params.ipname} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t ${ipfile} -c {inputfile} -f BAM -g {params.genomesize} {params.params} -n {param.ipname} --outdir ChIP/peaks/ --{params.peaktype}
        }} 2>&1 | tee -a "{log}"
        """
        
# rule merging_inputs:
    # input:
        # bamfiles = lambda wildcards: [ f"ChIP/mapped/final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__Input__{replicate}__{wildcards.ref_genome}.bam" 
                                      # for replicate in chip_input_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.ref_genome), []) ]
    # output:
        # mergefile = "ChIP/mapped/merged__{data_type}__{line}__{tissue}__Input__{ref_genome}.bam"
    # params:
        # sname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__Input__{wildcards.ref_genome}"
    # log:
        # temp(return_log_chip("{data_type}__{line}__{tissue}__Input__{ref_genome}", "merging", "merged"))
    # conda: CONDA_ENV
    # threads: workflow.cores
    # shell:
        # """
        # {{
        # printf "\nMerging replicates of {params.sname}\n"
		# samtools merge -@ {threads} ChIP/mapped/temp_{params.sname}.bam {input.bamfiles}
		# samtools sort -@ {threads} -o {output.mergefile} ChIP/mapped/temp_{params.sname}.bam
		# rm -f ChIP/mapped/temp_{params.sname}.bam
		# samtools index -@ {threads} {output.mergefile}
        # }} 2>&1 | tee -a "{log}"
        # """
        
rule merging_replicates:
    input:
        bamfiles = lambda wildcards: [ f"ChIP/mapped/final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.bam" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        mergefile = "ChIP/mapped/merged__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"
    params:
        sname = lambda wildcards: sample_name(wildcards, 'analysis')
    log:
        temp(return_log_chip("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging", "merged"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
		samtools merge -@ {threads} ChIP/mapped/temp_{params.sname}.bam {input.bamfiles}
		samtools sort -@ {threads} -o {output.mergefile} ChIP/mapped/temp_{params.sname}.bam
		rm -f ChIP/mapped/temp_{params.sname}.bam
		samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """

rule ChIP_all:
    input:
        lambda wildcards: define_final_chip_output(wildcards.ref_genome)
    output:
        touch = "ChIP/chkpts/ChIP_analysis__{ref_genome}.done"
    shell:
        """
        touch {output.touch}
        """