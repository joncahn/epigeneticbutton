# function to access logs more easily
def return_log_rna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"RNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

# def define_RNA_input_for_deg(ref_genome):
    # data_type = get_sample_info_from_name(sample_name, analysis_samples, 'data_type')
    # line = get_sample_info_from_name(sample_name, analysis_samples, 'line')
    # tissue = get_sample_info_from_name(sample_name, analysis_samples, 'tissue')
    # sample_type = get_sample_info_from_name(sample_name, analysis_samples, 'sample_type')
    # ref_genome = get_sample_info_from_name(sample_name, analysis_samples, 'ref_genome')
    # replicates = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])
    
    # return [ f"mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.deduplicated.CX_report.txt.gz"
                    # for replicate in replicates ]

def define_final_rna_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    map_files = []
    bigwig_files = []
    qc_files = []
    deg_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'RNA') & (samples['ref_genome'] == ref_genome) ]
    
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')        
        paired = get_sample_info_from_name(sname, samples, 'paired')
        if paired == "PE":
            map_files.append(f"RNA/logs/process_rna_pe_sample__{sname}.log")
            qc_files.append(f"RNA/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
            qc_files.append(f"RNA/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            qc_files.append(f"RNA/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
            qc_files.append(f"RNA/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
        else:
            map_files.append(f"RNA/logs/process_rna_se_sample__{sname}.log")
            qc_files.append(f"RNA/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file
            qc_files.append(f"RNA/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
        
        bigwig_files.append(f"RNA/tracks/{sname}__plus.bw")
        bigwig_files.append(f"RNA/tracks/{sname}__minus.bw")
        
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'RNA') & (analysis_samples['ref_genome'] == ref_genome) ]
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            bigwig_files.append(f"RNA/tracks/merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}__plus.bw")
            bigwig_files.append(f"RNA/tracks/merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}__minus.bw")
    
    # deg_files.append(f"RNA/DEGs/summary.txt")
    
    results = map_files
    
    if qc_option == "all":
        results += qc_files
        
    if analysis:
        results += bigwig_files
        # results += bigwig_files + deg_files
    return results
        
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/rna.yaml")

rule make_STAR_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa",
        gtf = "genomes/{ref_genome}/temp_{ref_genome}.gtf"
    output:
        indices = directory("genomes/{ref_genome}/STAR_index")
    params:
        star_index = config[config['species']]['star_index']
    log:
        temp(os.path.join(REPO_FOLDER,"logs","STAR_index_{ref_genome}.log"))
    conda: CONDA_ENV
    threads: config["resources"]["STAR_indices"]["threads"]
    resources:
        mem=config["resources"]["STAR_indices"]["mem"],
        tmp=config["resources"]["STAR_indices"]["tmp"]
    shell:
        """
        {{
        printf "\nBuilding STAR index directory for {wildcards.ref_genome}\n"
        mkdir "{output.indices}"
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir "{output.indices}" --genomeFastaFiles "{input.fasta}" --sjdbGTFfile "{input.gtf}" {params.star_index}
        }} 2>&1 | tee -a "{log}"
        """

rule STAR_map_pe:
    input:
        fastq1 = "RNA/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "RNA/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        bamfile = temp("RNA/mapped/star_pe__{sample_name}_Aligned.out.bam")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        file_order = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['file_order'],
        prefix = lambda wildcards: f"RNA/mapped/star_pe__{wildcards.sample_name}_"
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["STAR_map"]["threads"]
    resources:
        mem=config["resources"]["STAR_map"]["mem"],
        tmp=config["resources"]["STAR_map"]["tmp"]
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
        if [[ {params.fileorder} == "rampage"]]; then
            printf "Input file order for RAMPAGE (R2 R1)\n"
            input="{input.fastq2}" "{input.fastq1}"
        else
            printf "Input file order for RNAseq (R1 R2)\n"
            input="{input.fastq1}" "{input.fastq2}"
        fi
        STAR --version
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn ${{input}} --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 20 --quantMode GeneCounts
        mv "RNA/mapped/star_pe__{params.sample_name}_Log.final.out" "{output.metrics_map}"
        rm -f RNA/mapped/*"{params.sample_name}_Log"*
        }} 2>&1 | tee -a "{log}"
        """    

rule STAR_map_se:
    input:
        fastq0 = "RNA/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        bamfile = temp("RNA/mapped/star_se__{sample_name}_Aligned.out.bam")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"RNA/mapped/star_se__{wildcards.sample_name}_"
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["STAR_map"]["threads"]
    resources:
        mem=config["resources"]["STAR_map"]["mem"],
        tmp=config["resources"]["STAR_map"]["tmp"]
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
        STAR --version
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn "{input.fastq0}" --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --quantMode GeneCounts
        mv "RNA/mapped/star_se__{params.sample_name}_Log.final.out" "{output.metrics_map}"
        rm -f RNA/mapped/*"{params.sample_namparams.sample_name}_Log"*
        }} 2>&1 | tee -a "{log}"
        """
        
rule filter_rna_pe:
    input:
        bamfile = "RNA/mapped/star_pe__{sample_name}_Aligned.out.bam"
    output:
        mrkdup=temp("RNA/mapped/star_pe__{sample_name}_Processed.out.bam"),
        sorted_file=temp("RNA/mapped/star_pe__{sample_name}_Processed.sorted.out.bam"),
        metrics_flag = "RNA/reports/flagstat_pe__{sample_name}.txt",
        metrics_map = "RNA/reports/star_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_rna("{sample_name}", "filteringRNA", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["filter_rna"]["threads"]
    resources:
        mem=config["resources"]["filter_rna"]["mem"],
        tmp=config["resources"]["filter_rna"]["tmp"]
    shell:
        """
        {{
        ### Marking duplicates
        ## Errors can happen because of limitBAMsortRAM, which seem to happen when bam files are sorted by coordinates (now removed from mapping step). Might want parameters from sorting duplicates too.
        printf "\nMarking duplicates\n"
        STAR --runMode inputAlignmentsFromBAM --inputBAMfile "{input.bamfile}" --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix "RNA/mapped/star_pe__{params.sample_name}_"
        #### Indexing bam file
        printf "\nSorting bam file\n"
        samtools sort -@ {threads} "{output.mrkdup}" -o "{output.sorted_file}"
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "{output.sorted_file}"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.sorted_file}" > "{output.metrics_flag}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_rna_se:
    input:
        bamfile = "RNA/mapped/star_se__{sample_name}_Aligned.out.bam"
    output:
        sorted_file=temp("RNA/mapped/star_se__{sample_name}_Aligned.sorted.out.bam"),
        metrics_flag = "RNA/reports/flagstat_se_{sample_name}.txt",
        metrics_map = "RNA/reports/star_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_rna("{sample_name}", "filteringRNA", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["filter_rna"]["threads"]
    resources:
        mem=config["resources"]["filter_rna"]["mem"],
        tmp=config["resources"]["filter_rna"]["tmp"]
    shell:
        """
        {{
        #### Sorting bam file
        printf "\nSorting bam file\n"
        samtools sort -@ {threads} "{input.bamfile}" -o "{output.sorted_file}"
        #### Indexing bam file
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "{output.sorted_file}"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.sorted_file}" > "{output.metrics_flag}"
        }} 2>&1 | tee -a "{log}"
        """        

rule make_rna_stats_pe:
    input:
        metrics_trim = "RNA/reports/trim_pe__{sample_name}.txt",
        metrics_map = "RNA/reports/star_pe__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_rna(wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingSTAR", "filteringRNA"] ]
    output:
        stat_file = "RNA/reports/summary_RNA_PE_mapping_stats_{sample_name}.txt",
        log = "RNA/logs/process_rna_pe_sample__{sample_name}.log"
    params:
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Number of reads mapped to multiple loci" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Uniquely mapped reads number" "{input.metrics_map}" | awk '{{print $NF}}')
        allmap=$((multi+single))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """
        
rule make_rna_stats_se:
    input:
        metrics_trim = "RNA/reports/trim_se__{sample_name}.txt",
        metrics_map = "RNA/reports/star_se__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_rna(wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingSTAR", "filteringRNA"] ]
    output:
        stat_file = "RNA/reports/summary_RNA_SE_mapping_stats_{sample_name}.txt",
        log = "RNA/logs/process_rna_se_sample__{sample_name}.log"
    params:
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Number of reads mapped to multiple loci" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Uniquely mapped reads number" "{input.metrics_map}" | awk '{{print $NF}}')
        allmap=$((multi+single))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule pe_or_se_rna_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "filter_rna", "sorted_file")
    output:
        bam = "RNA/mapped/{sample_name}.bam",
        touch = "RNA/chkpts/map_rna__{sample_name}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        mv {input} {output.bam}
        mv {input}.bai {output.bam}.bai
        touch {output.touch} 
        """

rule merging_rna_replicates:
    input:
        bamfiles = lambda wildcards: [ f"RNA/mapped/{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.bam" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        mergefile = "RNA/mapped/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis')
    log:
        temp(return_log_rna("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging_rna_reps", ""))
    conda: CONDA_ENV
    threads: config["resources"]["merging_rna_replicates"]["threads"]
    resources:
        mem=config["resources"]["merging_rna_replicates"]["mem"],
        tmp=config["resources"]["merging_rna_replicates"]["tmp"]
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
		samtools merge -@ {threads} RNA/mapped/temp_{params.sname}.bam {input.bamfiles}
		samtools sort -@ {threads} -o {output.mergefile} RNA/mapped/temp_{params.sname}.bam
		rm -f RNA/mapped/temp_{params.sname}.bam
		samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """

rule make_rna_stranded_bigwigs:
    input: 
        bamfile = "RNA/mapped/{sample_name}.bam",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        bw_plus = "RNA/tracks/{sample_name}__plus.bw",
        bw_minus = "RNA/tracks/{sample_name}__minus.bw"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        param_bg = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['param_bg'],
        strandedness = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['strandedness'],
        multimap = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['multimap']
    conda: CONDA_ENV
    threads: config["resources"]["filter_rna"]["threads"]
    resources:
        mem=config["resources"]["filter_rna"]["mem"],
        tmp=config["resources"]["filter_rna"]["tmp"]
    shell:
        """
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --runThreadN {threads} --inputBAMfile "{output.sorted_file}" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "RNA/tracks/bg_{params.sample_name}_"
        ### Converting to bigwig files
        printf "\nConverting bedGraphs to bigWigs\n"
        if [[ {params.multimap} == "multiple" ]]; then
            bed1="RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str1.out.bg"
            bed2="RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str2.out.bg"
        elif [[ {params.multimap} == "unique" ]]; then
            bed1="RNA/tracks/bg_{params.sample_name}_Signal.Unique.str1.out.bg"
            bed2="RNA/tracks/bg_{params.sample_name}_Signal.Unique.str2.out.bg"
        fi        
        bedSort ${{bed1}} "RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg"
        bedSort ${{bed2}} "RNA/tracks/{params.sample_name}_Signal.sorted.str2.out.bg"
        if [[ "{params.strandedness}" == "forward" ]]; then
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.str2.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
        elif [[ "{params.strandedness}" == "reverse" ]]; then
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.str2.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
        fi
        rm -f RNA/tracks/*"{params.sample_name}_Signal"*
        rm -f RNA/tracks/*"{params.sample_name}_Log"*
        """

# rule call_all_DEGs:
    # input:
        # sample1 = lambda wildcards: define_DMR_samples(wildcards.sample1),
        # sample2 = lambda wildcards: define_DMR_samples(wildcards.sample2),
        # chrom_sizes = lambda wildcards: f"genomes/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}/chrom.sizes"
    # output:
        # deg_summary = "RNA/DEGs/summary.txt"
    # params:
        # script = os.path.join(REPO_FOLDER,"scripts/R_call_DMRs.R"),
        # context = config['mC_context'],
        # sample1 = lambda wildcards: wildcards.sample1,
        # sample2 = lambda wildcards: wildcards.sample2,
        # nb_sample1 = lambda wildcards: len(define_DMR_samples(wildcards.sample1)),
        # nb_sample2 = lambda wildcards: len(define_DMR_samples(wildcards.sample2))
    # log:
        # temp(return_log_mc("{sample1}__vs__{sample2}", "DMRs", ""))
    # conda: os.path.join(REPO_FOLDER,"envs/call_dmrs.yaml")
    # threads: config["resources"]["call_dmrs"]["threads"]
    # resources:
        # mem=config["resources"]["call_dmrs"]["mem"],
        # tmp=config["resources"]["call_dmrs"]["tmp"]
    # shell:
        # """
        # printf "running DMRcaller for {params.sample1} vs {params.sample2}\n"
        # Rscript "{params.script}" "{threads}" "{input.chrom_sizes}" "{params.context}" "{params.sample1}" "{params.sample2}" "{params.nb_sample1}" "{params.nb_sample2}" {input.sample1} {input.sample2}
        # """

rule all_rna:
    input:
        lambda wildcards: define_final_rna_output(wildcards.ref_genome)
    output:
        touch = "RNA/chkpts/RNA_analysis__{ref_genome}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        touch {output.touch}
        """        
