# function to access logs more easily
def return_log_rna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"RNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

# def define_RNA_samples(sample_name):
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
    final_files = []
    qc_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'RNA') & (samples['ref_genome'] == ref_genome) ]
    
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        if paired == "PE":
            final_files.append(f"RNA/logs/process_rna_pe_sample__{sname}.log")
            qc_files.append(f"RNA/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
            qc_files.append(f"RNA/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            qc_files.append(f"RNA/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
            qc_files.append(f"RNA/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
        else:
            final_files.append(f"RNA/logs/process_rna_se_sample__{sname}.log")
            qc_files.append(f"RNA/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file
            qc_files.append(f"RNA/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
        
    if qc_option == "all":
        return final_files + qc_files
    else:
        return final_files
        
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
        bamfile = temp("RNA/mapped/map_pe__{sample_name}_Aligned.out.bam"),
        touch = "RNA/chkpts/temp_pe__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"RNA/mapped/map_pe__{wildcards.sample_name}_"
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
        STAR --version
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn "{input.fastq1}" "{input.fastq2}" --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 20 --quantMode GeneCounts
        touch "{output.touch}"
        }} 2>&1 | tee -a "{log}"
        """    

rule STAR_map_se:
    input:
        fastq0 = "RNA/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        bamfile = temp("RNA/mapped/map_se__{sample_name}_Aligned.out.bam"),
        touch = "RNA/chkpts/temp_se__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"RNA/mapped/map_se__{wildcards.sample_name}_"
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
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """
        
rule filter_rna_pe:
    input:
        bamfile = "RNA/mapped/map_pe__{sample_name}_Aligned.out.bam",
        touch = "RNA/chkpts/temp_pe__{sample_name}.done",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        mrkdup=temp("RNA/mapped/map_pe__{sample_name}_Processed.out.bam"),
        sorted_file="RNA/mapped/map_pe__{sample_name}_Processed.sorted.out.bam",
        bw_plus = "RNA/tracks/{sample_name}_plus.bw",
        bw_minus = "RNA/tracks/{sample_name}_minus.bw",
        metrics_flag = "RNA/reports/flagstat_pe__{sample_name}.txt",
        metrics_map = "RNA/reports/star_pe__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        param_bg = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['param_bg'],
        strandedness = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['strandedness']
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
        STAR --runMode inputAlignmentsFromBAM --inputBAMfile "{input.bamfile}" --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix "RNA/mapped/map_pe__{params.sample_name}_"
        #### Indexing bam file
        printf "\nSorting bam file\n"
        samtools sort -@ {threads} "{output.mrkdup}" -o "{output.sorted_file}"
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "{output.sorted_file}"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.sorted_file}" > "{output.metrics_flag}"
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --runThreadN {threads} --inputBAMfile "{output.sorted_file}" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "RNA/tracks/bg_{params.sample_name}_"
        ### Converting to bigwig files
        printf "\nConverting bedGraphs to bigWigs\n"
        bedSort "RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str1.out.bg" "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str1.out.bg"
        bedSort "RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str2.out.bg" "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str2.out.bg"
        if [[ "{params.strandedness}" == "forward" ]]; then
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str1.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str2.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
        elif [[ "{params.strandedness}" == "reverse" ]]; then
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str1.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str2.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
        fi	
        mv "RNA/mapped/map_pe__{params.sample_name}_Log.final.out" "{output.metrics_map}"
        ### Cleaning up
        printf "\nCleaning up\n" ## Could add an option to clean all log output (e.g. Log.out, Log.progress.out, SJ.out.tab) or not.
        rm -f RNA/tracks/*"{params.sample_name}_Signal"*
        rm -f RNA/mapped/*"{params.sample_name}_Log"*
        rm -f RNA/tracks/*"{params.sample_name}_Log"*
        }} 2>&1 | tee -a "{log}"
        """

rule filter_rna_se:
    input:
        bamfile = "RNA/mapped/map_se__{sample_name}_Aligned.out.bam",
        touch = "RNA/chkpts/temp_se__{sample_name}.done",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        sorted_file="RNA/mapped/map_se__{sample_name}_Aligned.sorted.out.bam",
        bw_plus = "RNA/tracks/{sample_name}_plus.bw",
        bw_minus = "RNA/tracks/{sample_name}_minus.bw",
        metrics_flag = "RNA/reports/flagstat_se_{sample_name}.txt",
        metrics_map = "RNA/reports/star_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        param_bg = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['param_bg'],
        strandedness = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['strandedness']
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
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --runThreadN {threads} --inputBAMfile "{output.sorted_file}" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "RNA/tracks/bg_{params.sample_name}_"
        ### Converting to bigwig files
        printf "\nConverting bedGraphs to bigWigs\n"
        bedSort "RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str1.out.bg" "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str1.out.bg"
        bedSort "RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str2.out.bg" "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str2.out.bg"
        if [[ "{params.strandedness}" == "forward" ]]; then
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str1.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str2.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
        elif [[ "{params.strandedness}" == "reverse" ]]; then
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str1.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
            bedGraphToBigWig "RNA/tracks/{params.sample_name}_Signal.sorted.UniqueMultiple.str2.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
        fi	
        mv "RNA/mapped/map_se__{params.sample_name}_Log.final.out" "{output.metrics_map}"
        ### Cleaning up
        printf "\nCleaning up\n"
        rm -f RNA/tracks/*"{params.sample_name}_Signal"*
        rm -f RNA/mapped/*"{params.sample_namparams.sample_name}_Log"*
        rm -f RNA/tracks/*"{params.sample_name}_Log"*
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

rule dispatch_pair_map_only_rna:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "make_rna_stats", "log")
    output:
        touch = "RNA/chkpts/map__{sample_name}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        touch {output.touch}
        """

# rule call_all_DEGs:
    # input:
        # sample1 = lambda wildcards: define_DMR_samples(wildcards.sample1),
        # sample2 = lambda wildcards: define_DMR_samples(wildcards.sample2),
        # chrom_sizes = lambda wildcards: f"genomes/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}/chrom.sizes"
    # output:
        # dmr_summary = "mC/DMRs/summary__{sample1}__vs__{sample2}__DMRs.txt"
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
