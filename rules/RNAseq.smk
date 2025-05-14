# function to access logs more easily
def return_log_rna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"RNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")    

def define_final_rna_output(ref_genome):
    final_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'RNA') & (samples['ref_genome'] == ref_genome) ]
    
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        if paired == "PE":
            final_files.append(f"RNA/logs/process_pe_sample__{sname}.log")
        else:
            final_files.append(f"RNA/logs/process_pe_sample__{sname}.log")
        
    return final_files
        
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/rna.yaml")

rule stat_file_rna:
    output:
        stat_file = f"RNA/reports/summary_mapping_stats_{analysis_name}.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """

rule make_RNA_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa",
        gtf = "genomes/{ref_genome}/temp_{ref_genome}.gtf"
    output:
        indices = directory("genomes/{ref_genome}/STAR_index")
    params:
        star_index = config['star_index']
    log:
        os.path.join(REPO_FOLDER,"logs","STAR_index_{ref_genome}.log")
    conda: CONDA_ENV
    threads: workflow.cores
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
        bamfile = "RNA/mapped/map_pe__{sample_name}_Aligned.sortedByCoord.out.bam",
        touch = "RNA/chkpts/temp_pe__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"RNA/mapped/map_pe__{wildcards.sample_name}_"
    log:
        temp(return_log_rna("{sample_name}", "mapping", "PE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
        STAR --version
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn "{input.fastq1}" "{input.fastq2}" --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 20 --quantMode GeneCounts
        touch "{output.touch}"
        }} 2>&1 | tee -a "{log}"
        """    

rule STAR_map_se:
    input:
        fastq0 = "RNA/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        bamfile = "RNA/mapped/map_se__{sample_name}_Aligned.sortedByCoord.out.bam",
        touch = "RNA/chkpts/temp_se__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"RNA/mapped/map_se__{wildcards.sample_name}_"
    log:
        temp(return_log_rna("{sample_name}", "mapping", "SE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
        STAR --version
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn "{input.fastq0}" --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --quantMode GeneCounts
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """
        
rule filter_rna_pe:
    input:
        bamfile = "RNA/mapped/map_pe__{sample_name}_Aligned.sortedByCoord.out.bam",
        touch = "RNA/chkpts/temp_pe__{sample_name}.done",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
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
        temp(return_log_rna("{sample_name}", "filtering", "PE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        ### Marking duplicates
        printf "\nMarking duplicates\n"
        STAR --runMode inputAlignmentsFromBAM --inputBAMfile "{input.bamfile}" --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix "RNA/mapped/mrkdup_{params.sample_name}_"
        #### Indexing bam file
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "RNA/mapped/mrkdup_{params.sample_name}_Processed.out.bam"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "RNA/mapped/mrkdup_{params.sample_name}_Processed.out.bam" > "{output.metrics_flag}"
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --inputBAMfile "RNA/mapped/mrkdup_{params.sample_name}_Processed.out.bam" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "RNA/tracks/bg_{params.sample_name}_"
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
        bamfile = "RNA/mapped/map_se__{sample_name}_Aligned.sortedByCoord.out.bam",
        touch = "RNA/chkpts/temp_se__{sample_name}.done",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
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
        temp(return_log_rna("{sample_name}", "filtering", "SE"))
    conda: CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        #### Indexing bam file
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "{input.bamfile}"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{input.bamfile}" > "{output.metrics_flag}"
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --inputBAMfile "{input.bamfile}" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "RNA/tracks/bg_{params.sample_name}_"
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
        stat_file = f"RNA/reports/summary_mapping_stats_{analysis_name}.txt",
        metrics_trim = "RNA/reports/trim_pe__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        metrics_map = "RNA/reports/star_pe__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        logs = lambda wildcards: [ return_log_rna(sample_name(wildcards), step, get_sample_info(wildcards, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "RNA/logs/process_pe_sample__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Number of reads mapped to multiple loci" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Uniquely mapped reads number" "{input.metrics_map}" | awk '{{print $NF}}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={wildcards.line} -v t={wildcards.tissue} -v m={wildcards.sample_type} -v r={wildcards.replicate} -v g={wildcards.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{input.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """
        
rule make_rna_stats_se:
    input:
        stat_file = f"RNA/reports/summary_mapping_stats_{analysis_name}.txt",
        metrics_trim = "RNA/reports/trim_se__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        metrics_map = "RNA/reports/star_se__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.txt",
        logs = lambda wildcards: [ return_log_rna(sample_name(wildcards), step, get_sample_info(wildcards, 'paired')) for step in ["downloading", "trimming", "mapping", "filtering"] ]
    output:
        log = "RNA/logs/process_se_sample__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.log"
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Number of reads mapped to multiple loci" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Uniquely mapped reads number" "{input.metrics_map}" | awk '{{print $NF}}')
        allmap=$((multi+single))
        awk -v OFS="\t" -v l={wildcards.line} -v t={wildcards.tissue} -v m={wildcards.sample_type} -v r={wildcards.replicate} -v g={wildcards.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{input.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule dispatch_pair_map_only_rna:
    input:
        lambda wildcards: assign_mapping_paired(wildcards.sample_name, "make_rna_stats", "log")
    output:
        touch = "RNA/chkpts/map__{sample_name}.done"
    shell:
        """
        touch {output.touch}
        """

rule all_rna:
    input:
        lambda wildcards: define_final_rna_output(wildcards.ref_genome)
    output:
        touch = "RNA/chkpts/RNA_analysis__{ref_genome}.done"
    shell:
        """
        touch {output.touch}
        """        
