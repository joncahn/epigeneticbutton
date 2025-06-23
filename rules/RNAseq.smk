# function to access logs more easily
def return_log_rna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"RNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def define_RNA_input_for_degs(ref_genome):
    file_paths = []
    filtered_samples = samples[ (samples['data_type'] == 'RNAseq') & (samples['ref_genome'] == ref_genome) ].copy()
    return [f"RNA/DEG/counts__{sname}.tab" for sname in filtered_samples['sample_name']]

def define_final_rna_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    analysis_name = config['analysis_name']
    map_files = []
    bigwig_files = []
    qc_files = []
    deg_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'RNA') & (samples['ref_genome'] == ref_genome) ].copy()
    
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
        
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'RNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            bigwig_files.append(f"RNA/tracks/{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}__plus.bw")
            bigwig_files.append(f"RNA/tracks/{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}__minus.bw")
    
    filtered_analysis_samples2 = samples[ (samples['data_type'] == 'RNAseq') & (samples['ref_genome'] == ref_genome) ].copy()
    filtered_analysis_samples2['Sample'] = filtered_analysis_samples2['line'] + "__" + filtered_analysis_samples2['tissue']
    if len(filtered_analysis_samples2['Sample'].drop_duplicates()) >= 2:   
        deg_files.append(f"RNA/chkpts/calling_DEGs__{analysis_name}__{ref_genome}.done")
        deg_files.append(f"RNA/DEG/genes_rpkm__{analysis_name}__{ref_genome}.txt")
    elif len(filtered_analysis_samples2['Sample'].drop_duplicates()) == 1:
        deg_files.append(f"RNA/DEG/genes_rpkm__{analysis_name}__{ref_genome}.txt")
        
    results = map_files
    
    if qc_option == "all":
        results += qc_files
        
    if analysis:
        results += bigwig_files + deg_files

    return results
        
rule make_STAR_indices:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf"
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
        bamfile = temp("RNA/mapped/star_pe__{sample_name}_Aligned.out.bam"),
        count_file = temp("RNA/mapped/star_pe__{sample_name}_ReadsPerGene.out.tab"),
        metrics_map = "RNA/reports/star_pe__{sample_name}.txt"
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
        if [[ "{params.file_order}" == "rampage" ]]; then
            printf "Input file order for RAMPAGE (R2 R1)\n"
            input='"{input.fastq2}" "{input.fastq1}"'
        else
            printf "Input file order for RNAseq (R1 R2)\n"
            input='"{input.fastq1}" "{input.fastq2}"'
        fi
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
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
        bamfile = temp("RNA/mapped/star_se__{sample_name}_Aligned.out.bam"),
        count_file = temp("RNA/mapped/star_se__{sample_name}_ReadsPerGene.out.tab"),
        metrics_map = "RNA/reports/star_se__{sample_name}.txt"
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
        metrics_flag = "RNA/reports/flagstat_pe__{sample_name}.txt"
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
        metrics_flag = "RNA/reports/flagstat_se__{sample_name}.txt",
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
        bamfile = lambda wildcards: assign_mapping_paired(wildcards, "filter_rna", "sorted_file"),
        countfile = lambda wildcards: assign_mapping_paired(wildcards, "STAR_map", "count_file")
    output:
        bam_file = "RNA/mapped/final__{sample_name}.bam",
        count_file = "RNA/DEG/counts__{sample_name}.tab",
        touch = "RNA/chkpts/map_rna__{sample_name}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        mv {input.bamfile} {output.bam_file}
        mv {input.bamfile}.bai {output.bam_file}.bai
        mv {input.countfile} {output.count_file}
        touch {output.touch} 
        """

rule merging_rna_replicates:
    input:
        bamfiles = lambda wildcards: [ f"RNA/mapped/final__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.bam" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        mergefile = "RNA/mapped/merged__{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bam"
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
        bamfile = lambda wildcards: f"RNA/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'final'}__{wildcards.sample_name}.bam",
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
    log:
        temp(return_log_rna("{sample_name}", "making_bigiwig", ""))
    conda: CONDA_ENV
    threads: config["resources"]["filter_rna"]["threads"]
    resources:
        mem=config["resources"]["filter_rna"]["mem"],
        tmp=config["resources"]["filter_rna"]["tmp"]
    shell:
        """
        {{
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --runThreadN {threads} --inputBAMfile "{input.bamfile}" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "RNA/tracks/bg_{params.sample_name}_"
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
        }} 2>&1 | tee -a "{log}"
        """

rule prep_files_for_DEGs:
    input: 
        lambda wildcards: define_RNA_input_for_degs(wildcards.ref_genome)
    output:
        rna_samples = "RNA/DEG/samples__{analysis_name}__{ref_genome}.txt",
        rna_counts = "RNA/DEG/counts__{analysis_name}__{ref_genome}.txt"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_rna("{ref_genome}", "prep_for_DEGs", "{analysis_name}"))
    threads: config["resources"]["rna_degs"]["threads"]
    resources:
        mem=config["resources"]["rna_degs"]["mem"],
        tmp=config["resources"]["rna_degs"]["tmp"]
    run:
        filtered_samples = samples[ (samples['data_type'] == 'RNAseq') & (samples['ref_genome'] == {params.ref_genome}) ].copy()
        filtered_samples['Sample'] = filtered_samples['line'] + "__" + filtered_samples['tissue']
        filtered_samples['Replicate'] = filtered_samples['Sample'] + "__" + filtered_samples['replicate'].astype(str)
        
        RNA_samples = filtered_samples[['Replicate','Sample']].drop_duplicates()    
        RNA_samples = RNA_samples.sort_values(by=['Sample', 'Replicate'],ascending=[True, True]).reset_index(drop=True)
        RNA_samples['Color'] = pd.factorize(RNA_samples['Sample'])[0] + 1

        RNA_samples.to_csv(output.rna_samples, sep="\t", index=False)
        
        RNA_counts = None
        replicates = filtered_samples[['sample_name', 'Replicate']].drop_duplicates()
        for sname, rep in replicates.values:
            file_path = f"RNA/DEG/counts__{sname}.tab"
            temp = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1])
            temp.columns = ['GeneID', rep]

            if RNA_counts is None:
                RNA_counts = temp
            else:
                RNA_counts = pd.merge(RNA_counts, temp, on='GeneID', how='outer')
            
        replicate_order = RNA_samples['Replicate'].tolist()
        column_order = ['GeneID'] + replicate_order
        RNA_counts = RNA_counts[column_order]
        RNA_counts.to_csv(output.rna_counts, sep="\t", index=False)
    
rule call_all_DEGs:
    input:
        samples = "RNA/DEG/samples__{analysis_name}__{ref_genome}.txt",
        counts = "RNA/DEG/counts__{analysis_name}__{ref_genome}.txt"
    output:
        rdata = "RNA/DEG/ReadyToPlot_{analysis_name}__{ref_genome}.RData",
        touch = "RNA/chkpts/calling_DEGs__{analysis_name}__{ref_genome}.done"
    params:
        script = os.path.join(REPO_FOLDER,"scripts/R_call_DEGs.R"),
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        region_file = "combined/tracks/{ref_genome}__all_genes.bed"
    log:
        temp(return_log_rna("{ref_genome}", "call_DEGs", "{analysis_name}"))
    conda: os.path.join(REPO_FOLDER,"envs/call_degs.yaml")
    threads: config["resources"]["rna_degs"]["threads"]
    resources:
        mem=config["resources"]["rna_degs"]["mem"],
        tmp=config["resources"]["rna_degs"]["tmp"]
    shell:
        """
        printf "running edgeR for all samples in {params.ref_genome}\n"
        Rscript "{params.script}" "{input.counts}" "{input.samples}" "{params.analysis_name}" "{params.ref_genome}" "{params.region_file}"
        touch {output.touch}
        """

rule gather_gene_expression_rpkm:
    input:
        samples = "RNA/DEG/samples__{analysis_name}__{ref_genome}.txt",
        counts = "RNA/DEG/counts__{analysis_name}__{ref_genome}.txt"
    output:
        touch = "RNA/chkpts/gene_expression__{analysis_name}__{ref_genome}.done",
        rpkm = "RNA/DEG/genes_rpkm__{analysis_name}__{ref_genome}.txt"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        region_file = "combined/tracks/{ref_genome}__all_genes.bed"
    log:
        temp(return_log_rna("{ref_genome}", "gene_expression", "{analysis_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["region_file"]["threads"]
    resources:
        mem=config["resources"]["region_file"]["mem"],
        tmp=config["resources"]["region_file"]["tmp"]
    shell:
        """
        {{
        printf "Gathering gene expression levels for samples from {params.analysis_name} mapping to {params.ref_genome}\n"
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule all_rna:
    input:
        setup = "chkpts/directories_setup.done",
        final = lambda wildcards: define_final_rna_output(wildcards.ref_genome)
    output:
        touch = "RNA/chkpts/RNA_analysis__{analysis_name}__{ref_genome}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        touch {output.touch}
        """        
