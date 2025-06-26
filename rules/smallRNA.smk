ADD rules resources to config and clusters
dispatch_srna_fastq: *default
filter_structural_rna: *heavy
shortstack_map: *heavy

# function to access logs more easily
def return_log_smallrna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"smallRNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def define_input_file_for_shortstack(sample_name):
    paired = get_sample_info_from_name(sample_name, samples, 'paired')
    if paired == "se":
        return "filtered__{sample_name}__R0" if config['structural_rna_depletion'] else "trim__{sample_name}__R0"

rule make_bt2_indices:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gff = "genomes/{ref_genome}/{ref_genome}.gff",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        indices = directory("genomes/{ref_genome}/bt2_index")
    log:
        temp(os.path.join(REPO_FOLDER,"logs","bowtie_index_{ref_genome}.log"))
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
    
rule filter_structural_rna:
    input:
        fastq = "sRNA/fastq/trim__{sample_name}__R0.fastq.gz",
        fasta = config['structural_rna_fafile']
    output:
        filtered_fastq = "sRNA/fastq/filtered__{sample_name}__R0.fastq.gz"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "SE"))
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
                mkdir genomes/structural_RNAs/{params.ref_genome}/bt2_index
                gunzip {input.fasta} | bowtie2-build --threads {threads} -c "genomes/structural_RNAs/{params.ref_genome}/bt2_index"
            else
                printf "Generating bowtie2 index for structural RNAs of {params.ref_genome} using file: {input.fasta}\n"
                mkdir genomes/structural_RNAs/{params.ref_genome}/bt2_index
                bowtie2-build --threads {threads} "{input.fasta}" "genomes/structural_RNAs/{params.ref_genome}/bt2_index"
            fi
        fi
        bowtie2 --very-sensitive -p {threads} -x genomes/structural_RNAs/{params.ref_genome}/bt2_index -U {input.fastq} | samtools view -@ {threads} -f 0x4 | samtools fastq -@ {threads} > {output.filtered_fastq}
        pigz -p {threads} {output.filtered_fastq}
        }} 2>&1 | tee -a "{log}"
        """

rule dispatch_srna_fastq:
    input:
        fastq = lambda wildcards: f"sRNA/fastq/{define_input_file_for_shortstack(wildcards.sample_name)}.fastq.gz"
    output:
        fastq_file = temp("sRNA/fastq/clean_(sample_name).fastq.gz")
    conda: CONDA_ENV
    threads: config["resources"]["dispatch_srna_fastq"]["threads"]
    resources:
        mem=config["resources"]["dispatch_srna_fastq"]["mem"],
        tmp=config["resources"]["dispatch_srna_fastq"]["tmp"]
    shell:
        """
        cp {input.fastq} {output.fastq_file}
        """

rule shortstack_map:
    input:
        fastq = "sRNA/fastq/clean_{sample_name}.fastq.gz",
        fasta = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa"
    output:
        count_file = "sRNA/mapped/{sample_name}/ShortStack_All.gff3",
        bam_file = temp("RNA/mapped/clean_{sample_name).bam")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        srna_params = config['srna_mapping_params']
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["shortstack_map"]["threads"]
    resources:
        mem=config["resources"]["shortstack_map"]["mem"],
        tmp=config["resources"]["shortstack_map"]["tmp"]
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with Shortstack version:\n"
        ShortStack --version
        ShortStack --readfile {input.fastq} --genomefile {input.fasta} --bowtie_cores {threads} --sort_mem {resources.mem} {params.srna_params} --outdir sRNA/mapped/{params.sample_name}
        }} 2>&1 | tee -a "{log}"
        """
        
rule process_srna_sample:
    input:
        bamfile = "RNA/mapped/clean_{sample_name).bam"
    output:
        count_file = "sRNA/mapped/{sample_name}/ShortStack_All.gff3",
        bam_file = temp("RNA/mapped/clean_{sample_name).bam")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        srna_min = config['srna_min_size'],
        srna_max = config['srna_max_size']
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["shortstack_map"]["threads"]
    resources:
        mem=config["resources"]["shortstack_map"]["mem"],
        tmp=config["resources"]["shortstack_map"]["tmp"]
    shell:
        """
        {{
        samtools index -@ {threads} {input.bamfile}
        printf "Filtering only small RNA sizes ({params.srna_min}nt to {params.srna_max}nt) for {params.sample_name}\n"
        for ((nt={params.srna_min}; nt<={params.srna_max}; nt++)); do
            samtools view -h {input.bamfile} | awk -v m={params.srna_min} -v n={params.srna_max} '(length($10) >= m && length($10) <= n) || $1 ~ /^@/' | samtools view -bS - > mapped/${name}/sized_${name}.bam
        #### Getting stats of size distribution
        printf "\nGetting trimmed stats for ${name}\n"
        zcat fastq/trimmed_${name}.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c | awk -v OFS="\t" -v n=${name} '{print n,"trimmed",$2,$1}' > reports/sizes_trimmed_${name}.txt
        printf "\nGetting filtered stats for ${name}\n"
        zcat fastq/filtered_${name}.fastq.gz | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c | awk -v OFS="\t" -v n=${name} '{print n,"filtered",$2,$1}' > reports/sizes_filtered_${name}.txt
        printf "\nGetting mapped stats for ${name}\n"
        samtools view mapped/${name}/filtered_${name}.bam | awk '$2==0 || $2==16 {print length($10)}' | sort -n | uniq -c | awk -v OFS="\t" -v n=${name} '{print n,"mapped",$2,$1}' > reports/sizes_mapped_${name}.txt
        }} 2>&1 | tee -a "{log}"
        """

		# # printf "Filtering ${line}_${tissue}_Rep1 shRNA files for reads =24nt\n"
		# # samtools view -@ ${threads} -h ${pathtobamshrna}/${line}_${tissue}_shRNA_Rep1/filtered_${line}_${tissue}_shRNA_Rep1.bam | awk 'length($10) == 24 || $1 ~ /^@/' | samtools view -bS - > ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_Rep1.bam
		# # printf "Filtering ${line}_${tissue}_Rep2 shRNA files for reads =24nt\n"
		# # samtools view -@ ${threads} -h ${pathtobamshrna}/${line}_${tissue}_shRNA_Rep2/filtered_${line}_${tissue}_shRNA_Rep2.bam | awk 'length($10) == 24 || $1 ~ /^@/' | samtools view -bS - > ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_Rep2.bam
		# # printf "Merging ${line}_${tissue} Replicates\n"
		# # samtools merge -@ ${threads} ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_merged.bam ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_Rep*.bam
		# # samtools sort -@ ${threads} -o ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam ${pathtobamshrna}/temp_${line}_${tissue}_24RNA_merged.bam
		# # rm -f ${pathtobamshrna}/temp_${line}_${tissue}_24RNA*.bam
		# # samtools index -@ ${threads} ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam
		# # printf "Getting stranded coverage for ${line}_${tissue}\n"
		# # bamCoverage --filterRNAstrand forward -bs 1 -p ${threads} --normalizeUsing CPM -b ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam -o ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_plus.bw
		# # bamCoverage --filterRNAstrand reverse -bs 1 -p ${threads} --normalizeUsing CPM -b ${pathtobamshrna}/${line}_${tissue}_24RNA_merged.bam -o ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_minus.bw
		# # printf "Merging strands for ${line}_${tissue} 24RNA\n"
		# # bigWigMerge ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_plus.bw ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_minus.bw ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bg
		# # sort -k1,1 -k2,2n ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bg > ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum_sorted.bg
		# # bedGraphToBigWig ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum_sorted.bg ${ref_dir}/chrom.sizes ${pathtobwshrna}/${line}_${tissue}_24RNA_merged_sum.bw
		# # rm -f ${pathtobwshrna}/${line}_${tissue}_24RNA*.bg