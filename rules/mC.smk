# function to access logs more easily
def return_log_mc(sample_name, step, paired):
    logpath=os.path.join(REPO_FOLDER,"mC","logs",f"tmp__{sample_name}__{step}__{paired}.log")
    return logpath if config.get("debug_keep_logs", False) else temp(logpath)
    
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/mc.yaml")

def define_final_mC_output(ref_genome):
    qc_option = config["QC_option"]
    final_files = []
    qc_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'mC') & (samples['ref_genome'] == ref_genome) ]
    
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        final_files.append(f"mC/chkpts/bigwig__{sname}.done")
        if paired == "PE":
            final_files.append(f"mC/reports/final_reports_pe__{sname}.html")
            qc_files.append(f"mC/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
            qc_files.append(f"mC/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            qc_files.append(f"mC/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
            qc_files.append(f"mC/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
        else:
            final_files.append(f"mC/reports/final_reports_se__{sname}.html")
            qc_files.append(f"mC/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file
            qc_files.append(f"mC/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
        
    if qc_option == "all":
        return final_files + qc_files
    else:
        return final_files

rule make_bismark_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa"
    output:
        indices = directory("genomes/{ref_genome}/Bisulfite_Genome")
    params:
        limthreads = lambda wildcards, threads: threads // 2
    log:
        os.path.join(REPO_FOLDER,"logs","bismark_index_{ref_genome}.log")
    conda: CONDA_ENV
    threads: config["resources"]["bismark_indices"]["threads"]
    resources:
        mem=config["resources"]["bismark_indices"]["mem"],
        tmp=config["resources"]["bismark_indices"]["tmp"]
    shell:
        """
        {{
        printf "\nBuilding bismark index directory for {wildcards.ref_genome}\n"
        bismark_genome_preparation --parallel {params.limthreads} --bowtie2 --genomic_composition genomes/{wildcards.ref_genome}
        }} 2>&1 | tee -a "{log}"
        """
        
rule bismark_map_pe:
    input:
        fastq1 = "mC/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "mC/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/Bisulfite_Genome"
    output:
        temp_bamfile = temp("mC/mapped/{sample_name}/trim__{sample_name}__R1_bismark_bt2_pe.bam"),
        bamfile = "mC/mapped/{sample_name}/PE__{sample_name}.deduplicated.bam",
        cx_report = temp("mC/methylcall/PE__{sample_name}.deduplicated.CX_report.txt.gz"),
        metrics_alignement = temp("mC/mapped/{sample_name}/trim__{params.sample_name}__R1_bismark_bt2_PE_report.txt"),
        metrics_dedup = temp("mC/mapped/{sample_name}/trim__{params.sample_name}__R1_bismark_bt2_pe.deduplication_report.txt")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        mapping = config["mC_mapping"][config['mC_method']]['map_pe'],
        process = config["mC_mapping"][config['mC_method']]['process'],
        prefix = lambda wildcards: f"mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: threads // 4
    log:
        temp(return_log_mc("{sample_name}", "mapping", "PE"))
    conda: CONDA_ENV
    threads: config["resources"]["bismark_map"]["threads"]
    resources:
        mem=config["resources"]["bismark_map"]["mem"],
        tmp=config["resources"]["bismark_map"]["tmp"]
    shell:
        """
        {{
        printf "\nAligning {params.sample_name} with bismark/bowtie2\n"
        bismark --genome genomes/{params.ref_genome} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --gzip --nucleotide_coverage -1 {input.fastq1} -2 {input.fastq2}
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -p --output_dir {params.prefix}/ -o "PE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -p --comprehensive -o mC/methylcall/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder genomes/{params.ref_genome} {output.bamfile}
        rm -f mC/methylcall/C*context_PE__{params.sample_name}*
        rm -f mC/methylcall/PE__{params.sample_name}*bismark.cov*
        }} 2>&1 | tee -a "{log}"
        """

rule bismark_map_se:
    input:
        fastq0 = "mC/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/Bisulfite_Genome"
    output:
        temp_bamfile = temp("mC/mapped/{sample_name}/trim__{sample_name}__R0_bismark_bt2.bam"),
        bamfile = "mC/mapped/{sample_name}/SE__{sample_name}.deduplicated.bam",
        cx_report = temp("mC/methylcall/SE__{sample_name}.deduplicated.CX_report.txt.gz"),
        metrics_alignement = temp("mC/mapped/{sample_name}/trim__{params.sample_name}__R0_bismark_bt2_SE_report.txt"),
        metrics_dedup = temp("mC/mapped/{sample_name}/trim__{params.sample_name}__R0_bismark_bt2.deduplication_report.txt")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        mapping = config["mC_mapping"][config['mC_method']]['map_se'],
        process = config["mC_mapping"][config['mC_method']]['process'],
        prefix = lambda wildcards: f"mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: threads // 4
    log:
        temp(return_log_mc("{sample_name}", "mapping", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["bismark_map"]["threads"]
    resources:
        mem=config["resources"]["bismark_map"]["mem"],
        tmp=config["resources"]["bismark_map"]["tmp"]
    shell:
        """
        {{
        printf "\nAligning {params.sample_name} with bismark/bowtie2\n"
        bismark --genome genomes/{params.ref_genome} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --gzip --nucleotide_coverage {input.fastq0}
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -s --output_dir {params.prefix}/ -o "SE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -s --comprehensive -o mC/methylcall/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder genomes/{params.ref_genome} {output.bamfile}
        rm -f mC/methylcall/C*context_SE__{params.sample_name}*
        rm -f mC/methylcall/SE__{params.sample_name}*bismark.cov*
        }} 2>&1 | tee -a "{log}"
        """
        
rule make_mc_stats_pe:
    input:
        metrics_trim = "mC/reports/trim_pe__{sample_name}.txt",
        metrics_alignment = "mC/mapped/{sample_name}/trim__{params.sample_name}__R1_bismark_bt2_PE_report.txt",
        metrics_dedup = "mC/mapped/{sample_name}/trim__{params.sample_name}__R1_bismark_bt2_pe.deduplication_report.txt"
    output:
        stat_file = "mC/reports/summary_mC_PE_mapping_stats_{sample_name}.txt",
        reportfile = "mC/reports/final_reports_pe__{sample_name}.html"
    params:
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"mC/mapped/{wildcards.sample_name}"
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
        uniq=$(grep "Total count of deduplicated leftover sequences:" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{uniq}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_pe__{params.sample_name}.html" --dir mC/reports/ --alignment_report {input.metrics_alignment} --dedup_report {input.metrics_dedup} --splitting_report mC/methylcall/PE__{params.sample_name}.deduplicated_splitting_report.txt --mbias_report mC/methylcall/PE__{params.sample_name}.deduplicated.M-bias.txt --nucleotide_report {params.prefix}/trim__{params.sample_name}__R1_bismark_bt2_pe.nucleotide_stats.txt
        rm -f {input.logs}
        """
        
rule make_mc_stats_se:
    input:
        metrics_trim = "mC/reports/trim_se__{sample_name}.txt",
        metrics_alignment = "mC/mapped/{sample_name}/trim__{params.sample_name}__R0_bismark_bt2_SE_report.txt",
        metrics_dedup = "mC/mapped/{sample_name}/trim__{params.sample_name}__R0_bismark_bt2.deduplication_report.txt"
    output:
        stat_file = "mC/reports/summary_mC_SE_mapping_stats_{sample_name}.txt",
        reportfile = "mC/reports/final_reports_se__{sample_name}.html"
    params:
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"mC/mapped/{wildcards.sample_name}"
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
        uniq=$(grep "Total count of deduplicated leftover sequences:" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{uniq}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_se__{params.sample_name}.html" --dir mC/reports/ --alignment_report {input.metrics_alignment} --dedup_report {input.metrics_dedup} --splitting_report mC/methylcall/SE__{params.sample_name}.deduplicated_splitting_report.txt --mbias_report mC/methylcall/SE__{params.sample_name}.deduplicated.M-bias.txt --nucleotide_report {params.prefix}/trim__{params.sample_name}__R0_bismark_bt2.nucleotide_stats.txt
        rm -f {input.logs}
        """

rule pe_or_se_mc_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "bismark_map", "cx_report")
    output:
        cx_report = "mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        touch = "mC/chkpts/map__{sample_name}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        mv {input} {output.cx_report}
        touch {output.touch} 
        """

rule make_mc_bigwig_files:
    input:
        cx_report = "mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        bigwig = "mC/tracks/{sample_name}_CG.bw",
        touch = "mC/chkpts/bigwig__{sample_name}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        context = config['mC_context']
    log:
        temp(return_log_mc("{sample_name}", "bigwig", "both"))
    conda: CONDA_ENV
    threads: config["resources"]["mc_bigwig"]["threads"]
    resources:
        mem=config["resources"]["mc_bigwig"]["mem"],
        tmp=config["resources"]["mc_bigwig"]["tmp"]    
    shell:
        """
        {{
        if [ {params.context} == "All" ]; then
            zcat {input.cx_report} | awk -v OFS="\t" -v s={params.sample_name} '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "mC/methylcall/"s"_CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "mC/methylcall/"s"_CHG.bedGraph"; else print $1,$2-1,$2,$4/a*100 > "mC/methylcall/"s"_CG.bedGraph"}}'
            for strand in plus minus
            do
                case "${{strand}}" in 
                    plus)	sign="+";;
                    minus)	sign="-";;
                esac
                zcat {input.cx_report} | awk -v n=${{sign}} '$3==n' | awk -v OFS="\t" -v s={params.sample_name} -v d=${{strand}} '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "mC/methylcall/"s"_CHH_"d".bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "mC/methylcall/"s"_CHG_"d".bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 > "mC/methylcall/"s"_CG_"d".bedGraph"}}'
            done
            for context in CG CHG CHH
            do
                printf "\nMaking bigwig files of ${{context}} context for {params.sample_name}\n"
                LC_COLLATE=C sort -k1,1 -k2,2n mC/methylcall/{params.sample_name}_$){context}}.bedGraph > mC/methylcall/sorted_{params.sample_name}_${{context}}.bedGraph
                bedGraphToBigWig mC/methylcall/sorted_{params.sample_name}_${{context}}.bedGraph {input.chrom_sizes} mC/methylcall/{params.sample_name}_${{context}}.bw
                for strand in plus minus
                do
                    printf "\nMaking ${{strand}} strand bigwig files of ${{context}} context for {params.sample_name}\n"
                    LC_COLLATE=C sort -k1,1 -k2,2n mC/methylcall/{params.sample_name}_${{context}}_${{strand}}.bedGraph > mC/methylcall/sorted_{params.sample_name}_${{context}}_${{strand}}.bedGraph
                    bedGraphToBigWig mC/methylcall/sorted_{params.sample_name}_${{context}}_${{strand}}.bedGraph {input.chrom_sizes} mC/methylcall/{params.sample_name}_${{context}}_${{strand}}.bw
                done
            done
            rm -f mC/methylcall/*"{params.sample_name}"*bedGraph*
        elif [ {params.context} == "CG-only" ]; then
            printf "Script for CG-only not ready yet\n" ## To update for CG-only!
        else
            printf "Unknown sequence context selection! Check the config file and set 'mC_context' to either 'All' or 'CG-only'\n"
            exit 1
        fi
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule all_mC:
    input:
        lambda wildcards: define_final_mC_output(wildcards.ref_genome)
    output:
        touch = "mC/chkpts/mC_analysis__{ref_genome}.done"
    threads: 1
    resources:
        mem=32,
        tmp=32
    shell:
        """
        touch {output.touch}
        """        
