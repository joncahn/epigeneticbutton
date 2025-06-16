# function to access logs more easily
def return_log_mc(sample_name, step, paired):
    logpath=os.path.join(REPO_FOLDER,"mC","logs",f"tmp__{sample_name}__{step}__{paired}.log")
    return logpath if config.get("debug_keep_logs", False) else temp(logpath)
    
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/mc.yaml")

def parameters_for_mc(sample_name):
    temp = parse_sample_name(sample_name)['sample_type']
    options = {"WGBS", "Pico", "EMseq"}
    return temp if temp in options else "default"

def define_DMR_samples(sample_name):
    data_type = get_sample_info_from_name(sample_name, analysis_samples, 'data_type')
    line = get_sample_info_from_name(sample_name, analysis_samples, 'line')
    tissue = get_sample_info_from_name(sample_name, analysis_samples, 'tissue')
    sample_type = get_sample_info_from_name(sample_name, analysis_samples, 'sample_type')
    ref_genome = get_sample_info_from_name(sample_name, analysis_samples, 'ref_genome')
    replicates = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])
    
    return [ f"mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.deduplicated.CX_report.txt.gz"
                    for replicate in replicates ]

def define_final_mC_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    final_files = []
    dmr_files = []
    merged_files = []
    qc_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'mC') & (samples['ref_genome'] == ref_genome) ]
    
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        final_files.append(f"mC/chkpts/bigwig__{sname}.done")
        if paired == "PE":
            final_files.append(f"mC/reports/final_report_pe__{sname}.html")
            qc_files.append(f"mC/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
            qc_files.append(f"mC/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            qc_files.append(f"mC/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
            qc_files.append(f"mC/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
        else:
            final_files.append(f"mC/reports/final_report_se__{sname}.html")
            qc_files.append(f"mC/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file
            qc_files.append(f"mC/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
    
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'mC') & (analysis_samples['ref_genome'] == ref_genome) ]
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            merged_files.append(f"mC/chkpts/bigwig__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.done") # merged bigwig files
    
    for a, b in combinations(filtered_analysis_samples.itertuples(index=False), 2):
        a_dict = a._asdict()
        b_dict = b._asdict()
        sample1 = sample_name_str(a_dict, 'analysis')
        sample2 = sample_name_str(b_dict, 'analysis')
        dmr_files.append(f"mC/DMRs/summary__{sample1}__vs__{sample2}__dmrs.txt")
    
    if qc_option == "all":
        results = final_files + qc_files
    else:
        results = final_files
    
    if analysis:
        results += dmr_files + merged_files

    return results

rule make_bismark_indices:
    input:
        fasta = "genomes/{ref_genome}/temp_{ref_genome}.fa"
    output:
        indices = directory("genomes/{ref_genome}/Bisulfite_Genome")
    params:
        limthreads = lambda wildcards, threads: max(1, threads // 2)
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
        if [[ {params.limthreads} -gt 1 ]]; then
            bismark_genome_preparation --parallel {params.limthreads} --bowtie2 --genomic_composition genomes/{wildcards.ref_genome}
        else
            bismark_genome_preparation --bowtie2 --genomic_composition genomes/{wildcards.ref_genome}
        fi
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
        metrics_alignement = temp("mC/mapped/{sample_name}/trim__{sample_name}__R1_bismark_bt2_PE_report.txt"),
        metrics_dedup = temp("mC/mapped/{sample_name}/PE__{sample_name}.deduplication_report.txt")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome_path = lambda wildcards: os.path.join(REPO_FOLDER,"genomes",parse_sample_name(wildcards.sample_name)['ref_genome']),
        mapping = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['map_pe'],
        process = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['process_pe'],
        prefix = lambda wildcards: f"mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: max(1, threads // 3)
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
        bismark --genome {params.ref_genome_path} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --gzip --nucleotide_coverage -1 {input.fastq1} -2 {input.fastq2}
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -p --output_dir {params.prefix}/ -o "PE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -p --comprehensive -o mC/methylcall/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder {params.ref_genome_path} {output.bamfile}
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
        metrics_map = temp("mC/mapped/{sample_name}/trim__{sample_name}__R0_bismark_bt2_SE_report.txt"),
        metrics_dedup = temp("mC/mapped/{sample_name}/SE__{sample_name}.deduplication_report.txt")
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome_path = lambda wildcards: os.path.join(REPO_FOLDER,"genomes",parse_sample_name(wildcards.sample_name)['ref_genome']),
        mapping = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['map_se'],
        process = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['process_se'],
        prefix = lambda wildcards: f"mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: max(1, threads // 3)
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
        bismark --genome {params.ref_genome_path} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --gzip --nucleotide_coverage {input.fastq0}
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -s --output_dir {params.prefix} -o "SE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -s --comprehensive -o mC/methylcall/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder {params.ref_genome_path} {output.bamfile}
        rm -f mC/methylcall/C*context_SE__{params.sample_name}*
        rm -f mC/methylcall/SE__{params.sample_name}*bismark.cov*
        }} 2>&1 | tee -a "{log}"
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
        
rule make_mc_stats_pe:
    input:
        metrics_trim = "mC/reports/trim_pe__{sample_name}.txt",
        metrics_map = "mC/mapped/{sample_name}/trim__{sample_name}__R1_bismark_bt2_PE_report.txt",
        metrics_dedup = "mC/mapped/{sample_name}/PE__{sample_name}.deduplication_report.txt",
        cx_report = "mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        stat_file = "mC/reports/summary_mC_PE_mapping_stats_{sample_name}.txt",
        reportfile = "mC/reports/final_report_pe__{sample_name}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
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
        filt=$(grep "Sequences analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Sequences did not map uniquely" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Number of alignments with a unique best hit" "{input.metrics_map}" | awk '{{print $NF}}')
        uniq=$(grep "Total count of deduplicated leftover sequences" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > {output.stat_file}
        ## Can change the name of the plastid chromosome to calculate non-conversion rate
        zcat {input.cx_report} | awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v s={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v x=${{tot}} -v y=${{filt}} -v z=${{allmap}} -v u=${{uniq}} '{{a+=1; b=$4+$5; i+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {{m+=$4; n+=b;}}; if (b>0) {{c+=1; d+=b;}}; if (b>2) e+=1}} END {{if (n>0) {{o=m/n*100;}} else o="NA"; print l,t,s,r,g,x,y" ("y/x*100"%)",z" ("z/x*100"%)",u" ("u/x*100"%)",c/a*100,e/a*100,i/a,d/c,o}}' >> "{output.stat_file}"

        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_pe__{params.sample_name}.html" --dir mC/reports/ --alignment_report {input.metrics_map} --dedup_report {input.metrics_dedup} --splitting_report mC/methylcall/PE__{params.sample_name}.deduplicated_splitting_report.txt --mbias_report mC/methylcall/PE__{params.sample_name}.deduplicated.M-bias.txt --nucleotide_report {params.prefix}/trim__{params.sample_name}__R1_bismark_bt2_pe.nucleotide_stats.txt
        cp mC/methylcall/PE__"{params.sample_name}"*.txt mC/reports/
        cp {params.prefix}/trim__"{params.sample_name}"*.txt mC/reports/
        """
        
rule make_mc_stats_se:
    input:
        metrics_trim = "mC/reports/trim_se__{sample_name}.txt",
        metrics_map = "mC/mapped/{sample_name}/trim__{sample_name}__R0_bismark_bt2_SE_report.txt",
        metrics_dedup = "mC/mapped/{sample_name}/SE__{sample_name}.deduplication_report.txt",
        cx_report = "mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        stat_file = "mC/reports/summary_mC_SE_mapping_stats_{sample_name}.txt",
        reportfile = "mC/reports/final_report_se__{sample_name}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
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
        tot=$(grep "Total reads processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        filt=$(grep "Sequences analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Sequences did not map uniquely" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Number of alignments with a unique best hit" "{input.metrics_map}" | awk '{{print $NF}}')
        uniq=$(grep "Total count of deduplicated leftover sequences" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > {output.stat_file}
        ## Can change the name of the plastid chromosome to calculate non-conversion rate
        zcat {input.cx_report} | awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v s={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v x=${{tot}} -v y=${{filt}} -v z=${{allmap}} -v u=${{uniq}} '{{a+=1; b=$4+$5; i+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {{m+=$4; n+=b;}}; if (b>0) {{c+=1; d+=b;}}; if (b>2) e+=1}} END {{if (n>0) {{o=m/n*100;}} else o="NA"; print l,t,s,r,g,x,y" ("y/x*100"%)",z" ("z/x*100"%)",u" ("u/x*100"%)",c/a*100,e/a*100,i/a,d/c,o}}' >> "{output.stat_file}"

        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_se__{params.sample_name}.html" --dir mC/reports/ --alignment_report {input.metrics_map} --dedup_report {input.metrics_dedup} --splitting_report mC/methylcall/SE__{params.sample_name}.deduplicated_splitting_report.txt --mbias_report mC/methylcall/SE__{params.sample_name}.deduplicated.M-bias.txt --nucleotide_report {params.prefix}/trim__{params.sample_name}__R0_bismark_bt2.nucleotide_stats.txt
        mv mC/methylcall/SE__"{params.sample_name}"*.txt mC/reports/
        mv {params.prefix}/trim__"{params.sample_name}"*.txt mC/reports/
        """

rule merging_mc_replicates:
    input:
        report_files = lambda wildcards: [ f"mC/methylcall/{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.deduplicated.CX_report.txt.gz" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        bedfile = temp("mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bed"),
        tempmergefile = temp("mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.deduplicated.CX_report.txt"),
        mergefile = temp("mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.deduplicated.CX_report.txt.gz")
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis')
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging_reps", ""))
    conda: CONDA_ENV
    threads: config["resources"]["merging_mc_replicates"]["threads"]
    resources:
        mem=config["resources"]["merging_mc_replicates"]["mem"],
        tmp=config["resources"]["merging_mc_replicates"]["tmp"]
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
        zcat {input.report_files} | sort -k1,1 -k2,2n | awk -v OFS="\t" '{{print $1,$2-1,$2,$3,$4,$5,$6,$7}}' > {output.bedfile}
		bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i {output.bedfile} | awk -v OFS="\t" '{{print $1,$3,$4,$5,$6,$7,$8}}' > {output.tempmergefile}
        pigz -p {threads} "{output.tempmergefile}" -c > "{output.mergefile}"
        }} 2>&1 | tee -a "{log}"
        """    

rule make_mc_bigwig_files:
    input:
        cx_report = "mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        bigwig = "mC/tracks/{sample_name}__CG.bw",
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
        if [[ "{params.context}" == "all" ]]; then
            zcat {input.cx_report} | awk -v OFS="\t" -v s={params.sample_name} '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "mC/tracks/"s"__CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "mC/tracks/"s"__CHG.bedGraph"; else print $1,$2-1,$2,$4/a*100 > "mC/tracks/"s"__CG.bedGraph"}}'
            for strand in plus minus
            do
                case "${{strand}}" in 
                    plus)	sign="+";;
                    minus)	sign="-";;
                esac
                zcat {input.cx_report} | awk -v n=${{sign}} '$3==n' | awk -v OFS="\t" -v s={params.sample_name} -v d=${{strand}} '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "mC/tracks/"s"__CHH__"d".bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "mC/tracks/"s"__CHG__"d".bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 > "mC/tracks/"s"__CG__"d".bedGraph"}}'
            done
            for context in CG CHG CHH
            do
                printf "\nMaking bigwig files of ${{context}} context for {params.sample_name}\n"
                LC_COLLATE=C sort -k1,1 -k2,2n mC/tracks/{params.sample_name}__${{context}}.bedGraph > mC/tracks/sorted__{params.sample_name}__${{context}}.bedGraph
                bedGraphToBigWig mC/tracks/sorted__{params.sample_name}__${{context}}.bedGraph {input.chrom_sizes} mC/tracks/{params.sample_name}__${{context}}.bw
                for strand in plus minus
                do
                    printf "\nMaking ${{strand}} strand bigwig files of ${{context}} context for {params.sample_name}\n"
                    LC_COLLATE=C sort -k1,1 -k2,2n mC/tracks/{params.sample_name}__${{context}}__${{strand}}.bedGraph > mC/tracks/sorted__{params.sample_name}__${{context}}__${{strand}}.bedGraph
                    bedGraphToBigWig mC/tracks/sorted__{params.sample_name}__${{context}}__${{strand}}.bedGraph {input.chrom_sizes} mC/tracks/{params.sample_name}__${{context}}__${{strand}}.bw
                done
            done
            rm -f mC/tracks/*"{params.sample_name}"*bedGraph*
        elif [[ "{params.context}" == "CG-only" ]]; then
            printf "Script for CG-only not ready yet\n" ## To update for CG-only!
        else
            printf "Unknown sequence context selection! Check the config file and set 'mC_context' to either 'all' or 'CG-only'\n"
            exit 1
        fi
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule call_DMRs_pairwise:
    input:
        sample1 = lambda wildcards: define_DMR_samples(wildcards.sample1),
        sample2 = lambda wildcards: define_DMR_samples(wildcards.sample2),
        chrom_sizes = lambda wildcards: f"genomes/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}/chrom.sizes"
    output:
        dmr_summary = "mC/DMRs/summary__{sample1}__vs__{sample2}__DMRs.txt"
    params:
        script = os.path.join(REPO_FOLDER,"scripts/R_call_DMRs.R"),
        context = config['mC_context'],
        sample1 = lambda wildcards: wildcards.sample1,
        sample2 = lambda wildcards: wildcards.sample2,
        nb_sample1 = lambda wildcards: len(define_DMR_samples(wildcards.sample1)),
        nb_sample2 = lambda wildcards: len(define_DMR_samples(wildcards.sample2))
    log:
        temp(return_log_mc("{sample1}__vs__{sample2}", "DMRs", ""))
    conda: os.path.join(REPO_FOLDER,"envs/call_dmrs.yaml")
    threads: config["resources"]["call_dmrs"]["threads"]
    resources:
        mem=config["resources"]["call_dmrs"]["mem"],
        tmp=config["resources"]["call_dmrs"]["tmp"]
    shell:
        """
        printf "running DMRcaller for {params.sample1} vs {params.sample2}\n"
        Rscript "{params.script}" "{threads}" "{input.chrom_sizes}" "{params.context}" "{params.sample1}" "{params.sample2}" "{params.nb_sample1}" "{params.nb_sample2}" {input.sample1} {input.sample2}
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
