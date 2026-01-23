CONDA_ENV_MC=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_mc.yaml")

def return_log_mc(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"results","mC","logs",f"tmp__{sample_name}__{step}__{paired}.log")
     
def parameters_for_mc(sample_name):
    temp = parse_sample_name(sample_name)['sample_type']
    options = {"WGBS", "Pico", "EMseq", "ONT", "bedMethyl"}
    return temp if temp in options else "default"

def is_ont_sample(sample_name):
    """Check if a sample uses ONT direct methylation workflow."""
    return parse_sample_name(sample_name)['sample_type'] in ["ONT", "bedMethyl"]

def get_ont_input_type(sample_name):
    """Return the input type for ONT samples: 'bedMethyl' or 'modBAM'."""
    sample_type = parse_sample_name(sample_name)['sample_type']
    return "bedMethyl" if sample_type == "bedMethyl" else "modBAM"

def define_cx_report_input(wildcards):
    name = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}"
    if wildcards.replicate == "merged":
        return f"results/mC/methylcall/{name}.merged.CX_report.txt.gz"
    else:
        return f"results/mC/methylcall/{name}.deduplicated.CX_report.txt.gz"

def define_DMR_samples(sample_name):
    data_type = get_sample_info_from_name(sample_name, analysis_samples, 'data_type')
    line = get_sample_info_from_name(sample_name, analysis_samples, 'line')
    tissue = get_sample_info_from_name(sample_name, analysis_samples, 'tissue')
    sample_type = get_sample_info_from_name(sample_name, analysis_samples, 'sample_type')
    ref_genome = get_sample_info_from_name(sample_name, analysis_samples, 'ref_genome')
    replicates = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])
    
    return [ f"results/mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.deduplicated.CX_report.txt.gz"
                    for replicate in replicates ]
                    
def script_DMRs():
    script_dmrs = config['custom_script_dmrs']
    default = os.path.join(REPO_FOLDER,"workflow","scripts","R_call_DMRs.R")
    custom = os.path.join(REPO_FOLDER,"workflow","scripts","R_call_DMRs_custom.R")
    return custom if script_dmrs else default

def define_final_mC_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    trimmed_fastqs = config['trimmed_fastqs']
    map_files = []
    dmr_files = []
    bigwig_files = []
    qc_files = []
    ont_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'mC') & (samples['ref_genome'] == ref_genome) ].copy()

    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        sample_type = parse_sample_name(sname)['sample_type']

        # ONT samples use different workflow
        if sample_type in ["ONT", "bedMethyl"]:
            bigwig_files.append(f"results/mC/chkpts/bigwig__{sname}.done")
            ont_files.append(f"results/mC/ont/summary__{sname}.txt")  # modkit summary
        else:
            # Bismark workflow
            bigwig_files.append(f"results/mC/chkpts/bigwig__{sname}.done")
            if paired == "PE":
                map_files.append(f"results/mC/reports/final_report_pe__{sname}.html")
                qc_files.append(f"results/mC/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
                qc_files.append(f"results/mC/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
                if not trimmed_fastqs:
                    qc_files.append(f"results/mC/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
                    qc_files.append(f"results/mC/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            else:
                map_files.append(f"results/mC/reports/final_report_se__{sname}.html")
                qc_files.append(f"results/mC/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
                if not trimmed_fastqs:
                    qc_files.append(f"results/mC/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file

    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'mC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        if len(analysis_to_replicates[(row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome)]) >= 2:
            bigwig_files.append(f"results/mC/chkpts/bigwig__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.done") # merged bigwig files

    for a, b in combinations(filtered_analysis_samples.itertuples(index=False), 2):
        a_dict = a._asdict()
        b_dict = b._asdict()
        sample1 = sample_name_str(a_dict, 'analysis')
        sample2 = sample_name_str(b_dict, 'analysis')
        dmr_files.append(f"results/mC/DMRs/summary__{sample1}__vs__{sample2}__DMRs.txt")

    results = map_files + bigwig_files + ont_files

    if qc_option == "all":
        results += qc_files

    if analysis:
        results += dmr_files

    return results

rule make_bismark_indices:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    output:
        indices = directory("genomes/{ref_genome}/Bisulfite_Genome")
    params:
        limthreads = lambda wildcards, threads: max(1, threads // 2)
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","bismark_index_{ref_genome}.log"))
    conda: CONDA_ENV_MC
    threads: config["resources"]["make_bismark_indices"]["threads"]
    resources:
        mem_mb=config["resources"]["make_bismark_indices"]["mem_mb"],
        tmp_mb=config["resources"]["make_bismark_indices"]["tmp_mb"],
        qos=config["resources"]["make_bismark_indices"]["qos"]
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
        fastq1 = "results/mC/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "results/mC/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/Bisulfite_Genome"
    output:
        temp_bamfile = temp("results/mC/mapped/{sample_name}/trim__{sample_name}__R1_bismark_bt2_pe.bam"),
        bamfile = "results/mC/mapped/{sample_name}/PE__{sample_name}.deduplicated.bam",
        cx_report = temp("results/mC/mapped/PE__{sample_name}.deduplicated.CX_report.txt.gz"),
        metrics_alignement = temp("results/mC/mapped/{sample_name}/trim__{sample_name}__R1_bismark_bt2_PE_report.txt"),
        metrics_dedup = temp("results/mC/mapped/{sample_name}/PE__{sample_name}.deduplication_report.txt")
    wildcard_constraints:
        sample_name = r"(?!.*__(ONT|bedMethyl)__).*"  # Exclude ONT/bedMethyl samples
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome_path = lambda wildcards: os.path.join(REPO_FOLDER,"genomes",parse_sample_name(wildcards.sample_name)['ref_genome']),
        mapping = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['map_pe'],
        process = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['process_pe'],
        prefix = lambda wildcards: f"results/mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: max(1, threads // 3)
    log:
        temp(return_log_mc("{sample_name}", "mapping", "PE"))
    conda: CONDA_ENV_MC
    threads: config["resources"]["bismark_map_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["bismark_map_pe"]["mem_mb"],
        tmp_mb=config["resources"]["bismark_map_pe"]["tmp_mb"],
        qos=config["resources"]["bismark_map_pe"]["qos"]
    shell:
        """
        {{
        printf "\nAligning {params.sample_name} with bismark/bowtie2\n"
        bismark --genome {params.ref_genome_path} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --gzip --nucleotide_coverage -1 {input.fastq1} -2 {input.fastq2}
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -p --output_dir {params.prefix}/ -o "PE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -p --comprehensive -o results/mC/mapped/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder {params.ref_genome_path} {output.bamfile}
        rm -f results/mC/mapped/C*context_PE__{params.sample_name}*
        rm -f results/mC/mapped/PE__{params.sample_name}*bismark.cov*
        rm -f results/mC/mapped/PE__{params.sample_name}*bedGraph*
        }} 2>&1 | tee -a "{log}"
        """

rule bismark_map_se:
    input:
        fastq0 = "results/mC/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/Bisulfite_Genome"
    output:
        temp_bamfile = temp("results/mC/mapped/{sample_name}/trim__{sample_name}__R0_bismark_bt2.bam"),
        bamfile = "results/mC/mapped/{sample_name}/SE__{sample_name}.deduplicated.bam",
        cx_report = temp("results/mC/mapped/SE__{sample_name}.deduplicated.CX_report.txt.gz"),
        metrics_map = temp("results/mC/mapped/{sample_name}/trim__{sample_name}__R0_bismark_bt2_SE_report.txt"),
        metrics_dedup = temp("results/mC/mapped/{sample_name}/SE__{sample_name}.deduplication_report.txt")
    wildcard_constraints:
        sample_name = r"(?!.*__(ONT|bedMethyl)__).*"  # Exclude ONT/bedMethyl samples
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome_path = lambda wildcards: os.path.join(REPO_FOLDER,"genomes",parse_sample_name(wildcards.sample_name)['ref_genome']),
        mapping = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['map_se'],
        process = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['process_se'],
        prefix = lambda wildcards: f"results/mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: max(1, threads // 3)
    log:
        temp(return_log_mc("{sample_name}", "mapping", "SE"))
    conda: CONDA_ENV_MC
    threads: config["resources"]["bismark_map_se"]["threads"]
    resources:
        mem_mb=config["resources"]["bismark_map_se"]["mem_mb"],
        tmp_mb=config["resources"]["bismark_map_se"]["tmp_mb"],
        qos=config["resources"]["bismark_map_se"]["qos"]
    shell:
        """
        {{
        printf "\nAligning {params.sample_name} with bismark/bowtie2\n"
        bismark --genome {params.ref_genome_path} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --gzip --nucleotide_coverage {input.fastq0}
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -s --output_dir {params.prefix} -o "SE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -s --comprehensive -o results/mC/mapped/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder {params.ref_genome_path} {output.bamfile}
        rm -f results/mC/mapped/C*context_SE__{params.sample_name}*
        rm -f results/mC/mapped/SE__{params.sample_name}*bismark.cov*
        rm -f results/mC/mapped/SE__{params.sample_name}*bedGraph*
        }} 2>&1 | tee -a "{log}"
        """

rule pe_or_se_mc_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "bismark_map", "cx_report")
    output:
        cx_report = "results/mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        touch = "results/mC/chkpts/map_mC__{sample_name}.done"
    wildcard_constraints:
        sample_name = r"(?!.*__(ONT|bedMethyl)__).*"  # Exclude ONT/bedMethyl samples
    localrule: True
    shell:
        """
        mv {input} {output.cx_report}
        touch {output.touch} 
        """
        
rule make_mc_stats_pe:
    input:
        metrics_trim = "results/mC/reports/trim_pe__{sample_name}.txt",
        metrics_map = "results/mC/mapped/{sample_name}/trim__{sample_name}__R1_bismark_bt2_PE_report.txt",
        metrics_dedup = "results/mC/mapped/{sample_name}/PE__{sample_name}.deduplication_report.txt",
        cx_report = "results/mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        stat_file = "results/mC/reports/summary_mC_PE_mapping_stats_{sample_name}.txt",
        reportfile = "results/mC/reports/final_report_pe__{sample_name}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"results/mC/mapped/{wildcards.sample_name}",
        trimmed_fastq = config['trimmed_fastqs']
    log:
        temp(return_log_mc("{sample_name}", "making_stats", "PE"))
    conda: CONDA_ENV_MC
    threads: config["resources"]["make_mc_stats_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["make_mc_stats_pe"]["mem_mb"],
        tmp_mb=config["resources"]["make_mc_stats_pe"]["tmp_mb"],
        qos=config["resources"]["make_mc_stats_pe"]["qos"]
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        else
            tot=$(grep "Sequence pairs analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        fi
        filt=$(grep "Sequence pairs analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Sequence pairs did not map uniquely" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Number of paired-end alignments with a unique best hit" "{input.metrics_map}" | awk '{{print $NF}}')
        uniq=$(grep "Total count of deduplicated leftover sequences" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > {output.stat_file}
        ## Can change the name of the plastid chromosome to calculate non-conversion rate
        zcat {input.cx_report} | awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v s={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v x=${{tot}} -v y=${{filt}} -v z=${{allmap}} -v u=${{uniq}} '{{a+=1; b=$4+$5; i+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {{m+=$4; n+=b;}}; if (b>0) {{c+=1; d+=b;}}; if (b>2) e+=1}} END {{if (n>0) {{o=m/n*100;}} else o="NA"; print l,t,s,r,g,x,y" ("y/x*100"%)",z" ("z/x*100"%)",u" ("u/x*100"%)",c/a*100,e/a*100,i/a,d/c,o}}' >> "{output.stat_file}"

        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_pe__{params.sample_name}.html" --dir results/mC/reports/ --alignment_report {input.metrics_map} --dedup_report {input.metrics_dedup} --splitting_report results/mC/mapped/PE__{params.sample_name}.deduplicated_splitting_report.txt --mbias_report results/mC/mapped/PE__{params.sample_name}.deduplicated.M-bias.txt --nucleotide_report {params.prefix}/trim__{params.sample_name}__R1_bismark_bt2_pe.nucleotide_stats.txt
        cp results/mC/mapped/PE__"{params.sample_name}"*.txt results/mC/reports/
        cp {params.prefix}/trim__"{params.sample_name}"*.txt results/mC/reports/
        """
        
rule make_mc_stats_se:
    input:
        metrics_trim = "results/mC/reports/trim_se__{sample_name}.txt",
        metrics_map = "results/mC/mapped/{sample_name}/trim__{sample_name}__R0_bismark_bt2_SE_report.txt",
        metrics_dedup = "results/mC/mapped/{sample_name}/SE__{sample_name}.deduplication_report.txt",
        cx_report = "results/mC/methylcall/{sample_name}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        stat_file = "results/mC/reports/summary_mC_SE_mapping_stats_{sample_name}.txt",
        reportfile = "results/mC/reports/final_report_se__{sample_name}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        line = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        tissue = lambda wildcards: parse_sample_name(wildcards.sample_name)['tissue'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"results/mC/mapped/{wildcards.sample_name}",
        trimmed_fastq = config['trimmed_fastqs']
    log:
        temp(return_log_mc("{sample_name}", "making_stats", "SE"))
    conda: CONDA_ENV_MC
    threads: config["resources"]["make_mc_stats_se"]["threads"]
    resources:
        mem_mb=config["resources"]["make_mc_stats_se"]["mem_mb"],
        tmp_mb=config["resources"]["make_mc_stats_se"]["tmp_mb"],
        qos=config["resources"]["make_mc_stats_se"]["qos"]
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(grep "Total reads processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        else
            tot=$(grep "Sequences analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        fi
        filt=$(grep "Sequences analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Sequences did not map uniquely" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Number of alignments with a unique best hit" "{input.metrics_map}" | awk '{{print $NF}}')
        uniq=$(grep "Total count of deduplicated leftover sequences" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > {output.stat_file}
        ## Can change the name of the plastid chromosome to calculate non-conversion rate
        zcat {input.cx_report} | awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v s={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v x=${{tot}} -v y=${{filt}} -v z=${{allmap}} -v u=${{uniq}} '{{a+=1; b=$4+$5; i+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {{m+=$4; n+=b;}}; if (b>0) {{c+=1; d+=b;}}; if (b>2) e+=1}} END {{if (n>0) {{o=m/n*100;}} else o="NA"; print l,t,s,r,g,x,y" ("y/x*100"%)",z" ("z/x*100"%)",u" ("u/x*100"%)",c/a*100,e/a*100,i/a,d/c,o}}' >> "{output.stat_file}"

        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_se__{params.sample_name}.html" --dir results/mC/reports/ --alignment_report {input.metrics_map} --dedup_report {input.metrics_dedup} --splitting_report results/mC/mapped/SE__{params.sample_name}.deduplicated_splitting_report.txt --mbias_report results/mC/mapped/SE__{params.sample_name}.deduplicated.M-bias.txt --nucleotide_report {params.prefix}/trim__{params.sample_name}__R0_bismark_bt2.nucleotide_stats.txt
        mv results/mC/mapped/SE__"{params.sample_name}"*.txt results/mC/reports/
        mv {params.prefix}/trim__"{params.sample_name}"*.txt results/mC/reports/
        """

rule merging_mc_replicates:
    input:
        report_files = lambda wildcards: [ f"results/mC/methylcall/{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{replicate}__{wildcards.ref_genome}.deduplicated.CX_report.txt.gz" 
                                      for replicate in analysis_to_replicates.get((wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), []) ]
    output:
        bedfile = temp("results/mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.bed"),
        tempmergefile = temp("results/mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.merged.CX_report.txt"),
        mergefile = temp("results/mC/methylcall/{data_type}__{line}__{tissue}__{sample_type}__merged__{ref_genome}.merged.CX_report.txt.gz")
    params:
        sname = lambda wildcards: sample_name_str(wildcards, 'analysis')
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__{sample_type}__{ref_genome}", "merging_reps", ""))
    conda: CONDA_ENV_MC
    threads: config["resources"]["merging_mc_replicates"]["threads"]
    resources:
        mem_mb=config["resources"]["merging_mc_replicates"]["mem_mb"],
        tmp_mb=config["resources"]["merging_mc_replicates"]["tmp_mb"],
        qos=config["resources"]["merging_mc_replicates"]["qos"]
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
        cx_report = define_cx_report_input,
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        bigwigcg = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CG.bw",
        bigwigchg = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHG.bw",
        bigwigchh = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHH.bw",
        touch = "results/mC/chkpts/bigwig__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    wildcard_constraints:
        sample_type = r"(mC|WGBS|Pico|EMseq|merged)"  # Only Bismark sample types (ONT/bedMethyl handled by make_ont_bigwig_files)
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        ref_genome = lambda wildcards: wildcards.ref_genome,
        context = config['mC_context']
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "bigwig", ""))
    conda: CONDA_ENV_MC
    threads: config["resources"]["make_mc_bigwig_files"]["threads"]
    resources:
        mem_mb=config["resources"]["make_mc_bigwig_files"]["mem_mb"],
        tmp_mb=config["resources"]["make_mc_bigwig_files"]["tmp_mb"],
        qos=config["resources"]["make_mc_bigwig_files"]["qos"]
    shell:
        """
        {{
        if [[ "{params.context}" == "all" ]]; then
            zcat {input.cx_report} | awk -v OFS="\t" -v s={params.sample_name} '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "results/mC/tracks/"s"__CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "results/mC/tracks/"s"__CHG.bedGraph"; else print $1,$2-1,$2,$4/a*100 > "results/mC/tracks/"s"__CG.bedGraph"}}'
            for strand in plus minus; do
                case "${{strand}}" in 
                    plus)	sign="+";;
                    minus)	sign="-";;
                esac
                zcat {input.cx_report} | awk -v n=${{sign}} '$3==n' | awk -v OFS="\t" -v s={params.sample_name} -v d=${{strand}} '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "results/mC/tracks/"s"__CHH__"d".bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "results/mC/tracks/"s"__CHG__"d".bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 > "results/mC/tracks/"s"__CG__"d".bedGraph"}}'
            done
            for context in CG CHG CHH; do
                printf "\nMaking bigwig files of ${{context}} context for {params.sample_name}\n"
                LC_COLLATE=C sort -k1,1 -k2,2n results/mC/tracks/{params.sample_name}__${{context}}.bedGraph > results/mC/tracks/sorted__{params.sample_name}__${{context}}.bedGraph
                bedGraphToBigWig results/mC/tracks/sorted__{params.sample_name}__${{context}}.bedGraph {input.chrom_sizes} results/mC/tracks/{params.sample_name}__${{context}}.bw
                for strand in plus minus
                do
                    printf "\nMaking ${{strand}} strand bigwig files of ${{context}} context for {params.sample_name}\n"
                    LC_COLLATE=C sort -k1,1 -k2,2n results/mC/tracks/{params.sample_name}__${{context}}__${{strand}}.bedGraph > results/mC/tracks/sorted__{params.sample_name}__${{context}}__${{strand}}.bedGraph
                    bedGraphToBigWig results/mC/tracks/sorted__{params.sample_name}__${{context}}__${{strand}}.bedGraph {input.chrom_sizes} results/mC/tracks/{params.sample_name}__${{context}}__${{strand}}.bw
                done
            done
            rm -f results/mC/tracks/*"{params.sample_name}"*bedGraph*
        elif [[ "{params.context}" == "CG-only" ]]; then
            zcat {input.cx_report} | awk -v OFS="\t" '($4+$5)>0 {{a=$4+$5; print $1,$2-1,$2,$4/a*100}}' > "results/mC/tracks/"{params.sample_name}"__CG.bedGraph"
            for strand in plus minus; do
                case "${{strand}}" in 
                    plus)	sign="+";;
                    minus)	sign="-";;
                esac
                zcat {input.cx_report} | awk -v n=${{sign}} '$3==n' | awk -v OFS="\t" '($4+$5)>0 {{a=$4+$5; print $1,$2-1,$2,$4/a*100}}' > "results/mC/tracks/"{params.sample_name}"__CG__"${{strand}}".bedGraph"
            done
            printf "\nMaking bigwig files of CG context for {params.sample_name}\n"
            LC_COLLATE=C sort -k1,1 -k2,2n results/mC/tracks/{params.sample_name}__CG.bedGraph > results/mC/tracks/sorted__{params.sample_name}__CG.bedGraph
            bedGraphToBigWig results/mC/tracks/sorted__{params.sample_name}__CG.bedGraph {input.chrom_sizes} results/mC/tracks/{params.sample_name}__CG.bw
            for strand in plus minus
            do
                printf "\nMaking ${{strand}} strand bigwig files of CG context for {params.sample_name}\n"
                LC_COLLATE=C sort -k1,1 -k2,2n results/mC/tracks/{params.sample_name}__CG__${{strand}}.bedGraph > results/mC/tracks/sorted__{params.sample_name}__CG__${{strand}}.bedGraph
                bedGraphToBigWig results/mC/tracks/sorted__{params.sample_name}__CG__${{strand}}.bedGraph {input.chrom_sizes} results/mC/tracks/{params.sample_name}__CG__${{strand}}.bw
            done
            touch {output.bigwigchg} # they are required for downstream rules
            touch {output.bigwigchh} # they are required for downstream rules
            rm -f results/mC/tracks/*"{params.sample_name}"*bedGraph*
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
        dmr_summary = "results/mC/DMRs/summary__{sample1}__vs__{sample2}__DMRs.txt"
    params:
        script = script_DMRs(),
        context = config['mC_context'],
        sample1 = lambda wildcards: wildcards.sample1,
        sample2 = lambda wildcards: wildcards.sample2,
        nb_sample1 = lambda wildcards: len(define_DMR_samples(wildcards.sample1)),
        nb_sample2 = lambda wildcards: len(define_DMR_samples(wildcards.sample2))
    log:
        temp(return_log_mc("{sample1}__vs__{sample2}", "DMRs", ""))
    conda: CONDA_ENV_MC
    threads: config["resources"]["call_DMRs_pairwise"]["threads"]
    resources:
        mem_mb=config["resources"]["call_DMRs_pairwise"]["mem_mb"],
        tmp_mb=config["resources"]["call_DMRs_pairwise"]["tmp_mb"],
        qos=config["resources"]["call_DMRs_pairwise"]["qos"]
    shell:
        """
        {{
        printf "running DMRcaller for {params.sample1} vs {params.sample2}\n"
        Rscript "{params.script}" "{threads}" "{input.chrom_sizes}" "{params.context}" "{params.sample1}" "{params.sample2}" "{params.nb_sample1}" "{params.nb_sample2}" {input.sample1} {input.sample2}
        }} 2>&1 | tee -a "{log}"
        """    

rule all_mc:
    input:
        final = lambda wildcards: define_final_mC_output(wildcards.ref_genome)
    output:
        touch = "results/mC/chkpts/mC_analysis__{analysis_name}__{ref_genome}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """

################################################################################
# ONT Direct Methylation Rules
################################################################################

CONDA_ENV_ONT=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_ont.yaml")
MODKIT_VERSION = "0.6.1"
MODKIT_BIN = os.path.join(REPO_FOLDER, "workflow", "bin", "modkit")

rule download_modkit:
    """Download modkit binary from GitHub releases (not available via conda)."""
    output:
        binary = os.path.join(REPO_FOLDER, "workflow", "bin", "modkit")
    params:
        version = MODKIT_VERSION,
        bin_dir = os.path.join(REPO_FOLDER, "workflow", "bin")
    log:
        temp(os.path.join(REPO_FOLDER, "results", "logs", "download_modkit.log"))
    shell:
        """
        {{
        mkdir -p {params.bin_dir}

        # Download modkit release (u16 = Ubuntu 16 build, compatible with CentOS 7+)
        MODKIT_URL="https://github.com/nanoporetech/modkit/releases/download/v{params.version}/modkit_v{params.version}_u16_x86_64.tar.gz"
        printf "Downloading modkit v{params.version} from $MODKIT_URL\n"

        curl -fSL "$MODKIT_URL" -o /tmp/modkit.tar.gz
        tar -xzf /tmp/modkit.tar.gz -C {params.bin_dir}
        rm -f /tmp/modkit.tar.gz

        # Move modkit binary from extracted subdirectory to bin/
        # The tarball extracts to dist_modkit_v<version>_<hash>/modkit
        extracted_modkit=$(find {params.bin_dir} -name "modkit" -type f | head -1)
        if [[ -n "$extracted_modkit" && "$extracted_modkit" != "{output.binary}" ]]; then
            mv "$extracted_modkit" {output.binary}
            # Cleanup extracted directory
            rm -rf {params.bin_dir}/dist_modkit_*
        fi

        # Make executable
        chmod +x {output.binary}

        # Verify installation
        {output.binary} --version

        printf "modkit installed successfully at {output.binary}\n"
        }} > {log} 2>&1
        """

def get_ont_pileup_input(wildcards):
    """Determine input for modkit_pileup based on sample type."""
    sample_name = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}"
    sample_type = wildcards.sample_type
    if sample_type == "bedMethyl":
        # bedMethyl input - skip pileup, go directly to split
        return f"results/mC/ont/validated__{sample_name}.bedmethyl.gz"
    else:
        # modBAM input - need alignment check
        return f"results/mC/ont/aligned__{sample_name}.bam"

def define_ont_DMR_samples(sample_name):
    """Get bedMethyl files for ONT DMR analysis."""
    data_type = get_sample_info_from_name(sample_name, analysis_samples, 'data_type')
    line = get_sample_info_from_name(sample_name, analysis_samples, 'line')
    tissue = get_sample_info_from_name(sample_name, analysis_samples, 'tissue')
    sample_type = get_sample_info_from_name(sample_name, analysis_samples, 'sample_type')
    ref_genome = get_sample_info_from_name(sample_name, analysis_samples, 'ref_genome')
    replicates = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])

    return [ f"results/mC/ont/pileup__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bedmethyl.gz"
                    for replicate in replicates ]

rule make_modkit_context_beds:
    """Generate sorted, tabix-indexed CG/CHG/CHH context BED files for a reference genome.

    These are only needed for bedMethyl sample_type inputs. For ONT samples,
    modkit_pileup filters by context directly using --motif.
    """
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        modkit = MODKIT_BIN
    output:
        cg_bed = "genomes/{ref_genome}/modkit_CG.bed.gz",
        cg_tbi = "genomes/{ref_genome}/modkit_CG.bed.gz.tbi",
        chg_bed = "genomes/{ref_genome}/modkit_CHG.bed.gz",
        chg_tbi = "genomes/{ref_genome}/modkit_CHG.bed.gz.tbi",
        chh_bed = "genomes/{ref_genome}/modkit_CHH.bed.gz",
        chh_tbi = "genomes/{ref_genome}/modkit_CHH.bed.gz.tbi"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","modkit_context_beds_{ref_genome}.log"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["make_modkit_context_beds"]["threads"]
    resources:
        mem_mb=config["resources"]["make_modkit_context_beds"]["mem_mb"],
        tmp_mb=config["resources"]["make_modkit_context_beds"]["tmp_mb"],
        qos=config["resources"]["make_modkit_context_beds"]["qos"]
    shell:
        """
        {{
        printf "\nGenerating sorted, indexed methylation context BED files for {params.ref_genome}\n"

        # Generate CG context bed (sorted and bgzipped for tabix)
        printf "Generating CG context...\n"
        {input.modkit} motif bed {input.fasta} CG 0 | sort -k1,1 -k2,2n | bgzip -@ {threads} > {output.cg_bed}
        tabix -p bed {output.cg_bed}

        # Generate CHG context bed (CAG, CTG, CCG)
        printf "Generating CHG context...\n"
        {input.modkit} motif bed {input.fasta} CHG 0 | sort -k1,1 -k2,2n | bgzip -@ {threads} > {output.chg_bed}
        tabix -p bed {output.chg_bed}

        # Generate CHH context bed (CAA, CAT, CAC, CTA, CTT, CTC, CCA, CCT, CCC)
        printf "Generating CHH context...\n"
        {input.modkit} motif bed {input.fasta} CHH 0 | sort -k1,1 -k2,2n | bgzip -@ {threads} > {output.chh_bed}
        tabix -p bed {output.chh_bed}

        printf "\nContext BED files created and indexed successfully\n"
        }} 2>&1 | tee -a "{log}"
        """

rule get_modbam:
    """Acquire and validate a modBAM file with MM/ML methylation tags."""
    input:
        modbam = lambda wildcards: get_sample_info_from_name(
            f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
            samples, 'fastq_path'
        ),
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        validated_bam = "results/mC/ont/validated__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bam",
        validated_bai = "results/mC/ont/validated__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bam.bai"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        validate_script = os.path.join(REPO_FOLDER,"workflow","scripts","validate_ont_input.py")
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "get_modbam", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["get_modbam"]["threads"]
    resources:
        mem_mb=config["resources"]["get_modbam"]["mem_mb"],
        tmp_mb=config["resources"]["get_modbam"]["tmp_mb"],
        qos=config["resources"]["get_modbam"]["qos"]
    shell:
        """
        {{
        printf "\nValidating modBAM for {params.sample_name}\n"

        # Validate MM/ML tags exist
        python {params.validate_script} modBAM {input.modbam} {input.chrom_sizes}

        # Link or copy the validated BAM
        if [[ "{input.modbam}" == *.bam ]]; then
            ln -sf $(realpath {input.modbam}) {output.validated_bam}
        else
            # If not a .bam extension, copy it
            cp {input.modbam} {output.validated_bam}
        fi

        # Index the BAM
        samtools index -@ {threads} {output.validated_bam}

        printf "\nmodBAM validated and indexed\n"
        }} 2>&1 | tee -a "{log}"
        """

rule get_bedmethyl:
    """Acquire and validate a pre-computed bedMethyl file."""
    input:
        bedmethyl = lambda wildcards: get_sample_info_from_name(
            f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__bedMethyl__{wildcards.replicate}__{wildcards.ref_genome}",
            samples, 'fastq_path'
        ),
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        validated_bed = "results/mC/ont/validated__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}.bedmethyl.gz"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__bedMethyl__{wildcards.replicate}__{wildcards.ref_genome}",
        validate_script = os.path.join(REPO_FOLDER,"workflow","scripts","validate_ont_input.py")
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}", "get_bedmethyl", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["get_bedmethyl"]["threads"]
    resources:
        mem_mb=config["resources"]["get_bedmethyl"]["mem_mb"],
        tmp_mb=config["resources"]["get_bedmethyl"]["tmp_mb"],
        qos=config["resources"]["get_bedmethyl"]["qos"]
    shell:
        """
        {{
        printf "\nValidating bedMethyl for {params.sample_name}\n"

        # Validate bedMethyl format
        python {params.validate_script} bedMethyl {input.bedmethyl} {input.chrom_sizes}

        # Copy/compress the validated file
        if [[ "{input.bedmethyl}" == *.gz ]]; then
            cp {input.bedmethyl} {output.validated_bed}
        else
            pigz -p {threads} -c {input.bedmethyl} > {output.validated_bed}
        fi

        printf "\nbedMethyl validated\n"
        }} 2>&1 | tee -a "{log}"
        """

rule align_modbam:
    """Realign modBAM to reference genome if unaligned or reference mismatch."""
    input:
        modbam = "results/mC/ont/validated__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.bam",
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        aligned_bam = "results/mC/ont/aligned__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.bam",
        aligned_bai = "results/mC/ont/aligned__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.bam.bai"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__ONT__{wildcards.replicate}__{wildcards.ref_genome}",
        preset = config.get('ont_methylation', {}).get('alignment', {}).get('preset', 'lr:hqae')
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}", "align_modbam", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["align_modbam"]["threads"]
    resources:
        mem_mb=config["resources"]["align_modbam"]["mem_mb"],
        tmp_mb=config["resources"]["align_modbam"]["tmp_mb"],
        qos=config["resources"]["align_modbam"]["qos"]
    shell:
        """
        {{
        printf "\nChecking alignment status for {params.sample_name}\n"

        # Check if BAM is aligned by looking for @SQ headers
        has_sq=$(samtools view -H {input.modbam} | grep -c "^@SQ" || true)

        # Check if aligned to correct reference
        needs_realign=false

        if [[ "$has_sq" -eq 0 ]]; then
            printf "BAM is unaligned, will align to reference\n"
            needs_realign=true
        else
            # Check chromosome overlap with reference
            bam_chroms=$(samtools view -H {input.modbam} | grep "^@SQ" | cut -f2 | sed 's/SN://' | sort | head -20)
            ref_chroms=$(cut -f1 {input.chrom_sizes} | sort | head -20)
            n_ref_chroms=$(wc -l < {input.chrom_sizes})
            overlap=$(comm -12 <(echo "$bam_chroms") <(echo "$ref_chroms") | wc -l)
            # Require overlap of at least min(5, number of ref chromosomes)
            min_overlap=5
            if [[ "$n_ref_chroms" -lt "$min_overlap" ]]; then
                min_overlap=$n_ref_chroms
            fi
            if [[ "$overlap" -lt "$min_overlap" ]]; then
                printf "Low chromosome overlap ($overlap/$n_ref_chroms) with reference, will realign\n"
                needs_realign=true
            fi
        fi

        if [[ "$needs_realign" == "true" ]]; then
            printf "Aligning modBAM to {wildcards.ref_genome} with mm2plus\n"

            # Extract reads preserving MM/ML tags, align with mm2plus (accelerated minimap2), sort
            samtools fastq -T MM,ML {input.modbam} | \
                mm2plus -ax {params.preset} -t {threads} -y {input.fasta} - | \
                samtools sort -@ {threads} -o {output.aligned_bam} -

            samtools index -@ {threads} {output.aligned_bam}
        else
            printf "BAM is already aligned to compatible reference, linking\n"
            ln -sf $(realpath {input.modbam}) {output.aligned_bam}
            samtools index -@ {threads} {output.aligned_bam}
        fi

        printf "\nAlignment complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule modkit_summary:
    """Generate QC statistics from modBAM using modkit summary."""
    input:
        bam = "results/mC/ont/aligned__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.bam",
        modkit = MODKIT_BIN
    output:
        summary = "results/mC/ont/summary__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.txt"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__ONT__{wildcards.replicate}__{wildcards.ref_genome}"
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}", "modkit_summary", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["modkit_summary"]["threads"]
    resources:
        mem_mb=config["resources"]["modkit_summary"]["mem_mb"],
        tmp_mb=config["resources"]["modkit_summary"]["tmp_mb"],
        qos=config["resources"]["modkit_summary"]["qos"]
    shell:
        """
        {{
        printf "\nGenerating modkit summary for {params.sample_name}\n"

        {input.modkit} summary --threads {threads} {input.bam} > {output.summary}

        printf "\nSummary complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule modkit_summary_bedmethyl:
    """Generate placeholder summary for bedMethyl inputs."""
    input:
        bedmethyl = "results/mC/ont/validated__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}.bedmethyl.gz"
    output:
        summary = "results/mC/ont/summary__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}.txt"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__bedMethyl__{wildcards.replicate}__{wildcards.ref_genome}"
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}", "modkit_summary_bedmethyl", "ONT"))
    conda: CONDA_ENV_ONT
    threads: 1
    resources:
        mem_mb=config["resources"]["get_bedmethyl"]["mem_mb"],
        tmp_mb=config["resources"]["get_bedmethyl"]["tmp_mb"],
        qos=config["resources"]["get_bedmethyl"]["qos"]
    shell:
        """
        {{
        printf "\nGenerating summary for pre-computed bedMethyl {params.sample_name}\n"

        # Count lines and contexts
        total_sites=$(zcat {input.bedmethyl} | grep -v "^#" | wc -l)
        printf "Sample: {params.sample_name}\n" > {output.summary}
        printf "Input type: pre-computed bedMethyl\n" >> {output.summary}
        printf "Total methylation sites: $total_sites\n" >> {output.summary}

        # Count by context if name column contains context info
        printf "\nContext counts (from name column):\n" >> {output.summary}
        zcat {input.bedmethyl} | grep -v "^#" | cut -f4 | sort | uniq -c | sort -rn >> {output.summary} || true

        printf "\nSummary complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule modkit_pileup:
    """Generate context-specific bedMethyl files using modkit pileup with --motif filtering.

    This is more efficient than running a single pileup and splitting afterwards,
    as modkit filters during processing rather than requiring slow bedtools intersect.

    For CG context, generates both:
    - Strand-combined output (--cpg --combine-strands): one entry per CpG dinucleotide
    - Stranded output (--motif CG 0): separate entries for + and - strand cytosines
    """
    input:
        bam = "results/mC/ont/aligned__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.bam",
        bai = "results/mC/ont/aligned__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}.bam.bai",
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        modkit = MODKIT_BIN
    output:
        cg_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}__CG.bed.gz",
        cg_stranded_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}__CG_stranded.bed.gz",
        chg_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}__CHG.bed.gz",
        chh_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}__CHH.bed.gz"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__ONT__{wildcards.replicate}__{wildcards.ref_genome}",
        combine_mods = "--combine-mods" if config.get('ont_methylation', {}).get('pileup', {}).get('combine_mods', True) else ""
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__ONT__{replicate}__{ref_genome}", "modkit_pileup", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["modkit_pileup"]["threads"]
    resources:
        mem_mb=config["resources"]["modkit_pileup"]["mem_mb"],
        tmp_mb=config["resources"]["modkit_pileup"]["tmp_mb"],
        qos=config["resources"]["modkit_pileup"]["qos"]
    shell:
        """
        {{
        printf "\nRunning modkit pileup with context filtering for {params.sample_name}\n"

        # CG context pileup - strand-combined for primary output (one entry per CpG dinucleotide)
        # --modified-bases C with --combine-mods combines 5mC+5hmC into total cytosine methylation
        printf "Processing CG context (strand-combined)...\n"
        {input.modkit} pileup \
            --threads {threads} \
            --ref {input.fasta} \
            {params.combine_mods} \
            --modified-bases C \
            --filter-threshold 0.75 \
            --cpg --combine-strands \
            {input.bam} \
            /dev/stdout | pigz -p {threads} > {output.cg_bed}

        # CG context pileup - stranded for strand-specific bigwigs (separate + and - entries)
        printf "Processing CG context (stranded)...\n"
        {input.modkit} pileup \
            --threads {threads} \
            --ref {input.fasta} \
            {params.combine_mods} \
            --modified-bases C \
            --filter-threshold 0.75 \
            --motif CG 0 \
            {input.bam} \
            /dev/stdout | pigz -p {threads} > {output.cg_stranded_bed}

        # CHG context pileup (inherently stranded - not palindromic)
        printf "Processing CHG context...\n"
        {input.modkit} pileup \
            --threads {threads} \
            --ref {input.fasta} \
            {params.combine_mods} \
            --modified-bases C \
            --filter-threshold 0.75 \
            --motif CHG 0 \
            {input.bam} \
            /dev/stdout | pigz -p {threads} > {output.chg_bed}

        # CHH context pileup (inherently stranded - not palindromic)
        printf "Processing CHH context...\n"
        {input.modkit} pileup \
            --threads {threads} \
            --ref {input.fasta} \
            {params.combine_mods} \
            --modified-bases C \
            --filter-threshold 0.75 \
            --motif CHH 0 \
            {input.bam} \
            /dev/stdout | pigz -p {threads} > {output.chh_bed}

        printf "\nPileup complete for all contexts\n"
        }} 2>&1 | tee -a "{log}"
        """

rule copy_bedmethyl_for_pileup:
    """Copy pre-computed bedMethyl to pileup location for consistent downstream processing."""
    input:
        bedmethyl = "results/mC/ont/validated__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}.bedmethyl.gz"
    output:
        bedmethyl = "results/mC/ont/pileup__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}.bedmethyl.gz"
    localrule: True
    shell:
        """
        cp {input.bedmethyl} {output.bedmethyl}
        """

rule split_bedmethyl_by_context:
    """Split pre-computed bedMethyl file by methylation context using tabix for fast lookups.

    This rule is only used for bedMethyl sample_type inputs where we don't have the original
    modBAM to re-run modkit pileup with --motif filtering. For ONT samples, modkit_pileup
    handles context splitting directly.
    """
    input:
        bedmethyl = "results/mC/ont/pileup__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}.bedmethyl.gz",
        cg_bed = "genomes/{ref_genome}/modkit_CG.bed.gz",
        cg_tbi = "genomes/{ref_genome}/modkit_CG.bed.gz.tbi",
        chg_bed = "genomes/{ref_genome}/modkit_CHG.bed.gz",
        chg_tbi = "genomes/{ref_genome}/modkit_CHG.bed.gz.tbi",
        chh_bed = "genomes/{ref_genome}/modkit_CHH.bed.gz",
        chh_tbi = "genomes/{ref_genome}/modkit_CHH.bed.gz.tbi"
    output:
        cg_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}__CG.bed.gz",
        chg_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}__CHG.bed.gz",
        chh_bed = "results/mC/ont/context__{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}__CHH.bed.gz"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__bedMethyl__{wildcards.replicate}__{wildcards.ref_genome}"
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__bedMethyl__{replicate}__{ref_genome}", "split_bedmethyl", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["split_bedmethyl_by_context"]["threads"]
    resources:
        mem_mb=config["resources"]["split_bedmethyl_by_context"]["mem_mb"],
        tmp_mb=config["resources"]["split_bedmethyl_by_context"]["tmp_mb"],
        qos=config["resources"]["split_bedmethyl_by_context"]["qos"]
    shell:
        """
        {{
        printf "\nSplitting bedMethyl by context for {params.sample_name}\n"

        # Sort and index bedmethyl for tabix lookups
        printf "Preparing bedMethyl for context lookup...\n"
        sorted_bed="results/mC/ont/tmp_sorted_{params.sample_name}.bed.gz"
        zcat {input.bedmethyl} | sort -k1,1 -k2,2n | bgzip -@ {threads} > "$sorted_bed"
        tabix -p bed "$sorted_bed"

        # Use tabix intersect for fast lookups (context BEDs are pre-indexed)
        printf "Extracting CG context...\n"
        tabix "$sorted_bed" -R {input.cg_bed} 2>/dev/null | pigz -p {threads} > {output.cg_bed}

        printf "Extracting CHG context...\n"
        tabix "$sorted_bed" -R {input.chg_bed} 2>/dev/null | pigz -p {threads} > {output.chg_bed}

        printf "Extracting CHH context...\n"
        tabix "$sorted_bed" -R {input.chh_bed} 2>/dev/null | pigz -p {threads} > {output.chh_bed}

        # Cleanup temp files
        rm -f "$sorted_bed" "$sorted_bed.tbi"

        printf "\nContext splitting complete\n"
        }} 2>&1 | tee -a "{log}"
        """

def get_ont_bigwig_inputs(wildcards):
    """Get inputs for ONT bigwig generation, including stranded CG bed for ONT samples."""
    base_inputs = {
        "cg_bed": f"results/mC/ont/context__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}__CG.bed.gz",
        "chg_bed": f"results/mC/ont/context__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}__CHG.bed.gz",
        "chh_bed": f"results/mC/ont/context__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}__CHH.bed.gz",
        "chrom_sizes": f"genomes/{wildcards.ref_genome}/chrom.sizes"
    }
    # ONT samples have a dedicated stranded CG bed from modkit_pileup
    if wildcards.sample_type == "ONT":
        base_inputs["cg_stranded_bed"] = f"results/mC/ont/context__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}__CG_stranded.bed.gz"
    return base_inputs

rule make_ont_bigwig_files:
    """Generate bigwig files from ONT bedMethyl data (column 11 = percent modified).

    For CG context, generates both combined and strand-specific bigwigs (matching Bismark output).
    - Combined CG bigwig: from strand-combined bed (one value per CpG dinucleotide)
    - Strand-specific bigwigs: from stranded bed (ONT) or filtered by strand column (bedMethyl)
    """
    input:
        unpack(get_ont_bigwig_inputs)
    output:
        bigwig_cg = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CG.bw",
        bigwig_cg_plus = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CG__plus.bw",
        bigwig_cg_minus = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CG__minus.bw",
        bigwig_chg = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHG.bw",
        bigwig_chg_plus = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHG__plus.bw",
        bigwig_chg_minus = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHG__minus.bw",
        bigwig_chh = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHH.bw",
        bigwig_chh_plus = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHH__plus.bw",
        bigwig_chh_minus = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHH__minus.bw",
        touch = "results/mC/chkpts/bigwig__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
    wildcard_constraints:
        sample_type = r"(ONT|bedMethyl)"  # Only match ONT and bedMethyl samples
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        sample_type = lambda wildcards: wildcards.sample_type,
        context = config['mC_context'],
        # Pass stranded CG bed path as param since it's only available for ONT samples
        cg_stranded_bed = lambda wildcards: f"results/mC/ont/context__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}__CG_stranded.bed.gz" if wildcards.sample_type == "ONT" else ""
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "ont_bigwig", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["make_ont_bigwig_files"]["threads"]
    resources:
        mem_mb=config["resources"]["make_ont_bigwig_files"]["mem_mb"],
        tmp_mb=config["resources"]["make_ont_bigwig_files"]["tmp_mb"],
        qos=config["resources"]["make_ont_bigwig_files"]["qos"]
    shell:
        """
        {{
        printf "\nMaking bigwig files for {params.sample_name}\n"

        # Helper function to create bigwig from bedGraph
        make_bigwig() {{
            local input_bed="$1"
            local output_bw="$2"
            local strand_filter="$3"  # "+" or "-" or "" for all

            if [[ -n "$strand_filter" ]]; then
                # Filter by strand (column 6) then extract columns 1,2,3,11
                zcat "$input_bed" | awk -v OFS="\t" -v strand="$strand_filter" \
                    'NF>=11 && $6==strand {{print $1,$2,$3,$11}}' | \
                    LC_COLLATE=C sort -k1,1 -k2,2n > "${{output_bw}}.bedGraph"
            else
                # No strand filter
                zcat "$input_bed" | awk -v OFS="\t" 'NF>=11 {{print $1,$2,$3,$11}}' | \
                    LC_COLLATE=C sort -k1,1 -k2,2n > "${{output_bw}}.bedGraph"
            fi

            if [[ -s "${{output_bw}}.bedGraph" ]]; then
                bedGraphToBigWig "${{output_bw}}.bedGraph" {input.chrom_sizes} "$output_bw"
            else
                # Create empty bigwig placeholder if no data
                touch "$output_bw"
            fi
            rm -f "${{output_bw}}.bedGraph"
        }}

        if [[ "{params.context}" == "all" || "{params.context}" == "CG-only" ]]; then
            # CG combined bigwig (strand-combined, one value per CpG)
            printf "Creating CG combined bigwig...\n"
            make_bigwig "{input.cg_bed}" "{output.bigwig_cg}" ""

            # CG strand-specific bigwigs
            printf "Creating CG strand-specific bigwigs...\n"
            if [[ "{params.sample_type}" == "ONT" ]]; then
                # ONT: use dedicated stranded bed file (path passed via params)
                make_bigwig "{params.cg_stranded_bed}" "{output.bigwig_cg_plus}" "+"
                make_bigwig "{params.cg_stranded_bed}" "{output.bigwig_cg_minus}" "-"
            else
                # bedMethyl: filter combined bed by strand column
                make_bigwig "{input.cg_bed}" "{output.bigwig_cg_plus}" "+"
                make_bigwig "{input.cg_bed}" "{output.bigwig_cg_minus}" "-"
            fi
        fi

        if [[ "{params.context}" == "all" ]]; then
            # CHG bigwigs (combined and strand-specific)
            printf "Creating CHG bigwigs...\n"
            make_bigwig "{input.chg_bed}" "{output.bigwig_chg}" ""
            make_bigwig "{input.chg_bed}" "{output.bigwig_chg_plus}" "+"
            make_bigwig "{input.chg_bed}" "{output.bigwig_chg_minus}" "-"

            # CHH bigwigs (combined and strand-specific)
            printf "Creating CHH bigwigs...\n"
            make_bigwig "{input.chh_bed}" "{output.bigwig_chh}" ""
            make_bigwig "{input.chh_bed}" "{output.bigwig_chh_plus}" "+"
            make_bigwig "{input.chh_bed}" "{output.bigwig_chh_minus}" "-"
        elif [[ "{params.context}" == "CG-only" ]]; then
            # Create empty placeholder files for CHG/CHH
            touch {output.bigwig_chg} {output.bigwig_chg_plus} {output.bigwig_chg_minus}
            touch {output.bigwig_chh} {output.bigwig_chh_plus} {output.bigwig_chh_minus}
        else
            printf "Unknown sequence context selection! Check config 'mC_context'\n"
            exit 1
        fi

        touch {output.touch}
        printf "\nBigwig creation complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule call_DMRs_modkit:
    """Call DMRs between two ONT samples using modkit dmr pair."""
    input:
        sample1_beds = lambda wildcards: define_ont_DMR_samples(wildcards.sample1),
        sample2_beds = lambda wildcards: define_ont_DMR_samples(wildcards.sample2),
        fasta = lambda wildcards: f"genomes/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}.fa",
        modkit = MODKIT_BIN
    output:
        dmr_summary = "results/mC/DMRs/summary__{sample1}__vs__{sample2}__DMRs.txt"
    params:
        sample1 = lambda wildcards: wildcards.sample1,
        sample2 = lambda wildcards: wildcards.sample2,
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')
    log:
        temp(return_log_mc("{sample1}__vs__{sample2}", "modkit_DMRs", "ONT"))
    conda: CONDA_ENV_ONT
    threads: config["resources"]["call_DMRs_modkit"]["threads"]
    resources:
        mem_mb=config["resources"]["call_DMRs_modkit"]["mem_mb"],
        tmp_mb=config["resources"]["call_DMRs_modkit"]["tmp_mb"],
        qos=config["resources"]["call_DMRs_modkit"]["qos"]
    shell:
        """
        {{
        printf "\nRunning modkit DMR analysis: {params.sample1} vs {params.sample2}\n"

        # Merge sample bedmethyl files if multiple replicates
        sample1_merged="results/mC/DMRs/tmp__{params.sample1}.merged.bed"
        sample2_merged="results/mC/DMRs/tmp__{params.sample2}.merged.bed"

        zcat {input.sample1_beds} | sort -k1,1 -k2,2n > "$sample1_merged"
        zcat {input.sample2_beds} | sort -k1,1 -k2,2n > "$sample2_merged"

        # Run modkit dmr pair with segmentation
        {input.modkit} dmr pair \
            --threads {threads} \
            --ref {input.fasta} \
            --segment \
            -a "$sample1_merged" {params.sample1} \
            -b "$sample2_merged" {params.sample2} \
            -o results/mC/DMRs/dmr__{params.sample1}__vs__{params.sample2}

        # Create summary file
        printf "DMR Analysis Summary\n" > {output.dmr_summary}
        printf "Sample 1: {params.sample1}\n" >> {output.dmr_summary}
        printf "Sample 2: {params.sample2}\n" >> {output.dmr_summary}
        printf "Reference: {params.ref_genome}\n" >> {output.dmr_summary}
        printf "Method: modkit dmr pair --segment\n" >> {output.dmr_summary}
        printf "\n" >> {output.dmr_summary}

        # Add DMR counts
        if [[ -f "results/mC/DMRs/dmr__{params.sample1}__vs__{params.sample2}.bed" ]]; then
            dmr_count=$(wc -l < "results/mC/DMRs/dmr__{params.sample1}__vs__{params.sample2}.bed")
            printf "Total DMRs identified: $dmr_count\n" >> {output.dmr_summary}
            cat "results/mC/DMRs/dmr__{params.sample1}__vs__{params.sample2}.bed" >> {output.dmr_summary}
        else
            printf "No DMR output file found\n" >> {output.dmr_summary}
        fi

        # Cleanup temp files
        rm -f "$sample1_merged" "$sample2_merged"

        printf "\nDMR analysis complete\n"
        }} 2>&1 | tee -a "{log}"
        """        
