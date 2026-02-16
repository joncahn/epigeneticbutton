CONDA_ENV_MC=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_mc.yaml")

def return_log_mc(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"results","mC","logs",f"tmp__{sample_name}__{step}__{paired}.log")
     
def parameters_for_mc(sample_name):
    temp = parse_sample_name(sample_name)['sample_type']
    options = {"WGBS", "Pico", "EMseq", "dmC"}
    return temp if temp in options else "default"

def is_dmc_sample(sample_name):
    """Check if a sample uses direct methylation (dmC) workflow (vs bisulfite)."""
    return parse_sample_name(sample_name)['sample_type'] == "dmC"

def define_cx_report_input(wildcards):
    """Get CX_report path for a sample (used by bigwig generation and replicate merging).

    Routes to appropriate path based on sample_type:
    - Bismark samples: results/mC/methylcall/...deduplicated.CX_report.txt.gz
    - dmC (direct methylation): results/mC/dmc/cx_report__...CX_report.txt.gz
    - Merged replicates: results/mC/methylcall/...merged.CX_report.txt.gz
    """
    name = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}"
    if wildcards.replicate == "merged":
        return f"results/mC/methylcall/{name}.merged.CX_report.txt.gz"
    elif wildcards.sample_type == "dmC":
        # dmC samples: use converted CX_report from bedMethyl
        return f"results/mC/dmc/cx_report__{name}.CX_report.txt.gz"
    else:
        # Bismark samples: use deduplicated CX_report
        return f"results/mC/methylcall/{name}.deduplicated.CX_report.txt.gz"

def define_DMR_samples(sample_name):
    """Get CX_report files for DMR analysis.

    For Bismark samples: returns deduplicated CX_report files
    For dmC samples: returns converted CX_report files (from bedMethyl)
    """
    data_type = get_sample_info_from_name(sample_name, analysis_samples, 'data_type')
    line = get_sample_info_from_name(sample_name, analysis_samples, 'line')
    tissue = get_sample_info_from_name(sample_name, analysis_samples, 'tissue')
    sample_type = get_sample_info_from_name(sample_name, analysis_samples, 'sample_type')
    ref_genome = get_sample_info_from_name(sample_name, analysis_samples, 'ref_genome')

    # Return empty list if sample not found (prevents None in paths)
    if any(x is None for x in [data_type, line, tissue, sample_type, ref_genome]):
        return []

    replicates = analysis_to_replicates.get((data_type, line, tissue, sample_type, ref_genome), [])

    if sample_type == "dmC":
        # dmC samples: use converted CX_report files (unified format with all contexts)
        return [ f"results/mC/dmc/cx_report__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.CX_report.txt.gz"
                        for replicate in replicates ]
    else:
        # Bismark samples: use deduplicated CX_report files
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
    mC_context = config['mC_context']
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

        # dmC samples use direct methylation workflow
        if sample_type == "dmC":
            bigwig_files.append(f"results/mC/chkpts/bigwig__{sname}.done")
            ont_files.append(f"results/mC/dmc/summary__{sname}.txt")  # modkit summary
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
    
    # DMR analysis: all sample types use DMRcaller via unified CX_report format
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
        sample_name = r"(?!.*__dmC__).*"  # Exclude dmC (direct methylation) samples
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
        sample_name = r"(?!.*__dmC__).*"  # Exclude dmC (direct methylation) samples
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
        sample_name = r"(?!.*__dmC__).*"  # Exclude dmC (direct methylation) samples
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
    wildcard_constraints:
        sample_name = r"(?!.*__dmC__).*"  # Exclude dmC (direct methylation) samples
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
    wildcard_constraints:
        sample_name = r"(?!.*__dmC__).*"  # Exclude dmC (direct methylation) samples
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

def get_cx_reports_for_merging(wildcards):
    """Get CX_report paths for all replicates of a sample, for merging.

    Routes to appropriate paths based on sample_type:
    - Bismark samples: results/mC/methylcall/...deduplicated.CX_report.txt.gz
    - dmC (direct methylation): results/mC/dmc/cx_report__...CX_report.txt.gz
    """
    replicates = analysis_to_replicates.get(
        (wildcards.data_type, wildcards.line, wildcards.tissue, wildcards.sample_type, wildcards.ref_genome), [])

    if wildcards.sample_type == "dmC":
        return [f"results/mC/dmc/cx_report__{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{rep}__{wildcards.ref_genome}.CX_report.txt.gz"
                for rep in replicates]
    else:
        return [f"results/mC/methylcall/{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{rep}__{wildcards.ref_genome}.deduplicated.CX_report.txt.gz"
                for rep in replicates]

rule merging_mc_replicates:
    input:
        report_files = get_cx_reports_for_merging
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
    """Generate bigwig files from CX_report data.

    Handles all sample types (Bismark and dmC) via unified CX_report format.
    """
    input:
        cx_report = define_cx_report_input,
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        bigwigcg = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CG.bw",
        bigwigchg = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHG.bw",
        bigwigchh = "results/mC/tracks/{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}__CHH.bw",
        touch = "results/mC/chkpts/bigwig__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.done"
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
    """Call DMRs between two samples using DMRcaller.

    Works with both Bismark and dmC (direct methylation) samples:
    - Bismark: uses CX_report files from bismark_methylation_extractor
    - dmC: uses CX_report files converted from bedMethyl format
    """
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
# Direct Methylation (dmC) Rules
# Handles both modBAM (direct methylation basecalls) and pre-computed bedMethyl inputs
################################################################################

CONDA_ENV_DMC=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_dmc.yaml")
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

rule get_dmc_input:
    """Acquire and validate a direct methylation input file (modBAM or bedMethyl).

    Automatically detects the input type based on file extension and content:
    - modBAM: BAM files with MM/ML methylation tags
    - bedMethyl: pre-computed methylation calls in BED format

    Supports both:
    - Direct file paths: /path/to/sample.bam
    - Directory paths: /path/to/dir (finds file matching seq_id)

    When searching directories, bedMethyl files are preferred over modBAM if both
    exist with the same seq_id prefix (bedMethyl is pre-computed and ready to use).
    A warning is logged when a modBAM is skipped in favor of bedMethyl.

    Creates a marker file indicating the detected type for downstream rules.
    """
    input:
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt",
        validated = "results/mC/dmc/validated__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.input"
    wildcard_constraints:
        sample_type = r"dmC"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
        dmc_path = lambda wildcards: get_sample_info_from_name(
            f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
            samples, 'fastq_path'
        ),
        seq_id = lambda wildcards: get_sample_info_from_name(
            f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
            samples, 'seq_id'
        ),
        validate_script = os.path.join(REPO_FOLDER,"workflow","scripts","validate_dmc_input.py")
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "get_dmc_input", "dmC"))
    conda: CONDA_ENV_DMC
    threads: config["resources"]["get_modbam"]["threads"]
    resources:
        mem_mb=config["resources"]["get_modbam"]["mem_mb"],
        tmp_mb=config["resources"]["get_modbam"]["tmp_mb"],
        qos=config["resources"]["get_modbam"]["qos"]
    shell:
        """
        {{
        printf "\nDetecting and validating dmC input for {params.sample_name}\n"

        # Resolve input path: can be a file or a directory
        dmc_path="{params.dmc_path}"
        seq_id="{params.seq_id}"

        if [[ -f "$dmc_path" ]]; then
            # Direct file path provided
            input_file="$dmc_path"
            printf "Using direct file path: $input_file\n"
        elif [[ -d "$dmc_path" ]]; then
            # Directory path provided - find file matching seq_id
            # Uses *seq_id* pattern consistent with sample_download.smk
            # Prefers bedMethyl over modBAM if both exist (bedMethyl is pre-computed)
            printf "Searching directory for seq_id '$seq_id'...\n"

            # Check for bedMethyl files first (preferred)
            # Note: Use {{ || true; }} pattern to handle pipefail when ls finds no matches
            bedmethyl_file=""
            for ext in .bed.gz .bedmethyl .bed; do
                match_count=$( {{ ls -1 "$dmc_path"/*"$seq_id"*"$ext" 2>/dev/null || true; }} | wc -l)
                if [[ "$match_count" -eq 1 ]]; then
                    bedmethyl_file=$(ls "$dmc_path"/*"$seq_id"*"$ext")
                    break
                elif [[ "$match_count" -gt 1 ]]; then
                    printf "Error: Multiple bedMethyl files found matching seq_id '$seq_id' with extension '$ext':\n"
                    ls -1 "$dmc_path"/*"$seq_id"*"$ext"
                    printf "Please use a more specific seq_id or provide a direct file path.\n"
                    exit 1
                fi
            done

            # Check for modBAM files
            # Note: Use {{ || true; }} pattern to handle pipefail when ls finds no matches
            modbam_file=""
            modbam_count=$( {{ ls -1 "$dmc_path"/*"$seq_id"*.bam 2>/dev/null || true; }} | wc -l)
            if [[ "$modbam_count" -eq 1 ]]; then
                modbam_file=$(ls "$dmc_path"/*"$seq_id"*.bam)
            elif [[ "$modbam_count" -gt 1 ]]; then
                printf "Error: Multiple modBAM files found matching seq_id '$seq_id':\n"
                ls -1 "$dmc_path"/*"$seq_id"*.bam
                printf "Please use a more specific seq_id or provide a direct file path.\n"
                exit 1
            fi

            # Select input file: prefer bedMethyl over modBAM
            if [[ -n "$bedmethyl_file" ]]; then
                input_file="$bedmethyl_file"
                printf "Found bedMethyl: $input_file\n"
                if [[ -n "$modbam_file" ]]; then
                    printf "WARNING: modBAM also found ($modbam_file) but using bedMethyl (pre-computed calls preferred)\n"
                fi
            elif [[ -n "$modbam_file" ]]; then
                input_file="$modbam_file"
                printf "Found modBAM: $input_file\n"
            else
                printf "Error: No dmC file found matching seq_id '$seq_id' in '$dmc_path'\n"
                printf "Looked for patterns: *$seq_id*.bed.gz, *$seq_id*.bedmethyl, *$seq_id*.bed, *$seq_id*.bam\n"
                printf "Available files:\n"
                ls -la "$dmc_path"
                exit 1
            fi
        else
            printf "Error: Path does not exist: $dmc_path\n"
            exit 1
        fi

        # Auto-detect input type
        input_type=$(python {params.validate_script} detect "$input_file")
        printf "Detected input type: $input_type\n"

        # Write type marker for downstream rules
        echo "$input_type" > {output.type_marker}

        # Validate based on detected type
        python {params.validate_script} "$input_type" "$input_file" {input.chrom_sizes}

        # Create a symlink to the validated input
        ln -sf $(realpath "$input_file") {output.validated}

        printf "\nInput validated successfully\n"
        }} 2>&1 | tee -a "{log}"
        """

def get_dmc_input_type(wildcards):
    """Get the input type (modBAM or bedMethyl) for a dmC sample by reading the marker file."""
    sample_name = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}"
    marker_file = f"results/mC/dmc/input_type__{sample_name}.txt"
    # This function is called during DAG building, marker file may not exist yet
    # Return a checkpoint-compatible path
    return marker_file

checkpoint dmc_input_checkpoint:
    """Checkpoint to determine dmC input type for branching workflow."""
    input:
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt"
    output:
        touch = touch("results/mC/dmc/checkpoint__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.done")
    localrule: True

def get_pileup_input_for_dmc(wildcards):
    """Determine pileup input based on detected input type."""
    checkpoint_output = checkpoints.dmc_input_checkpoint.get(**wildcards).output[0]
    sample_name = f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}"
    type_marker = f"results/mC/dmc/input_type__{sample_name}.txt"
    with open(type_marker) as f:
        input_type = f.read().strip()
    if input_type == "modBAM":
        return f"results/mC/dmc/pileup_modbam__{sample_name}.bedmethyl.gz"
    else:
        return f"results/mC/dmc/pileup_bedmethyl__{sample_name}.bedmethyl.gz"

rule prepare_modbam_for_pileup:
    """Prepare modBAM input: index and optionally realign."""
    input:
        validated = "results/mC/dmc/validated__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.input",
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt",
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        aligned_bam = "results/mC/dmc/aligned__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bam",
        aligned_bai = "results/mC/dmc/aligned__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bam.bai"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
        preset = config.get('dmc_methylation', {}).get('alignment', {}).get('preset', 'lr:hqae')
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "prepare_modbam", "dmC"))
    conda: CONDA_ENV_DMC
    threads: config["resources"]["align_modbam"]["threads"]
    resources:
        mem_mb=config["resources"]["align_modbam"]["mem_mb"],
        tmp_mb=config["resources"]["align_modbam"]["tmp_mb"],
        qos=config["resources"]["align_modbam"]["qos"]
    shell:
        """
        {{
        # Check if input is modBAM
        input_type=$(cat {input.type_marker})
        if [[ "$input_type" != "modBAM" ]]; then
            printf "Skipping - input is not modBAM (is $input_type)\n"
            touch {output.aligned_bam} {output.aligned_bai}
            exit 0
        fi

        printf "\nChecking alignment status for {params.sample_name}\n"

        # Check if BAM is aligned by looking for @SQ headers
        has_sq=$(samtools view -H {input.validated} | grep -c "^@SQ" || true)

        # Check if aligned to correct reference
        needs_realign=false

        if [[ "$has_sq" -eq 0 ]]; then
            printf "BAM is unaligned, will align to reference\n"
            needs_realign=true
        else
            # Check chromosome overlap with reference
            bam_chroms=$(samtools view -H {input.validated} | grep "^@SQ" | cut -f2 | sed 's/SN://' | sort | head -20)
            ref_chroms=$(cut -f1 {input.chrom_sizes} | sort | head -20)
            n_ref_chroms=$(wc -l < {input.chrom_sizes})
            overlap=$(comm -12 <(echo "$bam_chroms") <(echo "$ref_chroms") | wc -l)
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
            samtools fastq -T MM,ML {input.validated} | \
                mm2plus -ax {params.preset} -t {threads} -y {input.fasta} - | \
                samtools sort -@ {threads} -o {output.aligned_bam} -
            samtools index -@ {threads} {output.aligned_bam}
        else
            printf "BAM is already aligned to compatible reference, linking\n"
            ln -sf $(realpath {input.validated}) {output.aligned_bam}
            samtools index -@ {threads} {output.aligned_bam}
        fi

        printf "\nAlignment complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule modkit_pileup_dmc:
    """Generate bedMethyl file from aligned modBAM using modkit pileup."""
    input:
        bam = "results/mC/dmc/aligned__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bam",
        bai = "results/mC/dmc/aligned__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bam.bai",
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt",
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        modkit = MODKIT_BIN
    output:
        bedmethyl = "results/mC/dmc/pileup_modbam__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bedmethyl.gz"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
        combine_mods = "--combine-mods" if config.get('dmc_methylation', {}).get('pileup', {}).get('combine_mods', True) else ""
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "modkit_pileup", "dmC"))
    conda: CONDA_ENV_DMC
    threads: config["resources"]["modkit_pileup"]["threads"]
    resources:
        mem_mb=config["resources"]["modkit_pileup"]["mem_mb"],
        tmp_mb=config["resources"]["modkit_pileup"]["tmp_mb"],
        qos=config["resources"]["modkit_pileup"]["qos"]
    shell:
        """
        {{
        printf "\nRunning modkit pileup for {params.sample_name}\n"

        {input.modkit} pileup \
            --threads {threads} \
            --ref {input.fasta} \
            {params.combine_mods} \
            --modified-bases C \
            --filter-threshold 0.75 \
            {input.bam} \
            /dev/stdout | pigz -p {threads} > {output.bedmethyl}

        printf "\nPileup complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule copy_bedmethyl_input:
    """Copy pre-computed bedMethyl to pileup location for consistent downstream processing."""
    input:
        validated = "results/mC/dmc/validated__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.input",
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt"
    output:
        bedmethyl = "results/mC/dmc/pileup_bedmethyl__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bedmethyl.gz"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}"
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "copy_bedmethyl", "dmC"))
    conda: CONDA_ENV_DMC
    threads: config["resources"]["get_bedmethyl"]["threads"]
    resources:
        mem_mb=config["resources"]["get_bedmethyl"]["mem_mb"],
        tmp_mb=config["resources"]["get_bedmethyl"]["tmp_mb"],
        qos=config["resources"]["get_bedmethyl"]["qos"]
    shell:
        """
        {{
        printf "\nCopying bedMethyl input for {params.sample_name}\n"

        # Copy/compress the validated file
        if [[ "{input.validated}" == *.gz ]]; then
            cp {input.validated} {output.bedmethyl}
        else
            pigz -p {threads} -c {input.validated} > {output.bedmethyl}
        fi

        printf "\nbedMethyl copied\n"
        }} 2>&1 | tee -a "{log}"
        """

rule merge_pileup_sources:
    """Create unified pileup path regardless of input type."""
    input:
        pileup = get_pileup_input_for_dmc
    output:
        bedmethyl = "results/mC/dmc/pileup__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bedmethyl.gz"
    localrule: True
    shell:
        """
        ln -sf $(basename {input.pileup}) {output.bedmethyl}
        """

rule modkit_summary_dmc:
    """Generate QC statistics from modBAM using modkit summary."""
    input:
        bam = "results/mC/dmc/aligned__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bam",
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt",
        modkit = MODKIT_BIN
    output:
        summary = "results/mC/dmc/summary__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}"
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "modkit_summary", "dmC"))
    conda: CONDA_ENV_DMC
    threads: config["resources"]["modkit_summary"]["threads"]
    resources:
        mem_mb=config["resources"]["modkit_summary"]["mem_mb"],
        tmp_mb=config["resources"]["modkit_summary"]["tmp_mb"],
        qos=config["resources"]["modkit_summary"]["qos"]
    shell:
        """
        {{
        input_type=$(cat {input.type_marker})

        if [[ "$input_type" == "modBAM" ]]; then
            printf "\nGenerating modkit summary for {params.sample_name}\n"
            {input.modkit} summary --threads {threads} {input.bam} > {output.summary}
        else
            printf "\nGenerating summary for pre-computed bedMethyl {params.sample_name}\n"
            printf "Sample: {params.sample_name}\n" > {output.summary}
            printf "Input type: pre-computed bedMethyl\n" >> {output.summary}
            printf "Note: Limited statistics available for bedMethyl input\n" >> {output.summary}
        fi

        printf "\nSummary complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule make_mc_stats_dmc:
    """Generate mapping stats for dmC samples in Bismark-compatible format.

    Produces stats file compatible with prepping_mapping_stats rule for combined reports.
    Uses CX_report file (unified format) for coverage statistics.
    """
    input:
        cx_report = "results/mC/dmc/cx_report__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.CX_report.txt.gz",
        type_marker = "results/mC/dmc/input_type__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    output:
        stat_file = "results/mC/reports/summary_mC_SE_mapping_stats_{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.txt"
    params:
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
        line = lambda wildcards: wildcards.line,
        tissue = lambda wildcards: wildcards.tissue,
        sample_type = "dmC",
        replicate = lambda wildcards: wildcards.replicate,
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "making_stats", "dmC"))
    conda: CONDA_ENV_DMC
    threads: config["resources"]["modkit_summary"]["threads"]
    resources:
        mem_mb=config["resources"]["modkit_summary"]["mem_mb"],
        tmp_mb=config["resources"]["modkit_summary"]["tmp_mb"],
        qos=config["resources"]["modkit_summary"]["qos"]
    shell:
        """
        {{
        printf "\nMaking mapping statistics summary for dmC sample {params.sample_name}\n"

        input_type=$(cat {input.type_marker})

        # Get genome size and coverage stats from CX_report (CG sites only for consistency)
        genome_size=$(awk '{{sum+=$2}} END {{print sum}}' {input.chrom_sizes})
        cov_stats=$(zcat {input.cx_report} | awk -v gs="$genome_size" '
            BEGIN {{total_cov=0; n_sites=0; n_sites_3x=0}}
            $6 == "CG" {{
                cov = $4 + $5;
                total_cov += cov;
                n_sites += 1;
                if (cov >= 3) n_sites_3x += 1;
            }}
            END {{
                if (n_sites > 0) {{
                    pct_cov = n_sites / gs * 100;
                    pct_cov_3x = n_sites_3x / gs * 100;
                    avg_cov_all = total_cov / gs;
                    avg_cov_covered = total_cov / n_sites;
                }} else {{
                    pct_cov = 0; pct_cov_3x = 0; avg_cov_all = 0; avg_cov_covered = 0;
                }}
                printf "%.4f\\t%.4f\\t%.4f\\t%.4f", pct_cov, pct_cov_3x, avg_cov_all, avg_cov_covered;
            }}
        ')

        pct_cov=$(echo "$cov_stats" | cut -f1)
        pct_cov_3x=$(echo "$cov_stats" | cut -f2)
        avg_cov_all=$(echo "$cov_stats" | cut -f3)
        avg_cov_covered=$(echo "$cov_stats" | cut -f4)

        # Write header
        printf "Line\\tTissue\\tSample\\tRep\\tReference_genome\\tTotal_reads\\tPassing_filtering\\tAll_mapped_reads\\tUniquely_mapped_reads\\tPercentage_covered\\tPercentage_covered_min3reads\\tAverage_coverage_all\\tAverage_coverage_covered\\tNon_conversion_rate(Pt/Lambda)\\n" > {output.stat_file}

        # For modBAM input, we can get read counts from the aligned BAM
        if [[ "$input_type" == "modBAM" ]]; then
            bam_file="results/mC/dmc/aligned__{params.sample_name}.bam"
            if [[ -f "$bam_file" ]]; then
                flagstat=$(samtools flagstat "$bam_file")
                tot=$(echo "$flagstat" | grep "in total" | awk '{{print $1}}')
                mapped=$(echo "$flagstat" | grep "mapped (" | head -1 | awk '{{print $1}}')
                uniq=$(samtools view -c -q 20 "$bam_file")
                if [ "$tot" -gt 0 ]; then
                    pct_mapped=$(awk "BEGIN {{printf \\"%.2f\\", $mapped/$tot*100}}")
                    pct_uniq=$(awk "BEGIN {{printf \\"%.2f\\", $uniq/$tot*100}}")
                else
                    pct_mapped="0.00"
                    pct_uniq="0.00"
                fi
                printf "{params.line}\\t{params.tissue}\\t{params.sample_type}\\t{params.replicate}\\t{params.ref_genome}\\t$tot\\t$tot (${{pct_mapped}}%%)\\t$mapped (${{pct_mapped}}%%)\\t$uniq (${{pct_uniq}}%%)\\t$pct_cov\\t$pct_cov_3x\\t$avg_cov_all\\t$avg_cov_covered\\tNA\\n" >> {output.stat_file}
            else
                printf "{params.line}\\t{params.tissue}\\t{params.sample_type}\\t{params.replicate}\\t{params.ref_genome}\\tNA\\tNA\\tNA\\tNA\\t$pct_cov\\t$pct_cov_3x\\t$avg_cov_all\\t$avg_cov_covered\\tNA\\n" >> {output.stat_file}
            fi
        else
            # For bedMethyl input, no BAM stats available
            printf "{params.line}\\t{params.tissue}\\t{params.sample_type}\\t{params.replicate}\\t{params.ref_genome}\\tNA\\tNA\\tNA\\tNA\\t$pct_cov\\t$pct_cov_3x\\t$avg_cov_all\\t$avg_cov_covered\\tNA\\n" >> {output.stat_file}
        fi

        printf "\ndmC stats complete\\n"
        }} 2>&1 | tee -a "{log}"
        """

rule convert_bedmethyl_to_cx_report:
    """Convert bedMethyl format to Bismark CX_report format.

    Determines context (CG/CHG/CHH) from the reference genome and outputs
    a single CX_report file matching the Bismark output format. This enables
    unified downstream processing with merging_mc_replicates, make_mc_bigwig_files,
    and call_DMRs_pairwise.

    When mC_context is 'CG-only', filters output to only include CG context.
    """
    input:
        bedmethyl = "results/mC/dmc/pileup__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.bedmethyl.gz",
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        fai = "genomes/{ref_genome}/{ref_genome}.fa.fai"
    output:
        cx_report = "results/mC/dmc/cx_report__{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}.CX_report.txt.gz"
    params:
        script = os.path.join(REPO_FOLDER, "workflow", "scripts", "bedmethyl_to_cx_report.py"),
        sample_name = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__dmC__{wildcards.replicate}__{wildcards.ref_genome}",
        context = config['mC_context']
    log:
        temp(return_log_mc("{data_type}__{line}__{tissue}__dmC__{replicate}__{ref_genome}", "bedmethyl_to_cx", "dmC"))
    conda: CONDA_ENV_DMC
    threads: 1
    resources:
        mem_mb=config["resources"]["convert_bedmethyl_to_cx_report"]["mem_mb"],
        tmp_mb=config["resources"]["convert_bedmethyl_to_cx_report"]["tmp_mb"],
        qos=config["resources"]["convert_bedmethyl_to_cx_report"]["qos"]
    shell:
        """
        {{
        printf "Converting bedMethyl to CX_report format for {params.sample_name}...\n"
        printf "Context mode: {params.context}\n"

        # Convert bedMethyl to CX_report (context determined from reference)
        python {params.script} {input.bedmethyl} {input.fasta} /dev/stdout > results/mC/dmc/tmp__{params.sample_name}.cx

        # Filter by context if CG-only mode
        if [[ "{params.context}" == "CG-only" ]]; then
            printf "Filtering to CG context only...\n"
            awk -F'\\t' '$6 == "CG"' results/mC/dmc/tmp__{params.sample_name}.cx > results/mC/dmc/tmp__{params.sample_name}_filtered.cx
            mv results/mC/dmc/tmp__{params.sample_name}_filtered.cx results/mC/dmc/tmp__{params.sample_name}.cx
        fi

        # Sort by chromosome and position, then compress
        printf "Sorting and compressing...\n"
        sort -k1,1 -k2,2n results/mC/dmc/tmp__{params.sample_name}.cx | pigz -p {threads} > {output.cx_report}
        rm -f results/mC/dmc/tmp__{params.sample_name}.cx

        printf "Conversion complete\n"
        }} 2>&1 | tee -a "{log}"
        """
