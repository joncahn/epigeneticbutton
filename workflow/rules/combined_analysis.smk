# function to access logs more easily
def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER,"results","combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_samplenames_per_env_and_ref(wildcards):
    names = []
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    if globenv == "all_chip":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:    
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        env2 = row.env
        file = f"results/{env2}/peaks/selected_peaks__{spname}.bedPeak"
        if env2 == "TF":
            label=f"{row.line}_{row.tissue}_{row.extra_info}"
        else:
            label=f"{row.line}_{row.tissue}_{row.sample_type}"
        
        names.append(f"{label}:{file}")
    
    return names

def define_peakfiles_per_env_and_ref(wildcards):
    files = []
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    if globenv == "all_chip":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:    
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        env2 = row.env
        files.append(f"results/{env2}/peaks/selected_peaks__{spname}.bedPeak")
    
    return files

def define_samplelabels_per_env_and_ref(wildcards):
    labels = []
    ref_genome = wildcards.ref_genome
    env = wildcards.env
    if env == "all_chip":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:    
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == env) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        if env == "TF":
            label=f"{row.line}_{row.tissue}_{row.extra_info}"
        else:
            label=f"{row.line}_{row.tissue}_{row.sample_type}"
        labels.append(label)
    
    return labels

def define_input_bedfile(wildcards):
    bedname = wildcards.bedname
    ref_genome = wildcards.ref_genome
    
    if bedname.startswith("combined_peakfiles"):
        file = f"results/combined/bedfiles/{bedname}__{ref_genome}.bed"
        
    return file

def define_sample_types_for_upset(wildcards):
    types = []
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    if globenv == "all_chip":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:    
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        if row.env == "ChIP":
            types.append(row.sample_type)
        elif row.env == "TF":
            types.append(row.extra_info)
        
    return types

def define_final_combined_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    analysis_name = config['analysis_name']
    text_files = []
    plot_files = []
    
    chip_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    if len(chip_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__ChIP__{analysis_name}__{ref_genome}.pdf")
    
    tf_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'TF') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    if len(tf_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__TF__{analysis_name}__{ref_genome}.pdf")
    
    if len(chip_analysis_samples) >=1 and len(tf_analysis_samples) >=1:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__all_chip__{analysis_name}__{ref_genome}.pdf")
    
    if analysis:
        results = plot_files + text_files
    else:
        results = []
    
    return results

###
# Rules to prep and then plot the mapping stats:
rule prepping_mapping_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"results/{wildcards.env}/reports/summary_{wildcards.env}_{get_sample_info_from_name(sample_name, samples, 'paired')}_mapping_stats_{sample_name}.txt"
                                                for sample_name in get_sample_names_by_env(wildcards.env, samples) ]
    output:
        temp_stat_file = temp("results/combined/reports/temp_summary_mapping_stats_{analysis_name}_{env}.txt"),
        stat_file = "results/combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "prep_mapping_stats"))
    conda: CONDA_ENV
    threads: config["resources"]["prepping_mapping_stats"]["threads"]
    resources:
        mem=config["resources"]["prepping_mapping_stats"]["mem"],
        tmp=config["resources"]["prepping_mapping_stats"]["tmp"]
    shell:
        """
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk -F "\t" -v OFS="\t" 'NR>1 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_mapping_stats:
    input:
        summary_stats = "results/combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    output:
        plot = "results/combined/plots/mapping_stats_{analysis_name}_{env}.pdf"
    params:
        analysis_name = lambda wildcards: wildcards.analysis_name,
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_mapping_stats.R")
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_mapping_stats"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_mapping_stats"]["threads"]
    resources:
        mem=config["resources"]["plotting_mapping_stats"]["mem"],
        tmp=config["resources"]["plotting_mapping_stats"]["tmp"]
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}" "{output.plot}"
        """
        
###
# Rules to prep and then plot the peak stats:
rule prepping_chip_peak_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"results/{wildcards.env}/reports/summary_{wildcards.env}_peak_stats_{sample_name}.txt" for sample_name in get_sample_names_by_env(wildcards.env, analysis_samples) ]
    output:
        temp_stat_file = temp("results/combined/reports/temp_summary_peak_stats_{analysis_name}_{env}.txt"),
        stat_file = "results/combined/reports/summary_peak_stats_{analysis_name}_{env}.txt"
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "prep_peak_stats"))
    conda: CONDA_ENV
    threads: config["resources"]["prepping_chip_peak_stats"]["threads"]
    resources:
        mem=config["resources"]["prepping_chip_peak_stats"]["mem"],
        tmp=config["resources"]["prepping_chip_peak_stats"]["tmp"]
    shell:
        """
        printf "Line\tTissue\tSample\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tPeaks_in_merged\tPeaks_in_pseudo_reps\tPeaks_in_idr\tSelected_peaks\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk 'NR>1' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_peaks_stats_chip_tf:
    input:
        summary_stats = "results/combined/reports/summary_peak_stats_{analysis_name}_{env}.txt"
    output:
        plot = "results/combined/plots/peak_stats_{analysis_name}_{env}.pdf"
    params:
        analysis_name = lambda wildcards: wildcards.analysis_name,
        env = lambda wildcards: wildcards.env,
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_peak_stats.R")
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_peak_stats"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_peaks_stats_chip_tf"]["threads"]
    resources:
        mem=config["resources"]["plotting_peaks_stats_chip_tf"]["mem"],
        tmp=config["resources"]["plotting_peaks_stats_chip_tf"]["tmp"]
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}" "{output.plot}" "{params.env}"
        """

###
# Rules to prep and then plot the sizes stats for sRNA
rule prepping_srna_sizes_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"results/sRNA/reports/sizes_stats__{sample_name}.txt"
                                                for sample_name in get_sample_names_by_env(wildcards.env, samples) ]
    output:
        temp_stat_file = temp("results/combined/reports/temp_summary_sizes_stats_{analysis_name}_{env}.txt"),
        stat_file = "results/combined/reports/summary_sizes_stats_{analysis_name}_{env}.txt"
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "prep_srna_sizes"))
    conda: CONDA_ENV
    threads: config["resources"]["prepping_srna_sizes_stats"]["threads"]
    resources:
        mem=config["resources"]["prepping_srna_sizes_stats"]["mem"],
        tmp=config["resources"]["prepping_srna_sizes_stats"]["tmp"]
    shell:
        """
        printf "Sample\tType\tSize\tCount\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk -F "\t" -v OFS="\t" 'NR>1' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_srna_sizes_stats:
    input:
        summary_stats = "results/combined/reports/summary_sizes_stats_{analysis_name}_{env}.txt"
    output:
        plot1 = "results/combined/plots/srna_sizes_stats_{analysis_name}_{env}.pdf",
        plot2 = "results/combined/plots/srna_sizes_stats_zoom_{analysis_name}_{env}.pdf"
    params:
        analysis_name = lambda wildcards: wildcards.analysis_name,
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_sizes_stats.R")
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_srna_sizes"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_srna_sizes_stats"]["threads"]
    resources:
        mem=config["resources"]["plotting_srna_sizes_stats"]["mem"],
        tmp=config["resources"]["plotting_srna_sizes_stats"]["tmp"]
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}"
        """


###
# Rules to prep and plot ChIP upset plots
rule combine_peakfiles:
    input:
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        peakfiles = lambda wildcards: define_peakfiles_per_env_and_ref(wildcards)
    output:
        temp1_file = temp("results/combined/bedfiles/temp1_combined_peakfiles__{env}__{analysis_name}__{ref_genome}.bed"),
        temp2_file = temp("results/combined/bedfiles/temp2_combined_peakfiles__{env}__{analysis_name}__{ref_genome}.bed"),
        merged_file = "results/combined/bedfiles/combined_peakfiles__{env}__{analysis_name}__{ref_genome}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samplenames_per_env_and_ref(wildcards),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_peaks_{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["combine_peakfiles"]["threads"]
    resources:
        mem=config["resources"]["combine_peakfiles"]["mem"],
        tmp=config["resources"]["combine_peakfiles"]["tmp"]
    shell:
        """
        {{
        printf "Merging peakfiles for {params.env} {params.analysis_name} {params.ref_genome}\n"
        for pair in {params.names}; do
            label=$(echo ${{pair}} | cut -d":" -f1)
            file=$(echo ${{pair}} | cut -d":" -f2)
            awk -v OFS="\t" -v l=${{label}} '{{print $1,$2,$3,l}}' ${{file}} >> {output.temp1_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        printf "Chr\tStart\tStop\tPeakID\tSamples\n" > {output.merged_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v e={params.env} -v a={params.analysis_name} '{{print $1,$2,$3,"combined_peak_"e"_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule get_annotations_for_bedfile:
    input:
        bedfile = lambda wildcards: define_input_bedfile(wildcards),
        region_file = lambda wildcards: f"results/combined/tracks/{wildcards.ref_genome}__all_genes.bed",
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes"
    output:
        temp_bedfile = temp("results/combined/bedfiles/temp__{bedname}__{ref_genome}.bed"),
        annotated_file = "results/combined/bedfiles/annotated__{bedname}__{ref_genome}.bed"
    params:
        bedname = lambda wildcards: wildcards.bedname
    log:
        temp(return_log_combined("{bedname}", "{ref_genome}", "annotate_bedfile"))
    conda: CONDA_ENV
    threads: config["resources"]["get_annotations_for_bedfile"]["threads"]
    resources:
        mem=config["resources"]["get_annotations_for_bedfile"]["mem"],
        tmp=config["resources"]["get_annotations_for_bedfile"]["tmp"]
    shell:
        """
        {{
        printf "Annotating {params.bedname} to the closest genes\n"
        # checking for presence of header
        read -r chrom start end _ < {input.bedfile}
        if [[ "${{start}}" =~ ^[0-9]+$ ]] && [[ "${{end}}" =~ ^[0-9]+$ ]]; then
            awk -v OFS="\t" -v n={params.bedname} '{{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile}
        else
            awk -v OFS="\t" -v n={params.bedname} 'NR>1 {{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile}
        fi
        printf "Chr\tStart\tStop\tPeakID\tDistance\tGene_strand\tGID\tCategory\n" > {output.annotated_file}
        bedtools closest -a {output.temp_bedfile} -b {input.region_file} -g {input.chrom_sizes} -D ref | awk -v OFS="\t" '{{if ($10=="+") print $1,$2,$3,$4,$11,$10,$8; else print $1,$2,$3,$4,-$11,$10,$8}}' | awk -F"[=;]" -v OFS="\t" '{{print $1,$2}}' | sed 's/gene://' | awk -v OFS="\t" '{{if ($5<-2000) {{d="Distal_downstream"}} else if ($5<0) {{d="Terminator"}} else if ($5==0) {{d="Gene_body"}} else if ($5>2000) {{d="Distal_upstream"}} else {{d="Promoter"}}; print $1,$2,$3,$4,$5,$6,$8,d}}' >> {output.annotated_file}
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_upset_peaks:
    input:
        mergedfile = "results/combined/bedfiles/combined_peakfiles__{env}__{analysis_name}__{ref_genome}.bed",
        annotatedfile = "results/combined/bedfiles/annotated__combined_peakfiles__{env}__{analysis_name}__{ref_genome}.bed"
    output:
        plot = "results/combined/plots/Upset_combined_peaks__{env}__{analysis_name}__{ref_genome}.pdf"
    params:
        env = lambda wildcards: wildcards.env,
        types = lambda wildcards: define_sample_types_for_upset(wildcards),
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot.R")
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_upset_{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_upset_peaks"]["threads"]
    resources:
        mem=config["resources"]["plotting_upset_peaks"]["mem"],
        tmp=config["resources"]["plotting_upset_peaks"]["tmp"]
    shell:
        """
        Rscript "{params.script}" "{input.mergedfile}" "{input.annotatedfile}" "{params.env}" "{params.types}" "{output.plot}"
        """

###
# # rule to plot heatmaps
# rule making_deeptools_matrix_on_targetfile:
    # input:
        # bigwigs = lambda wildcards: define_bigwigs_per_env_and_ref(wildcards),
        # target_file = lambda wildcards: define_rnaseq_target_file(wildcards)
    # output:
        # plot = "results/RNA/plots/Heatmap__{analysis_name}__{ref_genome}__{target_name}.pdf"
    # params:
        # analysis_name = config['analysis_name'],
        # ref_genome = lambda wildcards: wildcards.ref_genome,
        # target_name = lambda wildcards: wildcards.target_name,
        # labels = lambda wildcards: define_labels_per_env_and_ref(wildcards)
    # log:
        # temp(return_log_rna("{ref_genome}", "plot_expression_{target_name}", "{analysis_name}"))
    # conda: CONDA_ENV
    # threads: config["resources"]["making_deeptools_matrix_on_targetfile"]["threads"]
    # resources:
        # mem=config["resources"]["making_deeptools_matrix_on_targetfile"]["mem"],
        # tmp=config["resources"]["making_deeptools_matrix_on_targetfile"]["tmp"]
    # shell:
        # """
        # printf "running plot expression levels for {input.target_file} (from {params.analysis_name} and {params.ref_genome})\n"
        # """
        
# rule plotting_heatmap_on_targetfile:
    # input:
        # bigwigs = lambda wildcards: define_bigwigs_per_genome(wildcards),
        # target_file = lambda wildcards: define_rnaseq_target_file(wildcards)
    # output:
        # plot = "results/RNA/plots/Heatmap__{analysis_name}__{ref_genome}__{target_name}.pdf"
    # params:
        # analysis_name = config['analysis_name'],
        # ref_genome = lambda wildcards: wildcards.ref_genome,
        # target_name = lambda wildcards: wildcards.target_name,
        # labels = lambda wildcards: define_env_samplelabels_per_ref(wildcards)
    # log:
        # temp(return_log_rna("{ref_genome}", "plot_expression_{target_name}", "{analysis_name}"))
    # conda: CONDA_ENV
    # threads: config["resources"]["plotting_heatmap_on_targetfile"]["threads"]
    # resources:
        # mem=config["resources"]["plotting_heatmap_on_targetfile"]["mem"],
        # tmp=config["resources"]["plotting_heatmap_on_targetfile"]["tmp"]
    # shell:
        # """
        # printf "running plot expression levels for {input.target_file} (from {params.analysis_name} and {params.ref_genome})\n"        
        # """

###
# final rule
rule all_combined:
    input:
        expand("results/combined/plots/mapping_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","TF","RNA","mC"]]),
        expand("results/combined/plots/srna_sizes_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["sRNA"]]),
        expand("results/combined/plots/peak_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","TF"]]),
        final = lambda wildcards: define_final_combined_output(wildcards.ref_genome)
    output:
        touch = "results/combined/chkpts/combined_analysis__{analysis_name}__{ref_genome}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """        
