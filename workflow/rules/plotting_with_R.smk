###
# Rules to prep and then plot the mapping stats
rule prepping_mapping_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"results/{wildcards.env}/reports/summary_{wildcards.env}_{get_sample_info_from_name(sample_name, samples, 'paired')}_mapping_stats_{sample_name}.txt"
                                                for sample_name in get_sample_names_by_env(wildcards.env, samples) ]
    output:
        temp_stat_file = temp("results/combined/reports/temp_summary_mapping_stats_{analysis_name}_{env}.txt"),
        stat_file = "results/combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
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
        "results/{env}/logs/plotting_mapping_stats_{analysis_name}_{env}.log"
    conda: CONDA_ENV
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
        "results/{env}/logs/plotting_peaks_stats_{analysis_name}_{env}.log"
    conda: CONDA_ENV
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
        "results/sRNA/logs/plotting_srna_sizes_stats_{analysis_name}_{env}.log"
    conda: CONDA_ENV
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}"
        """