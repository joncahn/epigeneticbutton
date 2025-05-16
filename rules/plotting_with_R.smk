### R plots
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/r_plotting.yaml")

###
# Rules to prep and then plot the mapping stats
rule prepping_mapping_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"{wildcards.env}/reports/summary_{wildcards.env}_{get_sample_info_from_name(sample_name, samples, 'paired')}_mapping_stats_{sample_name}.txt"
                                                for sample_name in get_sample_names_by_env(wildcards.env, samples) ]
    output:
        stat_file = "combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt",
        temp_stat_file = temp("combined/reports/temp_summary_mapping_stats_{analysis_name}_{env}.txt")
    params:
        sample_file = config['sample_file'],
        analysis_name = lambda wildcards: wildcards.analysis_name,
        env = lambda wildcards: wildcards.env
    shell:
        """
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk 'NR>1' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"\
        """
    
rule plotting_mapping_stats_chip_rna:
    input:
        summary_stats = "combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    output:
        plot = "combined/plots/mapping_stats_{analysis_name}_{env}.pdf"
    params:
        analysisname = lambda wildcards: f"{wildcards.analysis_name}",
        script=os.path.join(REPO_FOLDER,"scripts/R_mapping_stats.R")
    log:
        "{env}/logs/plotting_mapping_stats_{analysis_name}_{env}.log"
    conda: CONDA_ENV
    script:
        {params.script}
        
###
# Rules to prep and then plot the peak stats
rule prepping_chip_peak_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"{wildcards.env}/reports/summary_{wildcards.env}_peak_stats_{sample_name}.txt" for sample_name in get_sample_names_by_env(wildcards.env, analysis_samples) ]
    output:
        stat_file = "combined/reports/summary_peak_stats_{analysis_name}_{env}.txt",
        tmp_stat_file = temp("combined/reports/temp_summary_peak_stats_{analysis_name}_{env}.txt")
    params:
        sample_file = lambda wildcards: f"{wildcards.analysis_name}__analysis_samplefile.txt",
        analysis_name = lambda wildcards: wildcards.analysis_name,
        env = lambda wildcards: wildcards.env
    shell:
        """
        printf "Line\tTissue\tSample\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk 'NR>1' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_peaks_stats_chip_tf:
    input:
        summary_stats = "combined/reports/summary_peak_stats_{analysis_name}_{env}.txt"
    output:
        plot = "combined/plots/peak_stats_{analysis_name}_{env}.pdf"
    params:
        analysisname = lambda wildcards: f"{wildcards.analysis_name}",
        env = lambda wildcards: f"{wildcards.env}",
        script=os.path.join(REPO_FOLDER,"scripts/R_peak_stats.R")
    log:
        "{env}/logs/plotting_peaks_stats_{analysis_name}_{env}.log"
    conda: CONDA_ENV
    script:
        {params.script}