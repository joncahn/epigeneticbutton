### R plots
CONDA_ENV=os.path.join(REPO_FOLDER,"envs/r_plotting.yaml")

###
# Rules to prep and then plot the mapping stats
rule prepping_mapping_stats:
    input:
        checkpoints = [ f"{env}/chkpts/map__{sample_name}.done"
                        for env in ["ChIP","RNA"]
                        if env in UNIQUE_ENVS
                        for sample_name in get_sample_names_by_env(env, samples) ],
        input_file = "{env}/reports/summary_{env}_mapping_stats_{analysis_name}.txt"
    output:
        stat_file = "combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    params:
        sample_file = config['sample_file'],
        analysis_name = lambda wildcards: wildcards.analysis_name,
        env = lambda wildcards: wildcards.env
    shell:
        """
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > "{output.stat_file}"
        while read data_type line tissue sample_type replicate ref_genome seq_id fastq_path paired env sample_name
        do
            if [[ "${{env}}" == "{params.env}" ]]; then
                awk -v a=${{line}} -v b=${{tissue}} -v c=${{sample_type}} -v d=${{replicate}} -v e=${{ref_genome}} '$1==a && $2==b && $3==c && $4==d && $5==e' {input.input_file} >> "combined/reports/temp_summary_mapping_stats_{params.analysis_name}_{params.env}.txt"
            fi
        done < "{params.sample_file}"
        sort combined/reports/temp_summary_mapping_stats_{params.analysis_name}_{params.env}.txt -u >> "{output.stat_file}"
        rm -f "combined/reports/temp_summary_mapping_stats_{params.analysis_name}_{params.env}.txt"
        """
    
rule plotting_mapping_stats_chip_rna:
    input:
        summary_stats = "combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    output:
        plot = "combined/plots/mapping_stats_{analysis_name}_{env}.pdf"
    params:
        analysisname = lambda wildcards: f"{wildcards.analysis_name}"
    conda: CONDA_ENV
    script:
        "scripts/R_mapping_stats.R"
        
###
# Rules to prep and then plot the peak stats
rule prepping_chip_peak_stats:
    input:
        checkpoints = [ f"{env}/logs/called_peaks__{sample_name}.log"
                        for env in ["ChIP", "TF"]
                        if env in UNIQUE_ENVS
                        for sample_name in get_sample_names_by_env("ChIP", analysis_samples) ],
        input_file = "{env}/reports/summary_{env}_peak_stats_{analysis_name}.txt"
    output:
        stat_file = "combined/reports/summary_peak_stats_{analysis_name}_{env}.txt"
    params:
        sample_file = lambda wildcards: f"{wildcards.analysis_name}__analysis_samplefile.txt",
        analysis_name = lambda wildcards: wildcards.analysis_name,
        env = lambda wildcards: wildcards.env
    shell:
        """
        printf "Line\tTissue\tSample\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tPeaks_in_merged\tPeaks_in_pseudo_reps\tSelected_peaks\n" > "{output.stat_file}"
        while read data_type line tissue sample_type replicate ref_genome seq_id fastq_path paired env sample_name
        do
            if [[ "${{env}}" == "{params.env}" ]]; then
                awk -v a=${{line}} -v b=${{tissue}} -v c=${{sample_type}} -v e=${{ref_genome}} '$1==a && $2==b && $3==c && $4==d' {input.input_file} >> "combined/reports/temp_summary_mapping_stats_{params.name}.txt"
            fi
        done < "{params.sample_file}"
        sort combined/reports/temp_summary_mapping_stats_{params.name}.txt -u >> "{output.stat_file}"
        rm -f "combined/reports/temp_summary_mapping_stats_{params.name}.txt"
        """
    
rule plotting_peaks_stats_chip_tf:
    input:
        summary_stats = "combined/reports/summary_peak_stats_{analysis_name}_{env}.txt"
    output:
        plot = "combined/plots/peak_stats_{analysis_name}_{env}.pdf"
    params:
        analysisname = lambda wildcards: f"{wildcards.analysis_name}",
        env = lambda wildcards: f"{wildcards.env}"
    conda: CONDA_ENV
    script:
        "scripts/R_peak_stats.R"