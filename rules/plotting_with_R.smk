### R plots

# Rules to prep and then plot the mapping stats
rule prepping_mapping_stats:
    input:
        [
            f"{env}/chkpts/map__{sample_name}.done"
            for env in UNIQUE_ENVS
            for sample_name in get_sample_names_by_env(env, samples)
        ],
        "{env}/reports/summary_{env}_mapping_stats_{analysis_name}.txt"
    output:
        stat_file = "combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    params:
        sample_file = config['sample_file'],
        name = lambda wildcards: f"{wildcards.analysis_name}_{wildcards.env}"
    shell:
        """
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > "{output.stat_file}"
        while read data_type line tissue sample_type replicate ref_genome seq_id fastq_path paired env sample_name
        do
            awk -v a=${{line}} -v b=${{tissue}} -v c=${{sample_type}} -v d=${{replicate}} -v e=${{ref_genome}} '$1==a && $2==b && $3==c && $4==d && $5==e' {input.input_file} >> "combined/reports/temp_summary_mapping_stats_{params.name}.txt"
        done < "{params.sample_file}"
        sort combined/reports/temp_summary_mapping_stats_{params.name}.txt -u >> "{output.stat_file}"
        rm -f "combined/reports/temp_summary_mapping_stats_{params.name}.txt"
        """
    
rule plotting_mapping_stats_chip_rna:
    input:
        summary_stats = "combined/reports/summary_mapping_stats_{analysis_name}_{env}.txt"
    output:
        plot = "combined/plots/mapping_stats_{analysis_name}_{env}.pdf"
    params:
        analysisname = lambda wildcards: f"{wildcards.analysis_name}"
    conda: "envs/r_plotting.yaml"
    script:
        "scripts/MaizeCode_R_mapping_stats.R"