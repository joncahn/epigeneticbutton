CONDA_ENV_ATAC=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_chip.yaml")

def define_final_atac_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    trimmed_fastqs = config['trimmed_fastqs']
    aligned_bams = config['aligned_bams']
    map_files = []
    stat_files = []
    qc_files = []
    peak_files = []
    bigwig_files = []

    filtered_rep_samples = samples[ (samples['env'] == 'ATAC') & (samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        env = "ATAC"
        if paired == "PE" and not aligned_bams:
            qc_files.append(f"results/{env}/reports/trim__{sname}__R1_fastqc.html")
            qc_files.append(f"results/{env}/reports/trim__{sname}__R2_fastqc.html")
            map_files.append(f"results/{env}/logs/process_chip_pe_sample__{sname}.log")
            if not trimmed_fastqs:
                qc_files.append(f"results/{env}/reports/raw__{sname}__R1_fastqc.html")
                qc_files.append(f"results/{env}/reports/raw__{sname}__R2_fastqc.html")
        elif paired == "SE" and not aligned_bams:
            qc_files.append(f"results/{env}/reports/trim__{sname}__R0_fastqc.html")
            map_files.append(f"results/{env}/logs/process_chip_se_sample__{sname}.log")
            if not trimmed_fastqs:
                qc_files.append(f"results/{env}/reports/raw__{sname}__R0_fastqc.html")

    # ATAC has no Input samples to filter out
    peaktype = config["atac_callpeaks"]["peaktype"]
    for _, row in filtered_rep_samples.iterrows():
        sname = sample_name_str(row, 'sample')
        paired = get_sample_info_from_name(sname, samples, 'paired')
        bigwig_files.append(f"results/ATAC/tracks/coverage__final__{sname}.bw")
        peak_files.append(f"results/ATAC/peaks/peaks_atac__final__{sname}_peaks.{peaktype}Peak")

    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ATAC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        peak_files.append(f"results/ATAC/peaks/selected_peaks__{spname}.bedPeak")
        reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
        if len(reps) >= 2:
            bigwig_files.append(f"results/ATAC/tracks/coverage__merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.bw")
            stat_files.append(f"results/ATAC/chkpts/idr__{spname}.done")

    results = map_files + bigwig_files

    if qc_option == "all":
        results += qc_files

    if analysis:
        results += peak_files + stat_files

    return results


rule atac_shift_and_bed:
    input:
        bamfile = "results/ATAC/mapped/{file_type}__{sample_name}.bam"
    output:
        bedfile = "results/ATAC/mapped/shifted_{file_type}__{sample_name}.bed"
    wildcard_constraints:
        file_type = "final|merged|pseudo1|pseudo2"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        file_type = lambda wildcards: wildcards.file_type
    log:
        temp(return_log_chip("ATAC","{sample_name}", "atac_shift_{file_type}", ""))
    conda: CONDA_ENV_ATAC
    threads: config["resources"]["atac_shift_and_bed"]["threads"]
    resources:
        mem_mb=config["resources"]["atac_shift_and_bed"]["mem_mb"],
        tmp_mb=config["resources"]["atac_shift_and_bed"]["tmp_mb"],
        qos=config["resources"]["atac_shift_and_bed"]["qos"]
    shell:
        """
        {{
        printf "\nApplying Tn5 shift and converting to BED for {params.file_type}__{params.sample_name}\n"
        if [[ ! -f "{input.bamfile}.bai" ]]; then
            samtools index -@ {threads} "{input.bamfile}"
        fi
        alignmentSieve --ATACshift -b {input.bamfile} -p {threads} --bam /dev/stdout | \
        bedtools bamtobed -i stdin > {output.bedfile}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_atac:
    input:
        bedfile = "results/ATAC/mapped/shifted_{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}.bed"
    output:
        peakfile = "results/ATAC/peaks/peaks_atac__{file_type}__{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}_peaks.narrowPeak"
    params:
        ipname = lambda wildcards: f"{wildcards.data_type}__{wildcards.line}__{wildcards.tissue}__{wildcards.sample_type}__{wildcards.replicate}__{wildcards.ref_genome}",
        filetype = lambda wildcards: wildcards.file_type,
        params = config["atac_callpeaks"]['params'],
        genomesize = lambda wildcards: config[config[wildcards.ref_genome]['species']]['genomesize']
    log:
        temp(return_log_chip("ATAC","{data_type}__{line}__{tissue}__{sample_type}__{replicate}__{ref_genome}", "{file_type}__narrowpeak_calling", "PE"))
    conda: CONDA_ENV_ATAC
    threads: config["resources"]["calling_peaks_atac"]["threads"]
    resources:
        mem_mb=config["resources"]["calling_peaks_atac"]["mem_mb"],
        tmp_mb=config["resources"]["calling_peaks_atac"]["tmp_mb"],
        qos=config["resources"]["calling_peaks_atac"]["qos"]
    shell:
        """
        {{
        printf "\nCalling narrow peaks for ATAC-seq {params.ipname} using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.bedfile} -f BED \
            -g {params.genomesize} {params.params} \
            -n peaks_atac__{params.filetype}__{params.ipname} \
            --outdir results/ATAC/peaks/
        }} 2>&1 | tee -a "{log}"
        """

rule make_coverage_atac:
    input:
        bamfile = "results/ATAC/mapped/{file_type}__{sample_name}.bam",
        bai = "results/ATAC/mapped/{file_type}__{sample_name}.bam.bai"
    output:
        bigwig = "results/ATAC/tracks/coverage__{file_type}__{sample_name}.bw"
    wildcard_constraints:
        file_type = "final|merged"
    params:
        binsize = config['atac_tracks']['binsize'],
        params = config['atac_tracks']['params']
    log:
        temp(return_log_chip("ATAC","{sample_name}", "making_coverage_{file_type}", ""))
    conda: CONDA_ENV_ATAC
    threads: config["resources"]["make_coverage_atac"]["threads"]
    resources:
        mem_mb=config["resources"]["make_coverage_atac"]["mem_mb"],
        tmp_mb=config["resources"]["make_coverage_atac"]["tmp_mb"],
        qos=config["resources"]["make_coverage_atac"]["qos"]
    shell:
        """
        {{
        printf "\nMaking coverage bigwig for ATAC-seq\n"
        bamCoverage -b {input.bamfile} -o {output.bigwig} -bs {params.binsize} -p {threads} {params.params}
        }} 2>&1 | tee -a "{log}"
        """

rule all_atac:
    input:
        final = lambda wildcards: define_final_atac_output(wildcards.ref_genome)
    output:
        touch = "results/ATAC/chkpts/ATAC_analysis__{analysis_name}__{ref_genome}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """
