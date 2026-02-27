CONDA_ENV_CHIP=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_chip.yaml")
CONDA_ENV_IDR=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_idr.yaml")

def return_log_chip(env, sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"results",env,"logs",f"tmp__{sample_name}__{step}__{paired}.log")

def get_bt2_option(env):
    """Select mapping strategy option based on environment."""
    if env == "ATAC":
        return config['atac_mapping_option']
    return config['chip_mapping_option']

def get_peaktype_for_env(sample_name_or_assay, env):
    """Return peaktype using the appropriate config for the environment.

    Accepts an Assay value (e.g. "ChIP_broad"), a sample_name/Sample_ID,
    or a legacy sample_type (e.g. "H3K9me2") -- in the last case the Assay
    is resolved via _resolve_assay.
    """
    if env == "ATAC":
        return config["atac_callpeaks"]["peaktype"]
    assay = _resolve_assay(sample_name_or_assay)
    return ASSAY_TO_PEAKTYPE.get(assay, "broad")

def _resolve_assay(sample_name_or_assay):
    """Resolve a sample_name or Assay value to its Assay.

    If the argument is already a valid Assay in ASSAY_TO_PEAKTYPE, return it.
    Otherwise look it up in samples or analysis_samples.
    """
    if sample_name_or_assay in ASSAY_TO_PEAKTYPE:
        return sample_name_or_assay
    match = samples.loc[samples["sample_name"] == sample_name_or_assay]
    if not match.empty:
        return match.iloc[0]["Assay"]
    match2 = analysis_samples.loc[analysis_samples["sample_name"] == sample_name_or_assay]
    if not match2.empty:
        return match2.iloc[0]["Assay"]
    return sample_name_or_assay

def assign_mapping_paired(wildcards, rulename, outputfile):
    sname = wildcards.sample_name
    env = get_sample_info_from_name(sname, samples, 'env')
    paired = get_sample_info_from_name(sname, samples, 'paired')
    aligned_bams = config['aligned_bams']
    if paired == "PE":
        rule_obj = getattr(rules, f"{rulename}_pe")
    elif paired == "SE":
        rule_obj = getattr(rules, f"{rulename}_se")

    return getattr(rule_obj.output, outputfile).format(sample_name=sname, env=env)

def assign_bam_file(wildcards):
    sname = wildcards.sample_name
    env = get_sample_info_from_name(sname, samples, 'env')
    aligned_bams = config['aligned_bams']
    new_bam = assign_mapping_paired(wildcards, "filter_chip", "bamfile")
    if aligned_bams:
        return f"results/{env}/mapped/copied__{sname}.bam"
    else:
        return new_bam

def assign_chip_input(wildcards):
    """Get the control sample name for a ChIP sample.

    For merged/pseudo file_types, finds any replicate's control.
    For per-replicate file_types, uses direct Control column lookup.
    """
    sname = wildcards.sample_name
    if hasattr(wildcards, 'file_type') and wildcards.file_type in ['merged', 'pseudo1', 'pseudo2']:
        # For merged/pseudo peak calling, find any replicate's control
        # The sample_name here is an analysis_name; look up replicate Sample_IDs
        rep_sids = get_replicate_sample_ids(sname, samples)
        if rep_sids:
            for sid in rep_sids:
                ctrl = get_control_sample_id(sid, samples)
                if ctrl:
                    return ctrl
        # Fallback: try direct lookup in case it's a per-replicate name
        ctrl = get_control_sample_id(sname, samples)
        if ctrl:
            return ctrl
        raise ValueError(f"No control found for {sname}")
    else:
        ctrl = get_control_sample_id(sname, samples)
        if ctrl is None:
            raise ValueError(f"No control found for sample '{sname}'")
        return ctrl

def assign_peak_files_for_idr(wildcards):
    """Return list of per-replicate peak files for IDR analysis."""
    sname = wildcards.sample_name  # analysis_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)
    rep_sids = get_replicate_sample_ids(sname, samples)
    if env == "ATAC":
        prefix = "peaks_atac"
    elif paired == "PE":
        prefix = "peaks_pe"
    else:
        prefix = "peaks_se"
    return [f"results/{env}/peaks/{prefix}__final__{sid}_peaks.{peaktype}Peak"
            for sid in rep_sids]

def input_peak_files_for_best_peaks(wildcards):
    """Return peak files for selecting best peaks via pseudorep consistency."""
    sname = wildcards.sample_name  # analysis_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)

    if env == "ATAC":
        prefix = "peaks_atac"
    elif paired == "PE":
        prefix = "peaks_pe"
    else:
        prefix = "peaks_se"

    rep_sids = get_replicate_sample_ids(sname, samples)

    if len(rep_sids) >= 2:
        result = [f"results/{env}/peaks/{prefix}__merged__{sname}_peaks.{peaktype}Peak",
                  f"results/{env}/peaks/{prefix}__pseudo1__{sname}_peaks.{peaktype}Peak",
                  f"results/{env}/peaks/{prefix}__pseudo2__{sname}_peaks.{peaktype}Peak",
                  f"results/{env}/peaks/idr_peaks__{sname}.bed"]
    else:
        one_sid = rep_sids[0]
        result = [f"results/{env}/peaks/{prefix}__final__{one_sid}_peaks.{peaktype}Peak",
                  f"results/{env}/peaks/{prefix}__pseudo1__{one_sid}_peaks.{peaktype}Peak",
                  f"results/{env}/peaks/{prefix}__pseudo2__{one_sid}_peaks.{peaktype}Peak",
                  "results/empty.txt"]

    return result

def get_replicate_name(wildcards, pos):
    """Return the peak file path for replicate at position `pos`."""
    sname = wildcards.sample_name  # analysis_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    if env == "ATAC":
        prefix = "peaks_atac"
    elif paired == "PE":
        prefix = "peaks_pe"
    else:
        prefix = "peaks_se"
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)
    rep_sids = get_replicate_sample_ids(sname, samples)

    if pos >= len(rep_sids):
        return "missingrep"
    else:
        return f"results/{env}/peaks/{prefix}__final__{rep_sids[pos]}_peaks.{peaktype}Peak"

def get_replicate_pairs(wildcards):
    """Return colon-separated pairs of replicate Sample_IDs for IDR."""
    sname = wildcards.sample_name  # analysis_name
    rep_sids = get_replicate_sample_ids(sname, samples)
    pairs = []
    for i in range(0, len(rep_sids)):
        for j in range(i + 1, len(rep_sids)):
            pairs.append(f"{rep_sids[i]}:{rep_sids[j]}")
    return pairs

def define_chipseq_target_file(wildcards):
    tarname = config['motif_target_file_label']
    env = wildcards.env
    peak_file = wildcards.peak_file
    parts = peak_file.split("__")
    file_type = parts[0]
    if file_type == "selected_peaks":
        # analysis_name is everything after "selected_peaks__"
        spname = "__".join(parts[1:])
        inputfile = f"results/{env}/peaks/{file_type}__{spname}.bedPeak"
        ref_genome = get_sample_info_from_name(spname, analysis_samples, 'ref_genome')
        fasta = f"genomes/{ref_genome}/{ref_genome}.fa"
        if any(analysis_samples['sample_name'] == spname):
            return [inputfile, fasta]
    elif file_type == "idr_peaks":
        spname = "__".join(parts[1:])
        inputfile = f"results/{env}/peaks/{file_type}__{spname}.bed"
        ref_genome = get_sample_info_from_name(spname, analysis_samples, 'ref_genome')
        fasta = f"genomes/{ref_genome}/{ref_genome}.fa"
        if any(analysis_samples['sample_name'] == spname):
            return [inputfile, fasta]
    elif file_type.startswith("peaks_"):
        # Format: peaks_{pe|se|atac}__{file_cat}__{sample_name}_peaks
        # parts[0] = "peaks_pe" (or se/atac), parts[1] = file_cat (final/merged/pseudo),
        # rest = sample_name + "_peaks" suffix
        file_cat = parts[1]
        # Rejoin remaining parts and strip "_peaks" suffix
        rest = "__".join(parts[2:])
        if rest.endswith("_peaks"):
            sname = rest[:-len("_peaks")]
        else:
            sname = rest
        assay = _resolve_assay(sname)
        peaktype = get_peaktype_for_env(assay, env)
        inputfile = f"results/{env}/peaks/{file_type}__{file_cat}__{sname}_peaks.{peaktype}Peak"
        parsed = parse_sample_name(sname)
        ref_genome = parsed['ref_genome']
        fasta = f"genomes/{ref_genome}/{ref_genome}.fa"
        if any(samples['sample_name'] == sname) or any(analysis_samples['sample_name'] == sname):
            return [inputfile, fasta]
    elif peak_file == tarname:
        ref_genome = config['motif_ref_genome']
        fasta = f"genomes/{ref_genome}/{ref_genome}.fa"
        inputfile = config['motif_target_file']
        return [inputfile, fasta]
    else:
        return ValueError(
            f"{wildcards.peak_file} is unknown."
            "Options are either peakfiles generated by the pipeline"
            "or the value of 'motifs_target_file_name' in the config file"
            )

def define_input_manorm(wildcards, string):
    """Resolve inputs for pairwise differential peak analysis (MAnorm).

    sample1 and sample2 are analysis_names (wildcards).
    """
    sample1 = wildcards.sample1
    env1 = get_sample_info_from_name(sample1, analysis_samples, 'env')
    paired1 = get_sample_info_from_name(sample1, analysis_samples, 'paired')

    sample2 = wildcards.sample2
    env2 = get_sample_info_from_name(sample2, analysis_samples, 'env')
    paired2 = get_sample_info_from_name(sample2, analysis_samples, 'paired')

    if env1 == "ChIP":
        params = config['diffpeaks_params']['chip_pe'] if paired1 == "PE" and paired2 == "PE" else config['diffpeaks_params']['chip_se']
    elif env1 == "ATAC":
        params = config['diffpeaks_params']['ATAC_pe'] if paired1 == "PE" and paired2 == "PE" else config['diffpeaks_params']['ATAC_se']
    else:
        # Fallback for any other env
        params = config['diffpeaks_params'].get('chip_se', '')

    assay1 = _resolve_assay(sample1)
    assay2 = _resolve_assay(sample2)
    peaktype1 = get_peaktype_for_env(assay1, env1)
    peaktype2 = get_peaktype_for_env(assay2, env2)
    if peaktype1 != peaktype2:
        raise ValueError(f"{sample1} and {sample2} have different peaktypes.")
    else:
        peaktype = peaktype1

    peakfile1 = f"results/{env1}/peaks/selected_peaks__{sample1}.bedPeak"
    peakfile2 = f"results/{env2}/peaks/selected_peaks__{sample2}.bedPeak"

    rep_sids1 = get_replicate_sample_ids(sample1, samples)
    add1 = "shifted_" if env1 == "ATAC" else ""
    if len(rep_sids1) >= 2:
        bamfile1 = f"results/{env1}/mapped/{add1}merged__{sample1}.bam"
    else:
        bamfile1 = f"results/{env1}/mapped/{add1}final__{rep_sids1[0]}.bam"

    rep_sids2 = get_replicate_sample_ids(sample2, samples)
    add2 = "shifted_" if env2 == "ATAC" else ""
    if len(rep_sids2) >= 2:
        bamfile2 = f"results/{env2}/mapped/{add2}merged__{sample2}.bam"
    else:
        bamfile2 = f"results/{env2}/mapped/{add2}final__{rep_sids2[0]}.bam"

    if string == "peaks1":
        return peakfile1
    elif string == "peaks2":
        return peakfile2
    elif string == "reads1":
        return bamfile1
    elif string == "reads2":
        return bamfile2
    elif string == "params":
        return params
    elif string == "format":
        return peaktype

def define_logs_final_input(wildcards):
    """Collect all log files needed for peak stats summary."""
    log_files = []
    sname = wildcards.sample_name  # analysis_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)
    rep_sids = get_replicate_sample_ids(sname, samples)

    for sid in rep_sids:
        log_files.append(return_log_chip(env, sid, f"final__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, sid, "making_bigwig_final", ""))
        if env != "ATAC":
            log_files.append(return_log_chip(env, sid, "making_fingerprint_final", ""))

    if len(rep_sids) >= 2:
        log_files.append(return_log_chip(env, sname, "IDR", ""))
        log_files.append(return_log_chip(env, sname, "merging_reps", ""))
        log_files.append(return_log_chip(env, sname, f"merged__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, sname, f"pseudo1__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, sname, f"pseudo2__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, sname, "splitting_pseudreps", ""))
    else:
        one_sid = rep_sids[0]
        log_files.append(return_log_chip(env, one_sid, "splitting_pseudreps", ""))
        log_files.append(return_log_chip(env, one_sid, f"pseudo1__{peaktype}peak_calling", paired))
        log_files.append(return_log_chip(env, one_sid, f"pseudo2__{peaktype}peak_calling", paired))

    log_files.append(return_log_chip(env, sname, "selecting_best_peaks", ""))

    return log_files

def define_final_chip_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    motifs = config['motifs']
    motifs_allreps = config['motifs_allreps']
    trimmed_fastqs = config['trimmed_fastqs']
    aligned_bams = config['aligned_bams']
    map_files = []
    stat_files = []
    qc_files = []
    peak_files = []
    bigwig_files = []
    motif_files = []
    allrep_files = []

    # All ChIP samples (including controls) for mapping/QC
    controls = identify_control_samples(samples)
    filtered_rep_samples = samples[(samples['env'] == 'ChIP') & (samples['ref_genome'] == ref_genome)].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = row['sample_name']
        paired = row['paired']
        env = row['env']
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

    # Non-control samples for peak calling and bigwigs
    filtered_rep_samples_no_input = filtered_rep_samples[~filtered_rep_samples['Sample_ID'].isin(controls)].copy()
    for _, row in filtered_rep_samples_no_input.iterrows():
        assay = row['Assay']
        peaktype = ASSAY_TO_PEAKTYPE.get(assay, "broad")
        sname = row['sample_name']
        paired = row['paired']
        env = row['env']
        bigwig_files.append(f"results/{env}/tracks/FC__final__{sname}.bw")
        stat_files.append(f"results/{env}/plots/Fingerprint__final__{sname}.png")
        if paired == "PE":
            peak_files.append(f"results/{env}/peaks/peaks_pe__final__{sname}_peaks.{peaktype}Peak")
        else:
            peak_files.append(f"results/{env}/peaks/peaks_se__final__{sname}_peaks.{peaktype}Peak")

    # Analysis-level outputs (merged, IDR, best peaks)
    filtered_analysis_samples = analysis_samples[(analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome)].copy()
    for _, row in filtered_analysis_samples.iterrows():
        aname = row['sample_name']  # analysis_name
        env = row['env']
        peak_files.append(f"results/{env}/peaks/selected_peaks__{aname}.bedPeak")
        rep_sids = get_replicate_sample_ids(aname, samples)
        if len(rep_sids) >= 2:
            bigwig_files.append(f"results/{env}/tracks/FC__merged__{aname}.bw")
            stat_files.append(f"results/{env}/chkpts/idr__{aname}.done")

    # Motif analysis for TF (now subsumed under ChIP env - only if data contains TF-like samples)
    # Note: TF env is eliminated; motif_files only apply if explicitly configured
    # Pairwise differential peaks (MAnorm) for ChIP
    chip_analysis_samples = analysis_samples[(analysis_samples["env"] == "ChIP") & (analysis_samples["ref_genome"] == ref_genome)].copy()
    for ip_target, group in chip_analysis_samples.groupby("sample_type"):
        if len(group) >= 2:
            for a, b in combinations(group.itertuples(index=False), 2):
                sample1 = a.sample_name
                sample2 = b.sample_name
                assay1 = a.Assay
                assay2 = b.Assay
                peaktype1 = ASSAY_TO_PEAKTYPE.get(assay1, "broad")
                peaktype2 = ASSAY_TO_PEAKTYPE.get(assay2, "broad")
                if peaktype1 == peaktype2:
                    peak_files.append(f"results/ChIP/peaks/{sample1}_vs_{sample2}/{sample1}_vs_{sample2}_all_MAvalues.xls")

    results = map_files + bigwig_files

    if qc_option == "all":
        results += qc_files

    if analysis:
        results += peak_files + stat_files

    if motifs:
        results += motif_files

    if motifs_allreps:
        results += allrep_files

    return results

rule make_bt2_indices:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    output:
        indices = directory("genomes/{ref_genome}/bt2_index")
    log:
        temp(os.path.join(REPO_FOLDER,"results","logs","bowtie_index_{ref_genome}.log"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["make_bt2_indices"]["threads"]
    resources:
        mem_mb=config["resources"]["make_bt2_indices"]["mem_mb"],
        tmp_mb=config["resources"]["make_bt2_indices"]["tmp_mb"],
        qos=config["resources"]["make_bt2_indices"]["qos"]
    shell:
        """
        {{
        printf "\nBuilding Bowtie2 index for {wildcards.ref_genome}\n"
        mkdir genomes/{wildcards.ref_genome}/bt2_index
        bowtie2-build --threads {threads} "{input.fasta}" "{output.indices}/{wildcards.ref_genome}"
        }} 2>&1 | tee -a "{log}"
        """

rule bowtie2_map_pe:
    input:
        fastq1 = "results/{env}/fastq/trim__{sample_name}__R1.fastq.gz",
        fastq2 = "results/{env}/fastq/trim__{sample_name}__R2.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/bt2_index"
    output:
        samfile = temp("results/{env}/mapped/mapped_pe__{sample_name}.sam"),
        metrics = "results/{env}/reports/bt2_pe__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: get_bt2_option(wildcards.env),
        mapping_params = lambda wildcards: config['bt2_mapping_strategy'][get_bt2_option(wildcards.env)]['map_pe']
    log:
        temp(return_log_chip("{env}","{sample_name}", "mappingBT2", "PE"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["bowtie2_map_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["bowtie2_map_pe"]["mem_mb"],
        tmp_mb=config["resources"]["bowtie2_map_pe"]["tmp_mb"],
        qos=config["resources"]["bowtie2_map_pe"]["qos"]
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} -x "{input.indices}/{params.ref_genome}" -1 "{input.fastq1}" -2 "{input.fastq2}" -S "{output.samfile}" 2>&1 | tee "{output.metrics}"
        }} 2>&1 | tee -a "{log}"
        """

rule bowtie2_map_se:
    input:
        fastq = "results/{env}/fastq/trim__{sample_name}__R0.fastq.gz",
        indices = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/bt2_index"
    output:
        samfile = temp("results/{env}/mapped/mapped_se__{sample_name}.sam"),
        metrics = "results/{env}/reports/bt2_se__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        map_option = lambda wildcards: get_bt2_option(wildcards.env),
        mapping_params = lambda wildcards: config['bt2_mapping_strategy'][get_bt2_option(wildcards.env)]['map_se']
    log:
        temp(return_log_chip("{env}","{sample_name}", "mappingBT2", "SE"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["bowtie2_map_se"]["threads"]
    resources:
        mem_mb=config["resources"]["bowtie2_map_se"]["mem_mb"],
        tmp_mb=config["resources"]["bowtie2_map_se"]["tmp_mb"],
        qos=config["resources"]["bowtie2_map_se"]["qos"]
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with {params.map_option} parameters with bowtie2 version:\n"
		bowtie2 --version
		bowtie2 -p {threads} {params.mapping_params} -x "{input.indices}/{params.ref_genome}" -U "{input.fastq}" -S "{output.samfile}" 2>&1 | tee "{output.metrics}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_chip_pe:
    input:
        samfile = "results/{env}/mapped/mapped_pe__{sample_name}.sam"
    output:
        bamfile = temp("results/{env}/mapped/mapped_pe__{sample_name}.bam"),
        metrics_dup = "results/{env}/reports/markdup_pe__{sample_name}.txt",
        metrics_flag = "results/{env}/reports/flagstat_pe__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        map_option = lambda wildcards: get_bt2_option(wildcards.env),
        filtering_params = lambda wildcards: config['bt2_mapping_strategy'][get_bt2_option(wildcards.env)]['filter']
    log:
        temp(return_log_chip("{env}","{sample_name}", "filteringChIP", "PE"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["filter_chip_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["filter_chip_pe"]["mem_mb"],
        tmp_mb=config["resources"]["filter_chip_pe"]["tmp_mb"],
        qos=config["resources"]["filter_chip_pe"]["qos"]
    shell:
        """
        {{
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {params.sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o "results/{params.env}/mapped/temp1_{params.sample_name}.bam" "{input.samfile}"
        samtools fixmate -@ {threads} -m "results/{params.env}/mapped/temp1_{params.sample_name}.bam" "results/{params.env}/mapped/temp2_{params.sample_name}.bam"
        samtools sort -@ {threads} -o "results/{params.env}/mapped/temp3_{params.sample_name}.bam" "results/{params.env}/mapped/temp2_{params.sample_name}.bam"
        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} "results/{params.env}/mapped/temp3_{params.sample_name}.bam" "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f results/{params.env}/mapped/temp*"_{params.sample_name}.bam"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_chip_se:
    input:
        samfile = "results/{env}/mapped/mapped_se__{sample_name}.sam"
    output:
        bamfile = temp("results/{env}/mapped/mapped_se__{sample_name}.bam"),
        metrics_dup = "results/{env}/reports/markdup_se__{sample_name}.txt",
        metrics_flag = "results/{env}/reports/flagstat_se__{sample_name}.txt"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        map_option = lambda wildcards: get_bt2_option(wildcards.env),
        filtering_params = lambda wildcards: config['bt2_mapping_strategy'][get_bt2_option(wildcards.env)]['filter']
    log:
        temp(return_log_chip("{env}","{sample_name}", "filteringChIP", "SE"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["filter_chip_se"]["threads"]
    resources:
        mem_mb=config["resources"]["filter_chip_se"]["mem_mb"],
        tmp_mb=config["resources"]["filter_chip_se"]["tmp_mb"],
        qos=config["resources"]["filter_chip_se"]["qos"]
    shell:
        """
        {{
        printf "\nRemoving low quality reads, secondary alignements and duplicates, sorting and indexing {params.sample_name} file using {params.map_option} with samtools version:\n"
        samtools --version
        samtools view -@ {threads} -b -h -q 10 -F 256 -o "results/{params.env}/mapped/temp1_{params.sample_name}.bam" "{input.samfile}"
        samtools sort -@ {threads} -o "results/{params.env}/mapped/temp2_{params.sample_name}.bam" "results/{params.env}/mapped/temp1_{params.sample_name}.bam"
        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} "results/{params.env}/mapped/temp2_{params.sample_name}.bam" "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f results/{params.env}/mapped/temp*"_{params.sample_name}.bam"
        }} 2>&1 | tee -a "{log}"
        """

rule make_chip_stats_pe:
    input:
        metrics_trim = "results/{env}/reports/trim_pe__{sample_name}.txt",
        metrics_map = "results/{env}/reports/bt2_pe__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.env, wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingBT2", "filteringChIP"] ]
    output:
        stat_file = "results/{env}/reports/summary_{env}_PE_mapping_stats_{sample_name}.txt",
        log = "results/{env}/logs/process_chip_pe_sample__{sample_name}.log"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        line = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'line'),
        tissue = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'tissue'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'sample_type'),
        replicate = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'replicate'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'ref_genome'),
        trimmed_fastq = config['trimmed_fastqs']
    threads: config["resources"]["make_chip_stats_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["make_chip_stats_pe"]["mem_mb"],
        tmp_mb=config["resources"]["make_chip_stats_pe"]["tmp_mb"],
        qos=config["resources"]["make_chip_stats_pe"]["qos"]
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(grep "Total read pairs processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        else
            tot=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        fi
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned concordantly >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned concordantly exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule make_chip_stats_se:
    input:
        metrics_trim = "results/{env}/reports/trim_se__{sample_name}.txt",
        metrics_map = "results/{env}/reports/bt2_se__{sample_name}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.env, wildcards.sample_name, step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingBT2", "filteringChIP"] ]
    output:
        stat_file = "results/{env}/reports/summary_{env}_SE_mapping_stats_{sample_name}.txt",
        log = "results/{env}/logs/process_chip_se_sample__{sample_name}.log"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        line = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'line'),
        tissue = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'tissue'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'sample_type'),
        replicate = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'replicate'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'ref_genome'),
        trimmed_fastq = config['trimmed_fastqs']
    threads: config["resources"]["make_chip_stats_se"]["threads"]
    resources:
        mem_mb=config["resources"]["make_chip_stats_se"]["mem_mb"],
        tmp_mb=config["resources"]["make_chip_stats_se"]["tmp_mb"],
        qos=config["resources"]["make_chip_stats_se"]["qos"]
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(grep "Total reads processed:" "{input.metrics_trim}" | awk '{{print $NF}}' | sed 's/,//g')
        else
            tot=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        fi
        filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        multi=$(grep "aligned >1 times" "{input.metrics_map}" | awk '{{print $1}}')
        single=$(grep "aligned exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        allmap=$((multi+single))
        printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        """

rule pe_or_se_chip_dispatch:
    input:
        assign_bam_file
    output:
        bam = "results/{env}/mapped/final__{sample_name}.bam",
        bai = "results/{env}/mapped/final__{sample_name}.bam.bai",
        touch = "results/{env}/chkpts/map_{env}__{sample_name}.done"
    wildcard_constraints:
        env = "ChIP|ATAC"
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["pe_or_se_chip_dispatch"]["threads"]
    resources:
        mem_mb=config["resources"]["pe_or_se_chip_dispatch"]["mem_mb"],
        tmp_mb=config["resources"]["pe_or_se_chip_dispatch"]["tmp_mb"],
        qos=config["resources"]["pe_or_se_chip_dispatch"]["qos"]
    shell:
        """
        cp {input} {output.bam}
        samtools index -@ {threads} "{output.bam}"
        touch {output.touch}
        """

rule make_coverage_chip:
    input:
        bamfile = "results/{env}/mapped/final__{sample_name}.bam"
    output:
        bigwigcov = "results/{env}/tracks/coverage__{sample_name}.bw"
    wildcard_constraints:
        env = "ChIP"
    params:
        binsize = config['chip_tracks']['binsize']
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["make_coverage_chip"]["threads"]
    resources:
        mem_mb=config["resources"]["make_coverage_chip"]["mem_mb"],
        tmp_mb=config["resources"]["make_coverage_chip"]["tmp_mb"],
        qos=config["resources"]["make_coverage_chip"]["qos"]
    shell:
        """
        bamCoverage -b {input.bamfile} -o {output.bigwigcov} -bs {params.binsize} -p {threads}
        """

rule make_bigwig_chip:
    input:
        ipfile = "results/{env}/mapped/{file_type}__{sample_name}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        bigwigfile = "results/{env}/tracks/FC__{file_type}__{sample_name}.bw"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged"
    params:
        ipname = lambda wildcards: f"{wildcards.file_type}__{wildcards.sample_name}",
        inputname = lambda wildcards: f"{wildcards.file_type}__{assign_chip_input(wildcards)}",
        binsize = config['chip_tracks']['binsize'],
        params = config['chip_tracks']['params']
    log:
        temp(return_log_chip("{env}","{sample_name}", "making_bigwig_{file_type}", ""))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["make_bigwig_chip"]["threads"]
    resources:
        mem_mb=config["resources"]["make_bigwig_chip"]["mem_mb"],
        tmp_mb=config["resources"]["make_bigwig_chip"]["tmp_mb"],
        qos=config["resources"]["make_bigwig_chip"]["qos"]
    shell:
        """
        {{
        printf "\nMaking bigwig files for {params.ipname} (vs {params.inputname}) with deeptools version:\n"
        deeptools --version
        bamCompare -b1 {input.ipfile} -b2 {input.inputfile} -o {output.bigwigfile} -p {threads} --binSize {params.binsize} {params.params}
        }} 2>&1 | tee -a "{log}"
        """

rule make_fingerprint_plot:
    input:
        ipfile = "results/{env}/mapped/{file_type}__{sample_name}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        pngplot = "results/{env}/plots/Fingerprint__{file_type}__{sample_name}.png"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged"
    params:
        ipname = lambda wildcards: f"{wildcards.file_type}__{wildcards.sample_name}",
        inputname = lambda wildcards: f"{wildcards.file_type}__{assign_chip_input(wildcards)}"
    log:
        temp(return_log_chip("{env}","{sample_name}", "making_fingerprint_{file_type}", ""))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["make_fingerprint_plot"]["threads"]
    resources:
        mem_mb=config["resources"]["make_fingerprint_plot"]["mem_mb"],
        tmp_mb=config["resources"]["make_fingerprint_plot"]["tmp_mb"],
        qos=config["resources"]["make_fingerprint_plot"]["qos"]
    shell:
        """
        {{
        printf "\nPlotting fingerprint for {params.ipname} (vs {params.inputname}) with deeptools version:\n"
        deeptools --version
        plotFingerprint -b {input.ipfile} {input.inputfile} -o {output.pngplot} -p {threads} -l {params.ipname} {params.inputname}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_macs2_pe:
    input:
        ipfile = "results/{env}/mapped/{file_type}__{sample_name}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        peakfile = "results/{env}/peaks/peaks_pe__{file_type}__{sample_name}_peaks.{peaktype}Peak"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged|pseudo1|pseudo2"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        inputname = lambda wildcards: assign_chip_input(wildcards),
        peaktype = lambda wildcards: ASSAY_TO_PEAKTYPE.get(_resolve_assay(wildcards.sample_name), "broad"),
        file_type = lambda wildcards: wildcards.file_type,
        env = lambda wildcards: wildcards.env,
        params = config["chip_callpeaks"]['params'],
        genomesize = lambda wildcards: config[config[parse_sample_name(wildcards.sample_name)['ref_genome']]['species']]['genomesize']
    log:
        temp(return_log_chip("{env}","{sample_name}", "{file_type}__{peaktype}peak_calling", "PE"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["calling_peaks_macs2_pe"]["threads"]
    resources:
        mem_mb=config["resources"]["calling_peaks_macs2_pe"]["mem_mb"],
        tmp_mb=config["resources"]["calling_peaks_macs2_pe"]["tmp_mb"],
        qos=config["resources"]["calling_peaks_macs2_pe"]["qos"]
    shell:
        """
        {{
        if [[ "{params.peaktype}" == "broad" ]]; then
            add="--broad"
        else
            add=""
        fi
        printf "\nCalling {params.peaktype} peaks for paired-end {params.sample_name} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAMPE -g {params.genomesize} {params.params} -n peaks_pe__{params.file_type}__{params.sample_name} --outdir results/{params.env}/peaks/ ${{add}}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_macs2_se:
    input:
        ipfile = "results/{env}/mapped/{file_type}__{sample_name}.bam",
        inputfile = lambda wildcards: f"results/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam"
    output:
        peakfile = "results/{env}/peaks/peaks_se__{file_type}__{sample_name}_peaks.{peaktype}Peak"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged|pseudo1|pseudo2"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        inputname = lambda wildcards: assign_chip_input(wildcards),
        peaktype = lambda wildcards: ASSAY_TO_PEAKTYPE.get(_resolve_assay(wildcards.sample_name), "broad"),
        file_type = lambda wildcards: wildcards.file_type,
        env = lambda wildcards: wildcards.env,
        params = config["chip_callpeaks"]['params'],
        genomesize = lambda wildcards: config[config[parse_sample_name(wildcards.sample_name)['ref_genome']]['species']]['genomesize']
    log:
        temp(return_log_chip("{env}","{sample_name}", "{file_type}__{peaktype}peak_calling", "SE"))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["calling_peaks_macs2_se"]["threads"]
    resources:
        mem_mb=config["resources"]["calling_peaks_macs2_se"]["mem_mb"],
        tmp_mb=config["resources"]["calling_peaks_macs2_se"]["tmp_mb"],
        qos=config["resources"]["calling_peaks_macs2_se"]["qos"]
    shell:
        """
        {{
        if [[ "{params.peaktype}" == "broad" ]]; then
            add="--broad"
        else
            add=""
        fi
        printf "\nCalling {params.peaktype} peaks for single-end {params.sample_name} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAM -g {params.genomesize} {params.params} -n peaks_se__{params.file_type}__{params.sample_name} --outdir results/{params.env}/peaks/ ${{add}}
        }} 2>&1 | tee -a "{log}"
        """

rule idr_analysis_replicates:
    input:
        peak_file = lambda wildcards: assign_peak_files_for_idr(wildcards)
    output:
        touch = "results/{env}/chkpts/idr__{sample_name}.done",
        merged = "results/{env}/peaks/idr_peaks__{sample_name}.bed"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        peaktype = lambda wildcards: get_peaktype_for_env(_resolve_assay(wildcards.sample_name), wildcards.env),
        paired = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'paired'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome'),
        env = lambda wildcards: wildcards.env,
        replicate_pairs = lambda wildcards: get_replicate_pairs(wildcards),
        rep_sids = lambda wildcards: get_replicate_sample_ids(wildcards.sample_name, samples)
    log:
        temp(return_log_chip("{env}","{sample_name}", "IDR", ""))
    conda: CONDA_ENV_IDR
    threads: config["resources"]["idr_analysis_replicates"]["threads"]
    resources:
        mem_mb=config["resources"]["idr_analysis_replicates"]["mem_mb"],
        tmp_mb=config["resources"]["idr_analysis_replicates"]["tmp_mb"],
        qos=config["resources"]["idr_analysis_replicates"]["qos"]
    shell:
        """
        {{
        printf "\nLooping over each unique pair of biological replicates for {params.sname} to perform IDR with:\n"
        idr --version
		if [[ "{params.env}" == "ATAC" ]]; then
            pre="atac"
        elif [[ "{params.paired}" == "PE" ]]; then
            pre="pe"
        else
            pre="se"
        fi
        temp="results/{params.env}/peaks/temp_idr_peaks__{params.sname}.bed"
        while read chr max; do
            printf "${{chr}}\t1\t${{max}}\n" >> "${{temp}}"
        done < "genomes/{params.ref_genome}/chrom.sizes"
        mkdir -p results/{params.env}/plots/
        for pair in {params.replicate_pairs}; do
            rep1_sid=$(echo ${{pair}} | cut -d":" -f1)
            rep2_sid=$(echo ${{pair}} | cut -d":" -f2)
            file1="results/{params.env}/peaks/peaks_${{pre}}__final__${{rep1_sid}}_peaks.{params.peaktype}Peak"
            file2="results/{params.env}/peaks/peaks_${{pre}}__final__${{rep2_sid}}_peaks.{params.peaktype}Peak"
            outfile="results/{params.env}/peaks/idr_${{pre}}__{params.sname}__${{rep1_sid}}_vs_${{rep2_sid}}_peaks.{params.peaktype}Peak"
            printf "\nPerforming IDR for ${{rep1_sid}} vs ${{rep2_sid}}\n"
            idr --input-file-type {params.peaktype}Peak --output-file-type {params.peaktype}Peak --samples ${{file1}} ${{file2}} -o ${{outfile}} -l results/{params.env}/reports/idr_{params.sname}.log --plot || true
            ## "|| true" avoids pipeline-breaking errors if no positive peaks were found
            mv "${{outfile}}.png" results/{params.env}/plots/
            filtered="${{outfile}}.filtered"
            awk -v OFS="\t" '$5>=540 {{print $1,$2,$3}}' ${{outfile}} | sort -k1,1 -k2,2n > "${{filtered}}"
            new="${{temp}}.new"
            bedtools intersect -a ${{temp}} -b ${{filtered}} > "${{new}}"
            mv "${{new}}" "${{temp}}"
        done
        cat ${{temp}} > {output.merged}
        bedtools intersect -a ${{file2}} -b ${{temp}} -u > "results/{params.env}/peaks/idr_peaks__{params.sname}.{params.peaktype}Peak"
        rm -f ${{temp}} ${{new}} ${{filtered}}
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule merging_chip_replicates:
    input:
        bamfiles = lambda wildcards: [
            f"results/{wildcards.env}/mapped/final__{sid}.bam"
            for sid in get_replicate_sample_ids(wildcards.sample_name, samples)
        ]
    output:
        temp_merge = temp("results/{env}/mapped/temp_merged__{sample_name}.bam"),
        mergefile = "results/{env}/mapped/merged__{sample_name}.bam",
        mergebai = "results/{env}/mapped/merged__{sample_name}.bam.bai"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env
    log:
        temp(return_log_chip("{env}","{sample_name}", "merging_reps", ""))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["merging_chip_replicates"]["threads"]
    resources:
        mem_mb=config["resources"]["merging_chip_replicates"]["mem_mb"],
        tmp_mb=config["resources"]["merging_chip_replicates"]["tmp_mb"],
        qos=config["resources"]["merging_chip_replicates"]["qos"]
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
		samtools merge -@ {threads} {output.temp_merge} {input.bamfiles}
		samtools sort -@ {threads} -o {output.mergefile} {output.temp_merge}
		samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """

rule making_pseudo_replicates:
    input:
        bamfile = lambda wildcards: f"results/{wildcards.env}/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'final'}__{wildcards.sample_name}.bam"
    output:
        temp_pseudo1 = temp("results/{env}/mapped/temp_pseudo1__{sample_name}.bam"),
        temp_pseudo2 = temp("results/{env}/mapped/temp_pseudo2__{sample_name}.bam"),
        pseudo1 = temp("results/{env}/mapped/pseudo1__{sample_name}.bam"),
        pseudo2 = temp("results/{env}/mapped/pseudo2__{sample_name}.bam")
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env
    log:
        temp(return_log_chip("{env}","{sample_name}", "splitting_pseudreps", ""))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["making_pseudo_replicates"]["threads"]
    resources:
        mem_mb=config["resources"]["making_pseudo_replicates"]["mem_mb"],
        tmp_mb=config["resources"]["making_pseudo_replicates"]["tmp_mb"],
        qos=config["resources"]["making_pseudo_replicates"]["qos"]
    shell:
        """
        {{
        printf "\nSplitting {params.sname} in two pseudo-replicates\n"
        samtools view -b -h -s 1.5 -@ {threads} -U {output.temp_pseudo2} -o {output.temp_pseudo1} {input.bamfile}
		samtools sort -@ {threads} -o {output.pseudo1} {output.temp_pseudo1}
		samtools sort -@ {threads} -o {output.pseudo2} {output.temp_pseudo2}
        }} 2>&1 | tee -a "{log}"
        """

rule create_empty_file:
    output:
        "results/empty.txt"
    localrule: True
    shell:
        "touch {output}"

rule best_peaks_pseudoreps:
    input:
        chrom_sizes = lambda wildcards: f"genomes/{get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome')}/chrom.sizes",
        peakfiles = input_peak_files_for_best_peaks
    output:
        bestpeaks = "results/{env}/peaks/selected_peaks__{sample_name}.bedPeak",
        stats_pseudoreps = temp("results/{env}/reports/stats_pseudoreps__{sample_name}.txt")
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        peaktype = lambda wildcards: get_peaktype_for_env(_resolve_assay(wildcards.sample_name), wildcards.env)
    log:
        temp(return_log_chip("{env}","{sample_name}", "selecting_best_peaks", ""))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["best_peaks_pseudoreps"]["threads"]
    resources:
        mem_mb=config["resources"]["best_peaks_pseudoreps"]["mem_mb"],
        tmp_mb=config["resources"]["best_peaks_pseudoreps"]["tmp_mb"],
        qos=config["resources"]["best_peaks_pseudoreps"]["qos"]
    shell:
        """
        {{
        printf "\nIntersecting merged peaks ({input.peakfiles[0]}), and both pseudo replicates to select the best peaks for {params.sname}\n"
        awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[0]} | sort -k1,1 -k2,2n -u > "results/{params.env}/peaks/temp_{params.sname}_merged.bed"
		awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[1]} | sort -k1,1 -k2,2n -u > "results/{params.env}/peaks/temp_{params.sname}_pseudo1.bed"
		awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[2]} | sort -k1,1 -k2,2n -u > "results/{params.env}/peaks/temp_{params.sname}_pseudo2.bed"
		bedtools intersect -a results/{params.env}/peaks/temp_{params.sname}_pseudo1.bed -b results/{params.env}/peaks/temp_{params.sname}_pseudo2.bed > "results/{params.env}/peaks/temp_{params.sname}_pseudos.bed"
		bedtools intersect -a results/{params.env}/peaks/temp_{params.sname}_merged.bed -b results/{params.env}/peaks/temp_{params.sname}_pseudo1.bed -u > "results/{params.env}/peaks/temp_{params.sname}_selected.bed"
		bedtools intersect -a {input.peakfiles[0]} -b results/{params.env}/peaks/temp_{params.sname}_selected.bed -u > "results/{params.env}/peaks/selected_peaks__{params.sname}.{params.peaktype}Peak"
        printf "\nGetting best quality peaks peaks\n"
        ## Note: If broadpeak, an additional "summit" column will be added for potential downstream processes, which only represent the middle of the peak, not its actual summit.
        sort -k1,1 -k2,2n -k5nr results/{params.env}/peaks/selected_peaks__{params.sname}.{params.peaktype}Peak | awk -v OFS="\t" -v t={params.peaktype} '{{a=$1":"$2":"$3; if (a!=n) {{if (t=="broad") $10=int(($3-$2)/2); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}} n=$1":"$2":"$3}}' | bedtools sort -g {input.chrom_sizes} > "{output.bestpeaks}"
        printf "\nExtracting peak stats for {params.sname}\n"
        merged=$(wc -l results/{params.env}/peaks/temp_{params.sname}_merged.bed | cut -d" " -f1)
		pseudos=$(awk '{{print $1,$2,$3}}' results/{params.env}/peaks/temp_{params.sname}_pseudos.bed | sort -k1,1 -k2,2n -u | wc -l)
		selected=$(cat results/{params.env}/peaks/temp_{params.sname}_selected.bed | sort -k1,1 -k2,2n -u | wc -l)
        if [[ "{input.peakfiles[3]}" == "results/empty.txt" ]]; then
            idr="0"
        else
            idr=$(wc -l {input.peakfiles[3]} | cut -d" " -f1)
        fi
		printf "Merged=${{merged}}\nPseudos=${{pseudos}}\nIDR=${{idr}}\nSelected=${{selected}}\n" > "{output.stats_pseudoreps}"
        rm -f "results/{params.env}/peaks/temp_{params.sname}"*
        }} 2>&1 | tee -a "{log}"
        """

rule make_peak_stats:
    input:
        logs = lambda wildcards: define_logs_final_input(wildcards),
        stats_pseudoreps = lambda wildcards: f"results/{wildcards.env}/reports/stats_pseudoreps__{wildcards.sample_name}.txt"
        ## maybe a better solution is to append a stat file with wc -l as they are generated, or to create a new stat file for each file, then accessible in bash by regex on the samplename
    output:
        stat_file = "results/{env}/reports/summary_{env}_peak_stats_{sample_name}.txt",
        log = "results/{env}/logs/called_peaks__{sample_name}.log"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        paired = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'paired'),
        peaktype = lambda wildcards: get_peaktype_for_env(_resolve_assay(wildcards.sample_name), wildcards.env),
        line = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'line'),
        tissue = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'tissue'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'sample_type'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome'),
        env = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'env'),
        tf_name = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'extra_info'),
        rep1 = lambda wildcards: get_replicate_name(wildcards, 0),
        rep2 = lambda wildcards: get_replicate_name(wildcards, 1)
    threads: config["resources"]["make_peak_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["make_peak_stats"]["mem_mb"],
        tmp_mb=config["resources"]["make_peak_stats"]["tmp_mb"],
        qos=config["resources"]["make_peak_stats"]["qos"]
    shell:
        """
        nrep1=$(awk '{{print $1,$2,$3}}' {params.rep1} | sort -k1,1 -k2,2n -u | wc -l)
        if [[ "{params.rep2}" == "missingrep" ]]; then
            nrep2=0
        else
            nrep2=$(awk '{{print $1,$2,$3}}' {params.rep2} | sort -k1,1 -k2,2n -u | wc -l)
        fi
        merged=$(grep "Merged" {input.stats_pseudoreps} | cut -d"=" -f2)
        pseudos=$(grep "Pseudos" {input.stats_pseudoreps} | cut -d"=" -f2)
        idr=$(grep "IDR" {input.stats_pseudoreps} | cut -d"=" -f2)
        selected=$(grep "Selected" {input.stats_pseudoreps} | cut -d"=" -f2)
        printf "Line\tTissue\tMark\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tCommon_peaks\tPeaks_in_merged\tPeaks_in_pseudo_reps\tPeaks_in_idr\tSelected_peaks\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.line} -v t={params.tissue} -v m={params.sample_type} -v r={params.ref_genome} -v a=${{nrep1}} -v b=${{nrep2}} -v c=${{merged}} -v d=${{pseudos}} -v e=${{idr}} -v f=${{selected}} 'BEGIN {{if (c==0) {{x=a}} else {{x=c}}; print l,t,m,r,a,b,c,d,e,f" ("f/x*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        """

rule find_motifs_in_file:
    input:
        define_chipseq_target_file
    output:
        temp_bed = temp("results/{env}/motifs/temp_regions_{peak_file}.bed"),
        temp_fa = temp("results/{env}/motifs/temp_regions_{peak_file}.fa"),
        touch = "results/{env}/chkpts/motifs__{peak_file}.done"
    wildcard_constraints:
        env = "ChIP"
    params:
        env = lambda wildcards: wildcards.env,
        peak_file = lambda wildcards: wildcards.peak_file,
        jaspar_db = config['jaspar_db']
    log:
        temp(return_log_chip("{env}","{peak_file}", "motifs", ""))
    conda: CONDA_ENV_IDR
    threads: config["resources"]["find_motifs_in_file"]["threads"]
    resources:
        mem_mb=config["resources"]["find_motifs_in_file"]["mem_mb"],
        tmp_mb=config["resources"]["find_motifs_in_file"]["tmp_mb"],
        qos=config["resources"]["find_motifs_in_file"]["qos"]
    shell:
        """
        {{
        peakfile="{input[0]}"
        ext=${{peakfile##*.}}
        if [[ "${{ext}}" == "narrowPeak" ]]; then
            printf "\nGetting peak fasta sequences around the summit for narrowPeak file: ${{peakfile}}\n"
            sort -k1,1 -k2,2n -k5,5nr ${{peakfile}} | awk -v OFS="\t" '{{a=$1":"$2":"$3; if (a!=n) {{s=$2+$10; print $1,s-200,s+200,$4;}} n=$1":"$2":"$3}}' > {output.temp_bed}
        elif [[ "${{ext}}" == "broadPeak" ]]; then
            printf "\nGetting peak fasta sequences around the middle for broadPeak file: ${{peakfile}}\n"
            sort -k1,1 -k2,2n -k5,5nr ${{peakfile}} | awk -v OFS="\t" '{{s=int(($2+$3)/2); t=($3-$2); if (t<500) print $1,$2,$3,$4; else print $1,s-200,s+200,$4}}' > {output.temp_bed}
        elif [[ "${{ext}}" == "bedPeak" ]]; then
            printf "\nGetting peak fasta sequences for bed file: ${{peakfile}}\n"
            cat ${{peakfile}} | awk -v OFS="\t" '{{s=$2+$10; print $1,s-200,s+200,$4}}' > {output.temp_bed}
        else
            printf "\nGetting peak fasta sequences for unknown file format: ${{peakfile}}\n"
            cat ${{peakfile}} | awk -v OFS="\t" '{{s=int(($2+$3)/2); t=($3-$2); if (t<500) print $1,$2,$3,$4; else print $1,s-200,s+200,$4}}' > {output.temp_bed}
        fi
        head {output.temp_bed}
        bedtools getfasta -name -fi {input[1]} -bed {output.temp_bed} > {output.temp_fa}
        head {output.temp_fa}
        printf "\nGetting motifs for {params.peak_file} with meme version:\n"
        meme -version
        meme-chip -oc results/{params.env}/motifs/{params.peak_file}/meme -meme-p {threads} -meme-nmotifs 10 -streme-nmotifs 10 {output.temp_fa}
        if [[ -s results/{params.env}/{params.peak_file}/meme/combined.meme ]]; then
            printf "\nLooking for similar motifs in JASPAR database with tomtom\n"
            tomtom -oc results/{params.env}/motifs/{params.peak_file}/tomtom/ results/{params.env}/motifs/{params.peak_file}/meme/combined.meme {params.jaspar_db}
        fi
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule perform_pairwise_diff_peaks:
    input:
        peak_file1 = lambda wildcards: define_input_manorm(wildcards, "peaks1"),
        peak_file2 = lambda wildcards: define_input_manorm(wildcards, "peaks2"),
        read_file1 = lambda wildcards: define_input_manorm(wildcards, "reads1"),
        read_file2 = lambda wildcards: define_input_manorm(wildcards, "reads2")
    output:
        result = "results/{env}/peaks/{sample1}_vs_{sample2}/{sample1}_vs_{sample2}_all_MAvalues.xls"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        diffpeaks = lambda wildcards: define_input_manorm(wildcards, "params"),
        peakformat = lambda wildcards: define_input_manorm(wildcards, "format"),
        output_folder = lambda wildcards: f"results/{wildcards.env}/peaks/{wildcards.sample1}_vs_{wildcards.sample2}"
    log:
        temp(return_log_chip("{env}","{sample1}_vs_{sample2}", "diff_peaks", ""))
    conda: CONDA_ENV_CHIP
    threads: config["resources"]["perform_pairwise_diff_peaks"]["threads"]
    resources:
        mem_mb=config["resources"]["perform_pairwise_diff_peaks"]["mem_mb"],
        tmp_mb=config["resources"]["perform_pairwise_diff_peaks"]["tmp_mb"],
        qos=config["resources"]["perform_pairwise_diff_peaks"]["qos"]
    shell:
        """
        {{
        if [[ "{params.peakformat}" == "narrow" ]]; then
            printf "\nComparing {wildcards.sample1} with {wildcards.sample2} with .narrowPeak files with MAnorm version:\n"
            manorm --version
            manorm --p1 {input.peak_file1} --p2 {input.peak_file2} --r1 {input.read_file1} --r2 {input.read_file2} --n1 {wildcards.sample1} --n2 {wildcards.sample2} -o {params.output_folder} --rf "bam" --pf "narrowpeak" {params.diffpeaks}
        elif [[ "{params.peakformat}" == "broad" ]]; then
            printf "\nComparing {wildcards.sample1} with {wildcards.sample2} with .broadPeak files with MAnorm version:\n"
            manorm --version
            manorm --p1 {input.peak_file1} --p2 {input.peak_file2} --r1 {input.read_file1} --r2 {input.read_file2} --n1 {wildcards.sample1} --n2 {wildcards.sample2} -o {params.output_folder} --rf "bam" --pf "broadpeak" {params.diffpeaks}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule all_chip:
    input:
        final = lambda wildcards: define_final_chip_output(wildcards.ref_genome)
    output:
        touch = "results/{env}/chkpts/{env}_analysis__{analysis_name}__{ref_genome}.done"
    wildcard_constraints:
        env = "ChIP"
    localrule: True
    shell:
        """
        touch {output.touch}
        """
