# function to access logs more easily
def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER,"results","combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_combined_target_file(wildcards):
    tname = config['combined_target_file_label']
    target_name = wildcards.target_name
    ref_genome = wildcards.ref_genome
    
    if target_name == tname:
        return config['combined_target_file']
    elif target_name.startswith("combined_peaks"):
        file = f"results/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("all_genes") or target_name.startswith("protein_coding_genes"):
        file = f"results/combined/tracks/{ref_genome}__{target_name}.bed"
    else:
        raise ValueError(   
            f"{target_name} does not match possible files." 
            "It can be 'combined_peaks', or the value of "
            "'combined_target_file_name' in the config file"
        )
    
    return file

def get_heatmap_param(matrix, key):
    override = config.get(key)
    if override is not None:
        return override

    return config['heatmaps'][matrix][key]

def get_matrix_inputs(wildcards):
    stranded_heatmaps = config['stranded_heatmaps']
    bedfile = define_combined_target_file(wildcards)
    prefix = f"results/combined/matrix/matrix_{wildcards.matrix_param}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}__{wildcards.target_name}"
    with checkpoints.is_stranded.get(bedfile=bedfile).output[0].open() as f:
        if f.read().strip() == "stranded" and stranded_heatmaps:
            return [ f"{prefix}__plus.gz", f"{prefix}__minus.gz" ]
        else:
            return [ f"{prefix}__unstranded.gz" ]

def define_sort_options(wildcards):
    sort_options = config['heatmaps_sort_options']
    matrix_param = wildcards.matrix_param
    env = wildcards.env
    analysis_name=config['analysis_name']
    ref_genome = wildcards.ref_genome
    target_name = wildcards.target_name
    if sort_options == "no":
        return "--sortRegions keep"
    elif sort_options == "mean":
        return "--sortRegions descend --sortUsing mean"
    elif sort_options == "median":
        return "--sortRegions descend --sortUsing median"
    else:
        print("unclear option: no sorting done")
        return "--sortRegions keep"

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
        file = f"results/{row.env}/peaks/selected_peaks__{spname}.bedPeak"
        if row.env == "TF":
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
        files.append(f"results/{row.env}/peaks/selected_peaks__{spname}.bedPeak")
    
    return files

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
    result = ":".join(types)     
    return result

def define_key_for_heatmaps(wildcards, string):
    bigwigs = []
    labels = []
    marks = []
    unique_tf = set()
    unique_chip = set()
    unique_rna = set()
    unique_srna = set()
    unique_mc = set()
    grouped_bw = defaultdict(list)
    grouped_labs = defaultdict(list)
    srna_sizes = config['srna_heatmap_sizes']
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    strand = wildcards.strand
    
    if globenv == "all":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] != "mC") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        prefix = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}"
        reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
        if row.env == "ChIP":
            merged = f"FC__merged__{prefix}__merged__{row.ref_genome}.bw"
            onerep = f"FC__final__{prefix}__{reps[0]}__{row.ref_genome}.bw"
            bw = f"results/{row.env}/tracks/{merged}" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}"
            label = f"{row.line}_{row.tissue}_{row.sample_type}"
            grouped_bw[f"chip_{row.sample_type}"].append(bw)
            grouped_labs[f"chip_{row.sample_type}"].append(label)
            unique_chip.add(row.sample_type)
        elif row.env == "TF":
            merged = f"FC__merged__{prefix}__merged__{row.ref_genome}.bw"
            onerep = f"FC__final__{prefix}__{reps[0]}__{row.ref_genome}.bw"
            bw = f"results/{row.env}/tracks/{merged}" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}"
            label = f"{row.line}_{row.tissue}_{row.extra_info}"
            grouped_bw[f"tf_{row.extra_info}"].append(bw)
            grouped_labs[f"tf_{row.extra_info}"].append(label)
            unique_tf.add(row.extra_info)
        elif row.env == "RNA":
            if strand == "unstranded":
                merged = f"{prefix}__merged__{row.ref_genome}"
                onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                bw1 = f"results/{row.env}/tracks/{merged}__plus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__plus.bw"
                bw2 = f"results/{row.env}/tracks/{merged}__minus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__minus.bw"
                label = f"{row.line}_{row.tissue}_{row.sample_type}"
                grouped_bw[f"{row.data_type}_plus"].append(bw1)
                grouped_bw[f"{row.data_type}_minus"].append(bw2)
                grouped_labs[f"{row.data_type}_plus"].append(f"{label}_plus")
                grouped_labs[f"{row.data_type}_minus"].append(f"{label}_minus")
                unique_rna.add(row.data_type)
            else:
                merged = f"{prefix}__merged__{row.ref_genome}"
                onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                bw = f"results/{row.env}/tracks/{merged}__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{strand}.bw"
                label = f"{row.line}_{row.tissue}_{row.sample_type}"
                grouped_bw[f"{row.data_type}_stranded"].append(bw)
                grouped_labs[f"{row.data_type}_stranded"].append(f"{label}")
                unique_rna.add(row.data_type)
        elif row.env == "sRNA":
            for size in srna_sizes:
                if strand == "unstranded":
                    merged = f"{prefix}__merged__{row.ref_genome}"
                    onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                    bw1 = f"results/{row.env}/tracks/{merged}__{size}nt__plus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__plus.bw"
                    bw2 = f"results/{row.env}/tracks/{merged}__{size}nt__minus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__minus.bw"
                    label = f"{row.line}_{row.tissue}_sRNA_{size}nt"
                    grouped_bw[f"sRNA_{size}_plus"].append(bw1)
                    grouped_bw[f"sRNA_{size}_minus"].append(bw2)
                    grouped_labs[f"sRNA_{size}_plus"].append(f"{label}_plus")
                    grouped_labs[f"sRNA_{size}_minus"].append(f"{label}_minus")
                    unique_srna.add(f"sRNA_{size}")
                else:
                    merged = f"{prefix}__merged__{row.ref_genome}"
                    onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                    bw = f"results/{row.env}/tracks/{merged}__{size}nt__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__{strand}.bw"
                    label = f"{row.line}_{row.tissue}_sRNA_{size}nt"
                    grouped_bw[f"sRNA_{size}_stranded"].append(bw)
                    grouped_labs[f"sRNA_{size}_stranded"].append(f"{label}")
                    unique_srna.add(f"sRNA_{size}")
        elif row.env == "mC":
            merged = f"{prefix}__merged__{row.ref_genome}"
            onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
            for context in ["CG","CHG","CHH"]:
                bw = f"results/{row.env}/tracks/{merged}__{context}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{context}.bw"
                label = f"{row.line}_{row.tissue}_m{context}"
                grouped_bw[f"m{context}"].append(bw)
                grouped_labs[f"m{context}"].append(f"{label}")
                unique_mc.add(f"m{context}")
                    
    bigwigs = (
        sum([grouped_bw.get(f"chip_{chip}", []) for chip in sorted(unique_chip)], []) + 
        sum([grouped_bw.get(f"tf_{tf}", []) for tf in sorted(unique_tf)], []) + 
        sum([grouped_bw.get(f"{rna}_plus", []) + grouped_bw.get(f"{rna}_minus", []) + grouped_bw.get(f"{rna}_stranded", []) for rna in sorted(unique_rna)], []) + 
        sum([grouped_bw.get(f"{srna}_plus", []) + grouped_bw.get(f"{srna}_minus", []) + grouped_bw.get(f"{arna}_stranded", []) for srna in sorted(unique_srna)], []) +
        sum([grouped_bw.get(f"{mc}", []) for mc in sorted(unique_mc), [])
    )
    labels = (
        sum([grouped_labs.get(f"chip_{chip}", []) for chip in sorted(unique_chip)], []) + 
        sum([grouped_labs.get(f"tf_{tf}", []) for tf in sorted(unique_tf)], []) + 
        sum([grouped_labs.get(f"{rna}_plus", []) + grouped_labs.get(f"{rna}_minus", []) + grouped_labs.get(f"{rna}_stranded", []) for rna in sorted(unique_rna)], []) + 
        sum([grouped_labs.get(f"{srna}_plus", []) + grouped_labs.get(f"{srna}_minus", []) + grouped_labs.get(f"{arna}_stranded", []) for srna in sorted(unique_srna)], []) +
        sum([grouped_labs.get(f"{mc}", []) for mc in sorted(unique_mc), [])
    )
    marks = ( sorted(unique_chip) + sorted(unique_tf) + [f"{rna}_plus" for rna in sorted(unique_rna)] + [f"{rna}_minus" for rna in sorted(unique_rna)] + [f"{srna}_plus" for srna in sorted(unique_srna)] + [f"{srna}_minus" for srna in sorted(unique_srna)] + sorted(unique_mc) ) if strand == "unstranded" else ( sorted(unique_chip) + sorted(unique_tf) + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )
    
    if string == "bigwigs":
        return bigwigs
    elif string == "labels":
        return labels
    elif string == "marks":
        return marks

def define_final_combined_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    analysis_name = config['analysis_name']
    text_files = []
    plot_files = []
    
    all_analysis_samples = analysis_samples[ analysis_samples['ref_genome'] == ref_genome ].copy()
    chip_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    tf_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'TF') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    mc_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'mC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    rna_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'RNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    srna_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'sRNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    
    if len(chip_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__ChIP__{analysis_name}__{ref_genome}.pdf")
    
    if len(tf_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__TF__{analysis_name}__{ref_genome}.pdf")
    
    if len(chip_analysis_samples) >=1 and len(tf_analysis_samples) >=1:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__all_chip__{analysis_name}__{ref_genome}.pdf")
    
    if len(mc_analysis_samples) >=1:
        if len(all_analysis_samples) > len(mc_analysis_samples):
            plot_files.append(f"results/combined/plots/Heatmap_sorted__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap_sorted__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap_sorted__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        else:
            plot_files.append(f"results/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        
        plot_files.append(f"results/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
    else:
        plot_files.append(f"results/combined/plots/Heatmap__regions__all__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Heatmap__tss__all__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Heatmap__tes__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    
    plot_files.append(f"results/combined/plots/Profile__regions__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    plot_files.append(f"results/combined/plots/Profile__tss__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    plot_files.append(f"results/combined/plots/Profile__tes__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    
    if analysis:
        results = plot_files + text_files
    else:
        results = []
    
    return results

###
# rules to look for header or strandedness of bedfile
rule has_header:
    input:
        bedfile = "{bedfile}"
    output:
        file = temp("{bedfile}.header")
    localrule: True
    run:
        with open(input.bedfile) as f:
            first_line = f.readline().strip().split('\t')
            try:
                res = "no" if (int(first_line[1]) >=0 and int(first_line[2]) >=0) else "yes"
            except (ValueError, IndexError):
                res = "yes"
        
        with open(output.file, "w") as out:
            out.write(res + "\n")
              
checkpoint is_stranded:
    input:
        bedfile = "{bedfile}",
        header = "{bedfile}.header"
    output:
        file = temp("{bedfile}.stranded")
    localrule: True
    run:
        with open(input.header) as h:
            header = h.read().strip()
        
        strand_values = set()
        with open(input.bedfile) as f:
            if header == "yes":
                next(f)
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    return False
                else:
                    strand_values.add(cols[5])
                    
        with open(output.file, "w") as out:
            if strand_values.issubset({"+","-"}):
                out.write("stranded" + "\n")
            else:
                out.write("unstranded" + "\n")

            
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
        temp1_file = temp("results/combined/bedfiles/temp1_combined_peaks__{env}__{analysis_name}__{ref_genome}.bed"),
        temp2_file = temp("results/combined/bedfiles/temp2_combined_peaks__{env}__{analysis_name}__{ref_genome}.bed"),
        merged_file = "results/combined/bedfiles/combined_peaks__{env}__{analysis_name}__{ref_genome}.bed"
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
        bedfile = lambda wildcards: define_combined_target_file(wildcards),
        region_file = lambda wildcards: f"results/combined/tracks/{wildcards.ref_genome}__all_genes.bed",
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        header = lambda wildcards: f"{define_combined_target_file(wildcards)}.header"
    output:
        temp_bedfile = temp("results/combined/bedfiles/temp__{target_name}__{ref_genome}.bed"),
        annotated_file = "results/combined/bedfiles/annotated__{target_name}__{ref_genome}.bed"
    params:
        target_name = lambda wildcards: wildcards.target_name,
    log:
        temp(return_log_combined("{target_name}", "{ref_genome}", "annotate_bedfile"))
    conda: CONDA_ENV
    threads: config["resources"]["get_annotations_for_bedfile"]["threads"]
    resources:
        mem=config["resources"]["get_annotations_for_bedfile"]["mem"],
        tmp=config["resources"]["get_annotations_for_bedfile"]["tmp"]
    shell:
        """
        {{
        printf "Annotating {params.target_name} to the closest genes\n"
        header=$(cat {input.header})
        if [[ "${{header}}" == "no" ]]; then
            awk -v OFS="\t" -v n={params.target_name} '{{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile}
        else
            awk -v OFS="\t" -v n={params.target_name} 'NR>1 {{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile}
        fi
        printf "Chr\tStart\tStop\tPeakID\tDistance\tGene_strand\tGID\tCategory\n" > {output.annotated_file}
        bedtools closest -a {output.temp_bedfile} -b {input.region_file} -g {input.chrom_sizes} -D ref | awk -v OFS="\t" '{{if ($10=="+") print $1,$2,$3,$4,$11,$10,$8; else print $1,$2,$3,$4,-$11,$10,$8}}' | awk -F"[=;]" -v OFS="\t" '{{print $1,$2}}' | sed 's/gene://' | awk -v OFS="\t" '{{if ($5<-2000) {{d="Distal_downstream"}} else if ($5<0) {{d="Terminator"}} else if ($5==0) {{d="Gene_body"}} else if ($5>2000) {{d="Distal_upstream"}} else {{d="Promoter"}}; print $1,$2,$3,$4,$5,$6,$8,d}}' >> {output.annotated_file}
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_upset_peaks:
    input:
        mergedfile = "results/combined/bedfiles/{target_name}__{env}__{analysis_name}__{ref_genome}.bed",
        annotatedfile = "results/combined/bedfiles/annotated__{target_name}__{env}__{analysis_name}__{ref_genome}.bed"
    output:
        plot = "results/combined/plots/Upset_{target_name}__{env}__{analysis_name}__{ref_genome}.pdf"
    params:
        env = lambda wildcards: wildcards.env,
        types = lambda wildcards: define_sample_types_for_upset(wildcards),
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot.R")
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_upset_{target_name}_{env}"))
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
# rule to plot heatmaps
rule making_stranded_matrix_on_targetfile:
    input:
        bigwigs = lambda wildcards: define_key_for_heatmaps(wildcards, "bigwigs"),
        target_file = lambda wildcards: define_combined_target_file(wildcards),
        header = lambda wildcards: f"{define_combined_target_file(wildcards)}.header"
    output:
        temp = temp("results/combined/matrix/temp_file_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}_{strand}.bed"),
        matrix = temp("results/combined/matrix/matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}__{strand}.gz")
    wildcard_constraints:
        strand = "plus|minus|unstranded"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        labels = lambda wildcards: define_key_for_heatmaps(wildcards, "labels"),
        marks = lambda wildcards: define_key_for_heatmaps(wildcards, "marks"),
        matrix = lambda wildcards: wildcards.matrix_param,
        strand = lambda wildcards: wildcards.strand,
        params = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'base'),
        bs = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'bs'),
        before = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'before'),
        after = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'after'),
        middle = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'middle')
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "making_matrix_{matrix_param}_{target_name}_{strand}"))
    conda: CONDA_ENV
    threads: config["resources"]["making_stranded_matrix_on_targetfile"]["threads"]
    resources:
        mem=config["resources"]["making_stranded_matrix_on_targetfile"]["mem"],
        tmp=config["resources"]["making_stranded_matrix_on_targetfile"]["tmp"]
    shell:
        """
        {{
        header="$(cat {input.header})"
        if [[ "{params.strand}" == "unstranded" ]]; then
            if [[ "${{header}}" == "no" ]]; then
                cat {input.target_file} > {output.temp}
            else
                awk 'NR>1' {input.target_file} > {output.temp}
            fi
        else
            case "{params.strand}" in
                plus)   sign="+";;
                minus)  sign="-";;
            esac
            awk -v s=${{sign}} '$6==s' {input.target_file} > {output.temp}
        fi
        echo "{params.labels}" | xargs -n1 > "results/combined/matrix/labels_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt"
        echo "{params.marks}" | xargs -n1 > "results/combined/matrix/marks_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt"
        printf "Making {params.strand} strand {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
        computeMatrix {params.params} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.labels} -bs {params.bs} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        }} 2>&1 | tee -a "{log}"
        """
                
rule merging_matrix:
    input:
        get_matrix_inputs
    output:
        "results/combined/matrix/final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.gz"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "merging_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["merging_matrix"]["threads"]
    resources:
        mem=config["resources"]["merging_matrix"]["mem"],
        tmp=config["resources"]["merging_matrix"]["tmp"]
    shell:
        """
        {{
        nfile=$(echo {input} | wc -w)
        if [[ ${{nfile}} -eq 2 ]]; then
            printf "\nMerging stranded matrices aligned by {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations rbind -m {input} -o {output}
        else
            cp {input} {output}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule computing_matrix_scales:
    input:
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.gz",
        target_file = lambda wildcards: define_combined_target_file(wildcards),
        header = lambda wildcards: f"{define_combined_target_file(wildcards)}.header"
    output:
        params_heatmap = "results/combined/matrix/params_heatmap_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        params_profile = "results/combined/matrix/params_profile_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        params_regions = "results/combined/matrix/params_regions_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        temp_values = "results/combined/matrix/temp_values_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        temp_profile = temp("results/combined/matrix/temp_profile_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf"),
        temp_profile_values = "results/combined/matrix/temp_profile_values_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        scales = config['heatmaps_scales'],
        profile = config['profiles_scale']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "getting_scales_matrix_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["computing_matrix_scales"]["threads"]
    resources:
        mem=config["resources"]["computing_matrix_scales"]["mem"],
        tmp=config["resources"]["computing_matrix_scales"]["tmp"]
    shell:
        """
        {{        
        header="$(cat {input.header})"
        count=$(wc -l {input.target_file} | cut -d' ' -f 1)
        if [[ "${{header}}" == "yes" ]]; then
            count=$((count-1))
        fi
        awk -v ORS="" -v r=${{count}} -v n={params.target_name} 'BEGIN {{print "--regionsLabel "n"("r")"}}' > {output.params_regions}

        if [[ "{params.scales}" == "default" ]]; then
            touch {output.params_heatmap}
            touch {output.params_profile}
            touch {output.temp_values}
            touch {output.temp_profile}
            touch {output.temp_profile_values}
            
        elif [[ "{params.scales}" == "type" ]]; then
            printf "Getting scales per type for {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations dataRange -m {input.matrix} > {output.temp_values}
            plotProfile -m {input.matrix} -out {output.temp_profile} --averageType {params.profile} --outFileNameData {output.temp_profile_values}
            
            mins=()
            maxs=()
            ymins=()
            ymaxs=()
            while read mark
            do
                zmini=$(grep "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=999999}} {{a=$5; if (a<m) m=a;}} END {{print m}}')
                zmaxi=$(grep "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=-999999}} {{a=$6; if (a>m) m=a;}} END {{print m}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" ]]; then
                    zmini="0"
                    zmaxi="0.005"
                fi
                
                ymini=$(grep "${{mark}}" {output.temp_profile_values} | awk '{{m=$3; for (i=3;i<=NF;i++) if ($i<m) m=$i; print m}}' | awk 'BEGIN {{m=99999}} {{if ($1<m) m=$1}} END {{if (m<0) a=m*1.2; else a=m*0.8; print a}}')
                ymaxi=$(grep "${{mark}}" {output.temp_profile_values} | awk '{{m=$3; for (i=3;i<=NF;i++) if ($i>m) m=$i; print m}}' | awk 'BEGIN {{m=-99999}} {{if ($1>m) m=$1}} END {{if (m<0) a=m*0.8; else a=m*1.2; print a}}')
                test=$(awk -v a=${{ymini}} -v b=${{ymaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ ${{test}} == "yes" ]]; then
                    ymini=("0")
                    ymaxi=("0.01")
                fi
                num=$(grep "${{mark}}" {output.temp_values} | wc -l)
                for i in $(seq 1 ${{num}})
                do
                    zmins+=("$zmini")
                    zmaxs+=("$zmaxi")
                    ymins+=("$ymini")
                    ymaxs+=("$ymaxi")
                done
            done < results/combined/matrix/marks_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt
            
            awk -v ORS="" -v a="${{zmins[*]}}" -v b="${{zmaxs[*]}}" 'BEGIN {{print "--zMin "a" --zMax "b}}' > {output.params_heatmap}
            awk -v ORS="" -v c="${{ymins[*]}}" -v d="${{ymaxs[*]}}" 'BEGIN {{print "--yMin "c" --yMax "d}}' > {output.params_profile}
        
        elif [[ "{params.scales}" == "sample" ]]; then
            printf "Getting scales per sample for {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations dataRange -m {input.matrix} > {output.temp_values}
            plotProfile -m {input.matrix} -out {output.temp_profile} --averageType {params.profile} --outFileNameData {output.temp_profile_values}
            
            zmins=()
            zmaxs=()
            ymins=()
            ymaxs=()
            while read sample
            do
                zmini=$(grep "${{sample}}" {output.temp_values} | awk '{{print $5}}')
                zmaxi=$(grep "${{sample}}" {output.temp_values} | awk '{{print $6}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" ]]; then
                    zmins+=("0")
                    zmaxs+=("0.005")
                else
                    zmins+=("$zmini")
                    zmaxs+=("$zmaxi")
                fi
                
                ymini=$(grep "${{sample}}" {output.temp_profile_values} | awk '{{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}}' | awk 'BEGIN {{m=99999}} {{if ($1<m) m=$1}} END {{if (m<0) a=m*1.2; else a=m*0.8; print a}}')
                ymaxi=$(grep "${{sample}}" {output.temp_profile_values} | awk '{{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}}' | awk 'BEGIN {{m=-99999}} {{if ($1>m) m=$1}} END {{if (m<0) a=m*0.8; else a=m*1.2; print a}}')
                test=$(awk -v a=${{ymini}} -v b=${{ymaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" ]]; then
                    ymins+=("0")
                    ymaxs+=("0.01")
                else
                    ymins+=("$ymini")
                    ymaxs+=("$ymaxi")
                fi
            done < results/combined/matrix/labels_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt
            
            awk -v ORS="" -v a="${{zmins[*]}}" -v b="${{zmaxs[*]}}" 'BEGIN {{print "--zMin "a" --zMax "b}}' > {output.params_heatmap}
            awk -v ORS="" -v c="${{ymins[*]}}" -v d="${{ymaxs[*]}}" 'BEGIN {{print "--yMin "c" --yMax "d}}' > {output.params_profile}
        else
            printf "{params.scales} unknown. Returning default\n"
            touch {output.params_heatmap}
            touch {output.params_profile}
            touch {output.temp_values}
            touch {output.temp_profile}
            touch {output.temp_profile_values}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_heatmap_on_targetfile:
    input:
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.gz",
        params_regions = "results/combined/matrix/params_regions_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        params_heatmap = "results/combined/matrix/params_heatmap_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        params_profile = "results/combined/matrix/params_profile_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"
    output:
        plot = "results/combined/plots/Heatmap__{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf",
        sorted_regions = "results/combined/matrix/Heatmap__{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}_sorted_regions.bed"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        env = lambda wildcards: wildcards.env,
        plot_params = lambda wildcards: config['heatmaps_plot_params'][wildcards.env],
        sort = lambda wildcards: define_sort_options(wildcards)
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "plot_heatmap_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_heatmap_on_targetfile"]["threads"]
    resources:
        mem=config["resources"]["plotting_heatmap_on_targetfile"]["mem"],
        tmp=config["resources"]["plotting_heatmap_on_targetfile"]["tmp"]
    shell:
        """
        new_params="$(cat {input.params_regions} {input.params_heatmap} {input.params_profile})"
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        printf "Plotting heatmap {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        plotHeatmap -m {input.matrix} -out {output.plot} {params.plot_params} {params.sort} ${{new_params}} ${{add}} --outFileSortedRegions {output.sorted_regions}
        """

rule sort_heatmap:
    input: 
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.gz",
        sorted_regions = "results/combined/matrix/Heatmap__{matrix_param}__all__{analysis_name}__{ref_genome}__{target_name}_sorted_regions.bed",
        params_regions = "results/combined/matrix/params_regions_final_matrix_{matrix_param}__all__{analysis_name}__{ref_genome}__{target_name}.txt"
    output:
        temp_matrix = temp("results/combined/matrix/temp_sorted_final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.gz"),
        matrix = "results/combined/matrix/sorted_final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.gz"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param
    log:
        temp(return_log_combined("{analysis_name}", "mC_{ref_genome}", "sort_heatmap_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["sort_heatmap"]["threads"]
    resources:
        mem=config["resources"]["sort_heatmap"]["mem"],
        tmp=config["resources"]["sort_heatmap"]["tmp"]
    shell:
        """
        printf "Sorting heatmap {params.matrix} for mC {params.target_name} on {params.ref_genome}\n"
        label="$(cat {input.params_regions} | cut -d" " -f 2)"
        computeMatrixOperations relabel -m {input.matrix} --groupLabels ${{label}} -o {output.temp_matrix}
        computeMatrixOperations sort -m {output.temp_matrix} -R {input.sorted_regions} -o {output.matrix}
        """

rule plotting_sorted_heatmap_on_targetfile:
    input:
        matrix = "results/combined/matrix/sorted_final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.gz",
        params_heatmap = "results/combined/matrix/params_heatmap_final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.txt",
        params_profile = "results/combined/matrix/params_profile_final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.txt"
    output:
        plot = "results/combined/plots/Heatmap_sorted__{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.pdf"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        plot_params = lambda wildcards: config['heatmaps_plot_params']['mC']
    log:
        temp(return_log_combined("{analysis_name}", "mC_{ref_genome}", "plot_sorted_heatmap_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_sorted_heatmap_on_targetfile"]["threads"]
    resources:
        mem=config["resources"]["plotting_sorted_heatmap_on_targetfile"]["mem"],
        tmp=config["resources"]["plotting_sorted_heatmap_on_targetfile"]["tmp"]
    shell:
        """
        new_params="$(cat {input.params_heatmap} {input.params_profile})"
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        printf "Plotting heatmap {params.matrix} for mC {params.target_name} on {params.ref_genome}\n"
        plotHeatmap -m {input.matrix} -out {output.plot} {params.plot_params} --sortRegions 'keep' ${{new_params}} ${{add}}
        """

rule plotting_profile_on_targetfile:
    input:
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.gz",
        params_regions = "results/combined/matrix/params_regions_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        params_profile = "results/combined/matrix/params_profile_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"
    output:
        plot1 = "results/combined/plots/Profile__{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf",
        plot2 = "results/combined/plots/Profile_pergroup__{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        env = lambda wildcards: wildcards.env,
        plot_params = config['profiles_plot_params']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "plot_profile_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_profile_on_targetfile"]["threads"]
    resources:
        mem=config["resources"]["plotting_profile_on_targetfile"]["mem"],
        tmp=config["resources"]["plotting_profile_on_targetfile"]["tmp"]
    shell:
        """
        {{
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        printf "Plotting profile {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        new_params="$(cat {input.params_regions} {input.params_profile})"
        printf "${{new_params}}"
        plotProfile -m {input.matrix} -out {output.plot1} {params.plot_params} ${{new_params}} ${{add}}
        
        printf "Plotting per group profile {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        ymin=$(cat {input.params_profile} | awk 'BEGIN {{y=99999}} {{for (i=1; i<=NF; i++) {{if ($i == "--yMin") {{for (j=i+1; j<=NF && $j !~ /^--/; j++) {{if ($j<y) y=$j}} break}} }} }} END {{print y}}' )
        ymax=$(cat {input.params_profile} | awk 'BEGIN {{y=-99999}} {{for (i=1; i<=NF; i++) {{if ($i == "--yMax") {{for (j=i+1; j<=NF && $j !~ /^--/; j++) {{if ($j>y) y=$j}} break}} }} }} END {{print y}}' )
        new_params="$(cat {input.params_regions})"
        plotProfile -m {input.matrix} -out {output.plot2} {params.plot_params} ${{new_params}} --yMin ${{ymin}} --yMax ${{ymax}} ${{add}} --perGroup
        }} 2>&1 | tee -a "{log}"
        """

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
