CONDA_ENV_UPSET=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_upset.yaml")

def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER,"results","combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_combined_target_file(wildcards):
    heatname = config['heatmap_target_file_label']
    browsername = config['browser_target_file_label']
    target_name = wildcards.target_name
    ref_genome = wildcards.ref_genome
    
    if target_name == heatname:
        return config['heatmap_target_file']
    elif target_name == browsername:
        return config['browser_target_file']
    elif target_name == "full_chromosomes":
        return f"results/combined/bedfiles/full_chromosomes__{ref_genome}.bed"
    elif target_name.startswith("combined_peaks"):
        file = f"results/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("combined_clusters"):
        file = f"results/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("combined_TSS"):
        file = f"results/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("all_genes") or target_name.startswith("protein_coding_genes"):
        file = f"results/combined/bedfiles/{ref_genome}__{target_name}.bed"
    elif target_name.startswith("all_TEs"):
        file = f"genomes/{ref_genome}/{ref_genome}__TE_file.bed"
    else:
        raise ValueError(   
            f"{target_name} does not match possible files. It can be 'combined_peaks', 'combined_clusters', 'all_genes', 'all_TEs'" 
            "'full_chromosomes' or the value of 'heatmap_target_file_label' or 'browser_target_file_label' in the config file"
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

def define_samples_for_upset(wildcards, string):
    names = []
    files = []
    types = set()
    ref_genome = wildcards.ref_genome
    srna_sizes = config['srna_heatmap_sizes']
    globenv = wildcards.env
    allreps = config['upset_allreps']
    if globenv == "all_chip":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF","ATAC"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    elif globenv == "RAMPAGE":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['sample_type'] == "RAMPAGE") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = sample_name_str(row, 'analysis')
        if row.env == "TF":
            if allreps:
                for replicate in analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), []):
                    peaktype = get_peaktype_for_env(row.sample_type, row.env)
                    paired = get_sample_info_from_name(spname, analysis_samples, 'paired')
                    if paired == "PE":
                        prefix = "peaks_pe"
                    else:
                        prefix = "peaks_se"
                    file = f"results/{row.env}/peaks/{prefix}__final__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{replicate}__{row.ref_genome}_peaks.{peaktype}Peak"
                    label = f"{row.line}_{row.tissue}_{row.extra_info}_{replicate}"
                    names.append(f"{label}:{file}")
                    files.append(file)
                    types.add(row.extra_info)
            else:
                file = f"results/{row.env}/peaks/selected_peaks__{spname}.bedPeak"
                label = f"{row.line}_{row.tissue}_{row.extra_info}"
                names.append(f"{label}:{file}")
                files.append(file)
                types.add(row.extra_info)
        elif row.env in ["ChIP", "ATAC"]:
            if allreps:
                for replicate in analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), []):
                    peaktype = get_peaktype_for_env(row.sample_type, row.env)
                    paired = get_sample_info_from_name(spname, analysis_samples, 'paired')
                    if row.env == "ATAC":
                        prefix = "peaks_atac"
                    elif paired == "PE":
                        prefix = "peaks_pe"
                    else:
                        prefix = "peaks_se"
                    file = f"results/{row.env}/peaks/{prefix}__final__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{replicate}__{row.ref_genome}_peaks.{peaktype}Peak"
                    label = f"{row.line}_{row.tissue}_{row.sample_type}_{replicate}"
                    names.append(f"{label}:{file}")
                    files.append(file)
                    types.add(row.sample_type)
            else:
                file = f"results/{row.env}/peaks/selected_peaks__{spname}.bedPeak"
                label = f"{row.line}_{row.tissue}_{row.sample_type}"
                names.append(f"{label}:{file}")
                files.append(file)
                types.add(row.sample_type)
        elif globenv == "RAMPAGE":
            if allreps:
                for replicate in analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), []):
                    file = f"results/RNA/TSS/TSS__final__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{replicate}__{row.ref_genome}_peaks.narrowPeak"
                    label = f"{row.line}_{row.tissue}_{replicate}"
                    names.append(f"{label}:{file}")
                    files.append(file)
                    types.add(f"{row.line}_{row.tissue}")
            else:
                file = f"results/RNA/TSS/TSS__merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}_peaks.narrowPeak"
                label = f"{row.line}_{row.tissue}"
                names.append(f"{label}:{file}")
                files.append(file)
                types.add(f"{row.line}_{row.tissue}")
        elif row.env == "sRNA":
            for replicate in analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), []):
                file = f"results/sRNA/mapped/{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{replicate}__{row.ref_genome}/clusters.bed"
                label = f"{row.line}_{row.tissue}_{replicate}"
                names.append(f"{label}:{file}")
                files.append(file)

    if globenv == "sRNA":
        srna_min = config['srna_min_size']
        srna_max = config['srna_max_size']
        types = [f"{s}nt" for s in range(srna_min, srna_max + 1)]
        types += ["MIRNA", "Others"]
        ordered = ":".join(types)
    else:
        ordered = ":".join(sorted(types))
    
    if string == "pairs":
        return names
    elif string == "files":
        return files
    elif string == "types":
        return ordered

# def define_bedfiles_per_env_and_ref(wildcards):
    # files = []
    # ref_genome = wildcards.ref_genome
    # globenv = wildcards.env
    # if globenv == "all_chip":
        # filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF","ATAC"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    # elif globenv == "RAMPAGE":
        # filtered_analysis_samples = analysis_samples[ (analysis_samples['sample_type'] == "RAMPAGE") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    # else:
        # filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    # for _, row in filtered_analysis_samples.iterrows():
        # spname = sample_name_str(row, 'analysis')
        # if globenv in ["all_chip", "ChIP", "TF", "ATAC"]:
            # files.append(f"results/{row.env}/peaks/selected_peaks__{spname}.bedPeak")
        # elif globenv == "sRNA":
            # files.extend(f"results/sRNA/mapped/{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{replicate}__{row.ref_genome}/clusters.bed"
                                      # for replicate in analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), []))
        # elif globenv == "RAMPAGE":
            # files.extend(f"results/RNA/TSS/TSS__final__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{replicate}__{row.ref_genome}_peaks.narrowPeak"
                                      # for replicate in analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), []))

    # return files

# def define_sample_types_for_upset(wildcards):
    # types = set()
    # ref_genome = wildcards.ref_genome
    # globenv = wildcards.env
    # if globenv == "all_chip":
        # filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","TF","ATAC"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    # elif globenv in ["ChIP", "TF", "ATAC"]:
        # filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    # elif globenv == "RAMPAGE":
        # filtered_analysis_samples = analysis_samples[ (analysis_samples['data_type'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    # else:
        # filtered_analysis_samples = None

    # if filtered_analysis_samples is not None:
        # for _, row in filtered_analysis_samples.iterrows():
            # if row.env in ["ChIP", "ATAC"]:
                # types.add(row.sample_type)
            # elif row.env == "TF":
                # types.add(row.extra_info)
            # elif row.sample_type == "RAMPAGE":
                # types.add(f"{row.line}__{row.tissue}")

    # if globenv == "sRNA":
        # srna_min = config['srna_min_size']
        # srna_max = config['srna_max_size']
        # types = [f"{s}nt" for s in range(srna_min, srna_max + 1)]
        # types += ["MIRNA", "Others"]
        # result = ":".join(types)
    # else:
        # result = ":".join(sorted(types))
    # return result

def define_upset_script(wildcards):
    globenv = wildcards.env
    if globenv in ["all_chip", "ChIP", "TF", "ATAC"]:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot_peaks.R")
    elif globenv == "sRNA":
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot_clusters.R")
    elif globenv == "RAMPAGE":
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot_TSS.R")

    return script

def assign_colors(keys, cmap_name="tab20"):
    cmap = plt.get_cmap(cmap_name)
    colors = {}
    for i,key in enumerate(sorted(set(keys))):
        colors[key] = mcolors.to_hex(cmap(i % cmap.N))
    return colors
    
def make_it_lighter(palette, factor):
    new_palette = {}
    for k, c in palette.items():
        rgb = mcolors.hex2color(c)
        h, l, s = colorsys.rgb_to_hls(*rgb)
        n = min(0.9, l*factor)
        new_palette[k] = mcolors.to_hex(colorsys.hls_to_rgb(h, n, s))
    return new_palette

def define_key_for_plots(wildcards, string):
    bigwigs = []
    labels = []
    marks = []
    unique_tf = set()
    unique_chip = set()
    unique_atac = set()
    unique_rna = set()
    unique_srna = set()
    unique_mc = set()
    grouped_bw = defaultdict(list)
    grouped_labs = defaultdict(list)
    label_to_mark = {}
    label_to_type = {}
    label_to_track = {}
    srna_sizes = config['srna_heatmap_sizes']
    plot_allreps = config['plot_allreps']
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    strand = getattr(wildcards, "strand", "unstranded")
    
    if globenv == "all":
        filtered_analysis_samples = analysis_samples[ analysis_samples['ref_genome'] == ref_genome ].copy()
    elif globenv == "most":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] != "mC") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        prefix = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}"
        reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
        if row.env == "ChIP":
            if not plot_allreps:
                merged = f"FC__merged__{prefix}__merged__{row.ref_genome}.bw"
                onerep = f"FC__final__{prefix}__{reps[0]}__{row.ref_genome}.bw"                
                bw = f"results/{row.env}/tracks/{merged}" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}"
                label = f"{row.line}_{row.tissue}_{row.sample_type}"
                grouped_bw[f"chip_{row.sample_type}"].append(bw)
                grouped_labs[f"chip_{row.sample_type}"].append(label)
                unique_chip.add(row.sample_type)
                label_to_mark[label] = row.sample_type
                label_to_type[label] = f"{row.line}_{row.tissue}"
                label_to_track[label] = row.sample_type
            else:
                for rep in reps:
                    bw = f"results/{row.env}/tracks/FC__final__{prefix}__{rep}__{row.ref_genome}.bw"
                    label = f"{row.line}_{row.tissue}_{row.sample_type}_{rep}"
                    grouped_bw[f"chip_{row.sample_type}"].append(bw)
                    grouped_labs[f"chip_{row.sample_type}"].append(label)
                    unique_chip.add(row.sample_type)
                    label_to_mark[label] = row.sample_type
                    label_to_type[label] = f"{row.line}_{row.tissue}"
                    label_to_track[label] = row.sample_type
            
        elif row.env == "TF":
            if not plot_allreps:
                merged = f"FC__merged__{prefix}__merged__{row.ref_genome}.bw"
                onerep = f"FC__final__{prefix}__{reps[0]}__{row.ref_genome}.bw"
                bw = f"results/{row.env}/tracks/{merged}" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}"
                label = f"{row.line}_{row.tissue}_{row.extra_info}"
                grouped_bw[f"tf_{row.extra_info}"].append(bw)
                grouped_labs[f"tf_{row.extra_info}"].append(label)
                unique_tf.add(row.extra_info)
                label_to_mark[label] = row.extra_info
                label_to_type[label] = f"{row.line}_{row.tissue}"
                label_to_track[label] = row.extra_info
            else:
                for rep in reps:
                    bw = f"results/{row.env}/tracks/FC__final__{prefix}__{rep}__{row.ref_genome}.bw"
                    label = f"{row.line}_{row.tissue}_{row.extra_info}_{rep}"
                    grouped_bw[f"tf_{row.extra_info}"].append(bw)
                    grouped_labs[f"tf_{row.extra_info}"].append(label)
                    unique_tf.add(row.extra_info)
                    label_to_mark[label] = row.extra_info
                    label_to_type[label] = f"{row.line}_{row.tissue}"
                    label_to_track[label] = row.extra_info
            
        elif row.env == "RNA":
            strandedness = config['rna_tracks'][row.data_type]['strandedness']
            if strandedness == "unstranded":
                if not plot_allreps:
                    merged = f"{prefix}__merged__{row.ref_genome}"
                    onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                    bw = f"results/{row.env}/tracks/{merged}__unstranded.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__unstranded.bw"
                    label = f"{row.line}_{row.tissue}_{row.sample_type}"
                    grouped_bw[f"{row.data_type}"].append(bw)
                    grouped_labs[f"{row.data_type}"].append(f"{label}")
                    unique_rna.add(row.data_type)
                    label_to_mark[label] = row.data_type
                    label_to_type[label] = f"{row.line}_{row.tissue}"
                    label_to_track[label] = row.data_type
                else:
                    for rep in reps:
                        bw = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__unstranded.bw"
                        label = f"{row.line}_{row.tissue}_{row.sample_type}_{rep}"
                        grouped_bw[f"{row.data_type}"].append(bw)
                        grouped_labs[f"{row.data_type}"].append(f"{label}")
                        unique_rna.add(row.data_type)
                        label_to_mark[label] = row.data_type
                        label_to_type[label] = f"{row.line}_{row.tissue}"
                        label_to_track[label] = row.data_type
            elif strand == "unstranded":
                if not plot_allreps:
                    merged = f"{prefix}__merged__{row.ref_genome}"
                    onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                    bw1 = f"results/{row.env}/tracks/{merged}__plus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__plus.bw"
                    bw2 = f"results/{row.env}/tracks/{merged}__minus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__minus.bw"
                    label = f"{row.line}_{row.tissue}_{row.sample_type}"
                    grouped_bw[f"{row.data_type}"].extend([bw1, bw2])
                    grouped_labs[f"{row.data_type}"].extend([f"{label}_plus", f"{label}_minus"])
                    unique_rna.add(row.data_type)
                    label_to_mark[f"{label}_plus"] = row.data_type
                    label_to_mark[f"{label}_minus"] = row.data_type
                    label_to_type[f"{label}_plus"] = f"{row.line}_{row.tissue}"
                    label_to_type[f"{label}_minus"] = f"{row.line}_{row.tissue}"
                    label_to_track[f"{label}_plus"] = f"{row.data_type}_plus"
                    label_to_track[f"{label}_minus"] = f"{row.data_type}_minus"
                else:
                    for rep in reps:
                        bw1 = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__plus.bw"
                        bw2 = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__minus.bw"
                        label = f"{row.line}_{row.tissue}_{row.sample_type}_{rep}"
                        grouped_bw[f"{row.data_type}"].extend([bw1, bw2])
                        grouped_labs[f"{row.data_type}"].extend([f"{label}_plus", f"{label}_minus"])
                        unique_rna.add(row.data_type)
                        label_to_mark[f"{label}_plus"] = row.data_type
                        label_to_mark[f"{label}_minus"] = row.data_type
                        label_to_type[f"{label}_plus"] = f"{row.line}_{row.tissue}"
                        label_to_type[f"{label}_minus"] = f"{row.line}_{row.tissue}"
                        label_to_track[f"{label}_plus"] = f"{row.data_type}_plus"
                        label_to_track[f"{label}_minus"] = f"{row.data_type}_minus"
            else:
                if not plot_allreps:
                    merged = f"{prefix}__merged__{row.ref_genome}"
                    onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                    bw = f"results/{row.env}/tracks/{merged}__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{strand}.bw"
                    label = f"{row.line}_{row.tissue}_{row.sample_type}"
                    grouped_bw[f"{row.data_type}"].append(bw)
                    grouped_labs[f"{row.data_type}"].append(f"{label}")
                    unique_rna.add(row.data_type)
                    label_to_mark[label] = row.data_type
                    label_to_type[label] = f"{row.line}_{row.tissue}"
                    label_to_track[label] = row.data_type
                else:
                    for rep in reps:
                        bw = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__{strand}.bw"
                        label = f"{row.line}_{row.tissue}_{row.sample_type}_{rep}"
                        grouped_bw[f"{row.data_type}"].append(bw)
                        grouped_labs[f"{row.data_type}"].append(f"{label}")
                        unique_rna.add(row.data_type)
                        label_to_mark[label] = row.data_type
                        label_to_type[label] = f"{row.line}_{row.tissue}"
                        label_to_track[label] = row.data_type
                        
        elif row.env == "sRNA":
            for size in srna_sizes:
                if strand == "unstranded":
                    if not plot_allreps:
                        merged = f"{prefix}__merged__{row.ref_genome}"
                        onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                        bw1 = f"results/{row.env}/tracks/{merged}__{size}nt__plus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__plus.bw"
                        bw2 = f"results/{row.env}/tracks/{merged}__{size}nt__minus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__minus.bw"
                        label = f"{row.line}_{row.tissue}_sRNA_{size}nt"
                        grouped_bw[f"sRNA_{size}nt"].extend([bw1, bw2])
                        grouped_labs[f"sRNA_{size}nt"].extend([f"{label}_plus", f"{label}_minus"])
                        unique_srna.add(f"sRNA_{size}nt")
                        label_to_mark[f"{label}_plus"] = f"sRNA_{size}nt"
                        label_to_mark[f"{label}_minus"] = f"sRNA_{size}nt"
                        label_to_type[f"{label}_plus"] = f"{row.line}_{row.tissue}"
                        label_to_type[f"{label}_minus"] = f"{row.line}_{row.tissue}"
                        label_to_track[f"{label}_plus"] = f"sRNA_{size}nt_plus"
                        label_to_track[f"{label}_minus"] = f"sRNA_{size}nt_minus"
                    else:
                        for rep in reps:
                            bw1 = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__{size}nt__plus.bw"
                            bw2 = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__{size}nt__minus.bw"
                            label = f"{row.line}_{row.tissue}_sRNA_{rep}_{size}nt"
                            grouped_bw[f"sRNA_{size}nt"].extend([bw1, bw2])
                            grouped_labs[f"sRNA_{size}nt"].extend([f"{label}_plus", f"{label}_minus"])
                            unique_srna.add(f"sRNA_{size}nt")
                            label_to_mark[f"{label}_plus"] = f"sRNA_{size}nt"
                            label_to_mark[f"{label}_minus"] = f"sRNA_{size}nt"
                            label_to_type[f"{label}_plus"] = f"{row.line}_{row.tissue}"
                            label_to_type[f"{label}_minus"] = f"{row.line}_{row.tissue}"
                            label_to_track[f"{label}_plus"] = f"sRNA_{size}nt_plus"
                            label_to_track[f"{label}_minus"] = f"sRNA_{size}nt_minus"
                else:
                    if not plot_allreps:
                        merged = f"{prefix}__merged__{row.ref_genome}"
                        onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                        bw = f"results/{row.env}/tracks/{merged}__{size}nt__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__{strand}.bw"
                        label = f"{row.line}_{row.tissue}_sRNA_{size}nt"
                        grouped_bw[f"sRNA_{size}nt"].append(bw)
                        grouped_labs[f"sRNA_{size}nt"].append(f"{label}")
                        unique_srna.add(f"sRNA_{size}nt")
                        label_to_mark[label] = f"sRNA_{size}nt"
                        label_to_type[label] = f"{row.line}_{row.tissue}"
                        label_to_track[label] = f"sRNA_{size}nt"
                    else:
                        for rep in reps:
                            bw = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__{size}nt__{strand}.bw"
                            label = f"{row.line}_{row.tissue}_sRNA_{rep}_{size}nt"
                            grouped_bw[f"sRNA_{size}nt"].append(bw)
                            grouped_labs[f"sRNA_{size}nt"].append(f"{label}")
                            unique_srna.add(f"sRNA_{size}nt")
                            label_to_mark[label] = f"sRNA_{size}nt"
                            label_to_type[label] = f"{row.line}_{row.tissue}"
                            label_to_track[label] = f"sRNA_{size}nt"
                        
        elif row.env == "mC":
            if not plot_allreps:
                merged = f"{prefix}__merged__{row.ref_genome}"
                onerep = f"{prefix}__{reps[0]}__{row.ref_genome}"
                for context in ["CG","CHG","CHH"]:
                    bw = f"results/{row.env}/tracks/{merged}__{context}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{context}.bw"
                    label = f"{row.line}_{row.tissue}_m{context}"
                    grouped_bw[f"m{context}"].append(bw)
                    grouped_labs[f"m{context}"].append(f"{label}")
                    unique_mc.add(f"m{context}")
                    label_to_mark[label] = f"m{context}"
                    label_to_type[label] = f"{row.line}_{row.tissue}"
                    label_to_track[label] = f"m{context}"
            else:
                for rep in reps:
                    for context in ["CG","CHG","CHH"]:
                        bw = f"results/{row.env}/tracks/{prefix}__{rep}__{row.ref_genome}__{context}.bw"
                        label = f"{row.line}_{row.tissue}_{rep}_m{context}"
                        grouped_bw[f"m{context}"].append(bw)
                        grouped_labs[f"m{context}"].append(f"{label}")
                        unique_mc.add(f"m{context}")
                        label_to_mark[label] = f"m{context}"
                        label_to_type[label] = f"{row.line}_{row.tissue}"
                        label_to_track[label] = f"m{context}"

        elif row.env == "ATAC":
            if not plot_allreps:
                merged = f"coverage__merged__{prefix}__merged__{row.ref_genome}.bw"
                onerep = f"coverage__final__{prefix}__{reps[0]}__{row.ref_genome}.bw"
                bw = f"results/ATAC/tracks/{merged}" if len(reps) >=2 else f"results/ATAC/tracks/{onerep}"
                label = f"{row.line}_{row.tissue}_{row.sample_type}"
                grouped_bw["atac"].append(bw)
                grouped_labs["atac"].append(label)
                unique_atac.add("ATAC")
                label_to_mark[label] = "ATAC"
                label_to_type[label] = f"{row.line}_{row.tissue}"
                label_to_track[label] = "ATAC"
            else:
                for rep in reps:
                    bw = f"results/ATAC/tracks/coverage__final__{prefix}__{rep}__{row.ref_genome}.bw"
                    label = f"{row.line}_{row.tissue}_{row.sample_type}_{rep}"
                    grouped_bw["atac"].append(bw)
                    grouped_labs["atac"].append(label)
                    unique_atac.add("ATAC")
                    label_to_mark[label] = "ATAC"
                    label_to_type[label] = f"{row.line}_{row.tissue}"
                    label_to_track[label] = "ATAC"

    bigwigs = (
        sum([grouped_bw.get(f"chip_{chip}", []) for chip in sorted(unique_chip)], []) +
        sum([grouped_bw.get(f"tf_{tf}", []) for tf in sorted(unique_tf)], []) +
        sum([grouped_bw.get("atac", [])], []) +
        sum([grouped_bw.get(f"{rna}", []) for rna in sorted(unique_rna)], []) +
        sum([grouped_bw.get(f"{srna}", []) for srna in sorted(unique_srna)], []) +
        sum([grouped_bw.get(f"{mc}", []) for mc in sorted(unique_mc)], [])
    )
    labels = (
        sum([grouped_labs.get(f"chip_{chip}", []) for chip in sorted(unique_chip)], []) +
        sum([grouped_labs.get(f"tf_{tf}", []) for tf in sorted(unique_tf)], []) +
        sum([grouped_labs.get("atac", [])], []) +
        sum([grouped_labs.get(f"{rna}", []) for rna in sorted(unique_rna)], []) +
        sum([grouped_labs.get(f"{srna}", []) for srna in sorted(unique_srna)], []) +
        sum([grouped_labs.get(f"{mc}", []) for mc in sorted(unique_mc)], [])
    )
    marks = ( sorted(unique_chip) + sorted(unique_tf) + sorted(unique_atac) + [f"{rna}_{strand}" for rna in sorted(unique_rna) for strand in ["plus", "minus"]] + [f"{srna}_{strand}" for srna in sorted(unique_srna) for strand in ["plus", "minus"]] + sorted(unique_mc) ) if strand == "unstranded" else ( sorted(unique_chip) + sorted(unique_tf) + sorted(unique_atac) + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )

    marksforbrowser = ( sorted(unique_chip) + sorted(unique_tf) + sorted(unique_atac) + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )
    
    types = sorted((filtered_analysis_samples["line"] + "_" + filtered_analysis_samples["tissue"]).tolist())
    
    back_palette = assign_colors(types, "tab20")
    track_palette = assign_colors(marksforbrowser, "Set2")
    plus_palette = make_it_lighter(track_palette, 1.1)
    minus_palette = make_it_lighter(track_palette, 1.4)
    for m in unique_rna | unique_srna:
        minus_palette[m] = plus_palette[m]
    
    backcolors = [back_palette[label_to_type[lab]] for lab in labels]
    trackcolors = [track_palette[label_to_mark[lab]] for lab in labels]
    fillcolorsplus = [plus_palette[label_to_mark[lab]] for lab in labels]
    fillcolorsminus = [minus_palette[label_to_mark[lab]] for lab in labels]
    alignedmarks = [label_to_track[lab] for lab in labels]
    
    if string == "bigwigs":
        return bigwigs
    elif string == "labels":
        return labels
    elif string == "marks":
        return marks
    elif string == "table":
        table_name = f"results/combined/matrix/sample_table__{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}.tab"
        os.makedirs(os.path.dirname(table_name), exist_ok=True)
        tab = pd.DataFrame({
            "bigwigs": bigwigs,
            "labels": labels,
            "backcolors": backcolors,
            "trackcolors": trackcolors,
            "fillcolorsplus": fillcolorsplus,
            "fillcolorsminus": fillcolorsminus,
            "marks": alignedmarks,
        })
        tab.to_csv(table_name, sep="\t", index=False, header=False)
        return table_name

def define_individual_browser_plots(wildcards):
    ref_genome = wildcards.ref_genome
    analysis_name = wildcards.analysis_name    
    env = wildcards.env
    target_name = wildcards.target_name
    target_file = define_combined_target_file(wildcards)
    
    checkpoints.is_stranded.get(bedfile=target_file)
    
    with open(target_file) as f:
        first_line = f.readline().strip().split("\t")
    try:
        header = "no" if (int(first_line[1]) >=0 and int(first_line[2]) >=0) else "yes"
    except (ValueError, IndexError):
        header = "yes"

    regions = pd.read_csv(target_file, sep="\t", header=0 if header=="yes" else None)
    
    files = []
    starter = 2 if header == "yes" else 1
    for row_number, _ in enumerate(regions.itertuples(index=False),start=starter):
        files.append(f"results/combined/plots/single_browser__{target_name}__line{row_number}__{env}__{analysis_name}__{ref_genome}.pdf")    
    return files

def define_input_for_pca(wildcards, string):
    tracks = []
    labels = []
    unique_group = set()
    label_to_group = {}
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    
    if globenv in ["mCG", "mCHG", "mCHH"]:
        filtered_samples = samples[ (samples['env'] == "mC") & (samples['ref_genome'] == ref_genome) ].copy()
        context = globenv[1:]
        for _, row in filtered_samples.iterrows():
            bw = f"results/mC/tracks/{sample_name_str(row, 'sample')}__{context}.bw"
            label = f"{row.line}_{row.tissue}_{row.replicate}"
            group = f"{row.line}_{row.tissue}"
            tracks.append(bw)
            labels.append(label)
            unique_group.add(group)
            label_to_group[label] = group
    elif globenv in ["TF", "ChIP", "ATAC"]:
        filtered_samples = samples[ (samples['env'] == globenv) & (samples['ref_genome'] == ref_genome) ].copy()
        for _, row in filtered_samples.iterrows():
            bam = f"results/{globenv}/mapped/final__{sample_name_str(row, 'sample')}.bam"
            label = f"{row.sample_type}_{row.line}_{row.tissue}_{row.data_type}_{row.replicate}"
            group = f"{row.sample_type}_{row.line}_{row.tissue}"
            tracks.append(bam)
            labels.append(label)
            unique_group.add(group)
            label_to_group[label] = group
    elif globenv == "all_chip":
        filtered_samples = samples[ (samples['env'].isin(["ChIP","TF","ATAC"])) & (samples['ref_genome'] == ref_genome) ].copy()
        for _, row in filtered_samples.iterrows():
            bam = f"results/{row.env}/mapped/final__{sample_name_str(row, 'sample')}.bam"
            label = f"{row.sample_type}_{row.line}_{row.tissue}_{row.data_type}_{row.replicate}"
            group = f"{row.sample_type}_{row.line}_{row.tissue}"
            tracks.append(bam)
            labels.append(label)
            unique_group.add(group)
            label_to_group[label] = group
    
    palette = assign_colors(unique_group, "tab20")
    colors = [palette[label_to_group[lab]] for lab in labels]
    
    if string == "tracks": 
        return tracks
    elif string == "labels":
        return labels
    elif string == "colors":
        return colors

def define_final_stats_output():
    aligned_bams = config['aligned_bams']
    stat_files = []    
    if not aligned_bams:
        stat_files += expand("results/combined/plots/mapping_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","TF","ATAC"]])

    stat_files += expand("results/combined/plots/mapping_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["RNA","mC"]])
    stat_files += expand("results/combined/plots/srna_sizes_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["sRNA"]])
    stat_files += expand("results/combined/plots/peak_stats_{analysis_name}_{env}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","TF","ATAC"]])

    return stat_files

def define_final_combined_output(ref_genome):
    qc_option = config["QC_option"]
    full_analysis = config['full_analysis']
    te_analysis = config['te_analysis']
    analysis_name = config['analysis_name']
    mc_sort = config['heatmap_sort_mc_after_others']    
    mc_context = config['mC_context']
    plot_files = []
    te_plots = []
    
    all_analysis_samples = analysis_samples[ analysis_samples['ref_genome'] == ref_genome ].copy()
    chip_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    tf_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'TF') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    atac_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ATAC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    mc_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'mC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    rna_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'RNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    srna_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'sRNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    rampage_analysis_samples = analysis_samples[ (analysis_samples['data_type'] == 'RAMPAGE') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
        
    if len(all_analysis_samples) >=1:
        plot_files.append(f"results/combined/plots/Browser_full_chromosomes__all__{analysis_name}__{ref_genome}.pdf")
    
    if len(chip_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__ChIP__{analysis_name}__{ref_genome}.pdf")
    
    if len(tf_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__TF__{analysis_name}__{ref_genome}.pdf")
    
    peak_envs = sum(1 for s in [chip_analysis_samples, tf_analysis_samples, atac_analysis_samples] if len(s) >= 1)
    if peak_envs >= 2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__all_chip__{analysis_name}__{ref_genome}.pdf")
    
    if len(srna_analysis_samples) >=1:
        plot_files.append(f"results/combined/plots/Upset_combined_clusters__sRNA__{analysis_name}__{ref_genome}.pdf")

    if len(atac_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_peaks__ATAC__{analysis_name}__{ref_genome}.pdf")
    
    if len(rampage_analysis_samples) >=2:
        plot_files.append(f"results/combined/plots/Upset_combined_TSS__RAMPAGE__{analysis_name}__{ref_genome}.pdf")

    if len(mc_analysis_samples) >=1:
        if len(all_analysis_samples) > len(mc_analysis_samples) and mc_sort:
            plot_files.append(f"results/combined/plots/Heatmap_sorted__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap_sorted__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap_sorted__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")            
            te_plots.append(f"results/combined/plots/Heatmap_sorted__regions__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"results/combined/plots/Heatmap_sorted__tss__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"results/combined/plots/Heatmap_sorted__tes__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            
            plot_files.append(f"results/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            te_plots.append(f"results/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"results/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"results/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        else:
            plot_files.append(f"results/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            te_plots.append(f"results/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"results/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"results/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        
        plot_files.append(f"results/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"results/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"results/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"results/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        
    else:
        plot_files.append(f"results/combined/plots/Heatmap__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Heatmap__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Heatmap__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"results/combined/plots/Heatmap__regions__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"results/combined/plots/Heatmap__tss__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"results/combined/plots/Heatmap__tes__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
    
        plot_files.append(f"results/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"results/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"results/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"results/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"results/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_TEs.pdf")

    mc_rep_samples = samples[ (samples['env'] == "mC") & (samples['ref_genome'] == ref_genome) ].copy()
    if len(mc_rep_samples) >=3:
        plot_files.append(f"results/combined/plots/PCA__mCG__{analysis_name}__{ref_genome}.pdf")
        if mc_context == "all":
            plot_files.append(f"results/combined/plots/PCA__mCHG__{analysis_name}__{ref_genome}.pdf")
            plot_files.append(f"results/combined/plots/PCA__mCHH__{analysis_name}__{ref_genome}.pdf")
    
    all_chip_samples = samples[ (samples['env'].isin(["ChIP","TF","ATAC"])) & (samples['ref_genome'] == ref_genome) ].copy()
    if len(all_chip_samples) >=3:
            plot_files.append(f"results/combined/plots/PCA__all_chip__{analysis_name}__{ref_genome}.pdf")
    
    for env in [e for e in UNIQUE_ENVS if e in ["TF", "ChIP", "ATAC"]]:
        env_rep_samples = samples[ (samples['env'] == env) & (samples['ref_genome'] == ref_genome) ].copy()
        if len(env_rep_samples) >=3:
            plot_files.append(f"results/combined/plots/PCA__{env}__{analysis_name}__{ref_genome}.pdf")

    results = []
    
    if full_analysis:
        results += plot_files
        
    if te_analysis:
        results += te_plots   
    
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
        
        has_strand = False
        strand_values = set()
        with open(input.bedfile) as f:
            if header == "yes":
                next(f)
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 6:
                    has_strand = True
                    strand_values.add(cols[5])

        with open(output.file, "w") as out:
            if has_strand and strand_values.issubset({"+", "-"}):
                out.write("stranded\n")
            else:
                out.write("unstranded\n")

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
    threads: config["resources"]["prepping_mapping_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["prepping_mapping_stats"]["mem_mb"],
        tmp_mb=config["resources"]["prepping_mapping_stats"]["tmp_mb"],
        qos=config["resources"]["prepping_mapping_stats"]["qos"]
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
        mem_mb=config["resources"]["plotting_mapping_stats"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_mapping_stats"]["tmp_mb"],
        qos=config["resources"]["plotting_mapping_stats"]["qos"]
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
    threads: config["resources"]["prepping_chip_peak_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["prepping_chip_peak_stats"]["mem_mb"],
        tmp_mb=config["resources"]["prepping_chip_peak_stats"]["tmp_mb"],
        qos=config["resources"]["prepping_chip_peak_stats"]["qos"]
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
        mem_mb=config["resources"]["plotting_peaks_stats_chip_tf"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_peaks_stats_chip_tf"]["tmp_mb"],
        qos=config["resources"]["plotting_peaks_stats_chip_tf"]["qos"]
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
    threads: config["resources"]["prepping_srna_sizes_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["prepping_srna_sizes_stats"]["mem_mb"],
        tmp_mb=config["resources"]["prepping_srna_sizes_stats"]["tmp_mb"],
        qos=config["resources"]["prepping_srna_sizes_stats"]["qos"]
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
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_sizes_stats.R"),
        zoommin=config['srna_min_size'],
        zoommax=config['srna_max_size']
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_srna_sizes"))
    conda: CONDA_ENV
    threads: config["resources"]["plotting_srna_sizes_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["plotting_srna_sizes_stats"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_srna_sizes_stats"]["tmp_mb"],
        qos=config["resources"]["plotting_srna_sizes_stats"]["qos"]
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}" "{params.zoommin}" "{params.zoommax}"
        """

###
# Rules to prep and plot ChIP upset plots
rule combine_clusterfiles:
    input:
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        clusterfiles = lambda wildcards: define_samples_for_upset(wildcards, "files")
    output:
        temp1_file = temp("results/combined/bedfiles/temp1_combined_clusters__{env}__{analysis_name}__{ref_genome}.bed"),
        temp2_file = temp("results/combined/bedfiles/temp2_combined_clusters__{env}__{analysis_name}__{ref_genome}.bed"),
        merged_file = "results/combined/bedfiles/combined_clusters__{env}__{analysis_name}__{ref_genome}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samples_for_upset(wildcards, "pairs"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_clusters_{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["combine_clusterfiles"]["threads"]
    resources:
        mem_mb=config["resources"]["combine_clusterfiles"]["mem_mb"],
        tmp_mb=config["resources"]["combine_clusterfiles"]["tmp_mb"],
        qos=config["resources"]["combine_clusterfiles"]["qos"]
    shell:
        """
        {{
        printf "Merging peakfiles for {params.env} {params.analysis_name} {params.ref_genome}\n"
        for pair in {params.names}; do
            label=$(echo ${{pair}} | cut -d":" -f1)
            file=$(echo ${{pair}} | cut -d":" -f2)
            awk -v OFS="\t" -v l=${{label}} '{{print $1,$2,$3,l"_"$4}}' ${{file}} >> {output.temp1_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        printf "Chr\tStart\tStop\tClusterID\tSamples\n" > {output.merged_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v e={params.env} -v a={params.analysis_name} '{{print $1,$2,$3,"combined_clusters_"e"_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule combine_peakfiles:
    input:
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        peakfiles = lambda wildcards: define_samples_for_upset(wildcards, "files")
    output:
        temp1_file = temp("results/combined/bedfiles/temp1_combined_peaks__{env}__{analysis_name}__{ref_genome}.bed"),
        temp2_file = temp("results/combined/bedfiles/temp2_combined_peaks__{env}__{analysis_name}__{ref_genome}.bed"),
        merged_file = "results/combined/bedfiles/combined_peaks__{env}__{analysis_name}__{ref_genome}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samples_for_upset(wildcards, "pairs"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_peaks_{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["combine_peakfiles"]["threads"]
    resources:
        mem_mb=config["resources"]["combine_peakfiles"]["mem_mb"],
        tmp_mb=config["resources"]["combine_peakfiles"]["tmp_mb"],
        qos=config["resources"]["combine_peakfiles"]["qos"]
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
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v e={params.env} -v a={params.analysis_name} '{{print $1,$2,$3,"combined_peaks_"e"_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule combine_TSS:
    input:
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        peakfiles = lambda wildcards: define_samples_for_upset(wildcards, "files")
    output:
        temp1_file = temp("results/combined/bedfiles/temp1_combined_TSS__{env}__{analysis_name}__{ref_genome}.bed"),
        temp2_file = temp("results/combined/bedfiles/temp2_combined_TSS__{env}__{analysis_name}__{ref_genome}.bed"),
        merged_file = "results/combined/bedfiles/combined_TSS__{env}__{analysis_name}__{ref_genome}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samples_for_upset(wildcards, "pairs"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_TSS_{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["combine_TSS"]["threads"]
    resources:
        mem_mb=config["resources"]["combine_TSS"]["mem_mb"],
        tmp_mb=config["resources"]["combine_TSS"]["tmp_mb"],
        qos=config["resources"]["combine_TSS"]["qos"]
    shell:
        """
        {{
        printf "Merging TSS files for {params.env} {params.analysis_name} {params.ref_genome}\n"
        for pair in {params.names}; do
            label=$(echo ${{pair}} | cut -d":" -f1)
            file=$(echo ${{pair}} | cut -d":" -f2)
            awk -v OFS="\t" -v l=${{label}} '{{print $1,$2,$3,l}}' ${{file}} >> {output.temp1_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        printf "Chr\tStart\tStop\tTSSID\tSamples\n" > {output.merged_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v a={params.analysis_name} '{{print $1,$2,$3,"combined_TSS_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """

rule get_annotations_for_bedfile:
    input:
        bedfile = lambda wildcards: define_combined_target_file(wildcards),
        region_file = lambda wildcards: f"results/combined/bedfiles/{wildcards.ref_genome}__all_genes.bed",
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
        mem_mb=config["resources"]["get_annotations_for_bedfile"]["mem_mb"],
        tmp_mb=config["resources"]["get_annotations_for_bedfile"]["tmp_mb"],
        qos=config["resources"]["get_annotations_for_bedfile"]["qos"]
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
        printf "Chr\tStart\tStop\tRegionID\tDistance\tGene_strand\tGID\tCategory\n" > {output.annotated_file}
        bedtools closest -a {output.temp_bedfile} -b {input.region_file} -g {input.chrom_sizes} -D ref | awk -v OFS="\t" '{{if ($10=="+") print $1,$2,$3,$4,$11,$10,$8; else print $1,$2,$3,$4,-$11,$10,$8}}' | awk -F"[=;]" -v OFS="\t" '{{print $1,$2}}' | sed 's/gene://' | awk -v OFS="\t" '{{if ($5<-2000) {{d="Distal_downstream"}} else if ($5<0) {{d="Terminator"}} else if ($5==0) {{d="Gene_body"}} else if ($5>2000) {{d="Distal_upstream"}} else {{d="Promoter"}}; print $1,$2,$3,$4,$5,$6,$8,d}}' >> {output.annotated_file}
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_upset_regions:
    input:
        mergedfile = "results/combined/bedfiles/{target_name}__{env}__{analysis_name}__{ref_genome}.bed",
        annotatedfile = "results/combined/bedfiles/annotated__{target_name}__{env}__{analysis_name}__{ref_genome}.bed"
    output:
        plot = "results/combined/plots/Upset_{target_name}__{env}__{analysis_name}__{ref_genome}.pdf"
    params:
        env = lambda wildcards: wildcards.env,
        types = lambda wildcards: define_samples_for_upset(wildcards, "types"),
        script = lambda wildcards: define_upset_script(wildcards)
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_upset_{target_name}_{env}"))
    conda: CONDA_ENV_UPSET
    threads: config["resources"]["plotting_upset_regions"]["threads"]
    resources:
        mem_mb=config["resources"]["plotting_upset_regions"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_upset_regions"]["tmp_mb"],
        qos=config["resources"]["plotting_upset_regions"]["qos"]
    shell:
        """
        Rscript "{params.script}" "{input.mergedfile}" "{input.annotatedfile}" "{params.env}" "{params.types}" "{output.plot}"
        """

###
# rule to plot heatmaps
rule making_stranded_matrix_on_targetfile:
    input:
        bigwigs = lambda wildcards: define_key_for_plots(wildcards, "bigwigs"),
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
        labels = lambda wildcards: define_key_for_plots(wildcards, "labels"),
        marks = lambda wildcards: define_key_for_plots(wildcards, "marks"),
        matrix = lambda wildcards: wildcards.matrix_param,
        strand = lambda wildcards: wildcards.strand,
        base = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'base'),
        bs = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'bs'),
        base_mc = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'base_mc'),
        bs_mc = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'bs_mc'),
        before = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'before'),
        after = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'after'),
        middle = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'middle')
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "making_matrix_{matrix_param}_{target_name}_{strand}"))
    conda: CONDA_ENV
    threads: config["resources"]["making_stranded_matrix_on_targetfile"]["threads"]
    resources:
        mem_mb=config["resources"]["making_stranded_matrix_on_targetfile"]["mem_mb"],
        tmp_mb=config["resources"]["making_stranded_matrix_on_targetfile"]["tmp_mb"],
        qos=config["resources"]["making_stranded_matrix_on_targetfile"]["qos"]
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
        if [[ "{params.env}" == "mC" ]]; then
            computeMatrix {params.base_mc} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.labels} -bs {params.bs_mc} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        else
            computeMatrix {params.base} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.labels} -bs {params.bs} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        fi
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
        mem_mb=config["resources"]["merging_matrix"]["mem_mb"],
        tmp_mb=config["resources"]["merging_matrix"]["tmp_mb"],
        qos=config["resources"]["merging_matrix"]["qos"]
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
        params_heatmap = temp("results/combined/matrix/params_heatmap_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"),
        params_profile = temp("results/combined/matrix/params_profile_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"),
        params_regions = temp("results/combined/matrix/params_regions_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"),
        temp_values = temp("results/combined/matrix/temp_values_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"),
        temp_profile = temp("results/combined/matrix/temp_profile_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf"),
        temp_profile_values = temp("results/combined/matrix/temp_profile_values_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        scales = config['heatmaps_scales'],
        profile = config['profiles_scale'],
        cg_scale = config['heat_mcg'],
        chg_scale = config['heat_mchg'],
        chh_scale = config['heat_mchh']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "getting_scales_matrix_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["computing_matrix_scales"]["threads"]
    resources:
        mem_mb=config["resources"]["computing_matrix_scales"]["mem_mb"],
        tmp_mb=config["resources"]["computing_matrix_scales"]["tmp_mb"],
        qos=config["resources"]["computing_matrix_scales"]["qos"]
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
                zmini=$(grep "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=999999}} {{a=$5; if (a != "nan" && a<m) m=a;}} END {{if (m != 999999) print m; else print 0}}')
                zmaxi=$(grep "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=-999999}} {{a=$6; if (a != "nan" && a>m) m=a;}} END {{if (m != -999999) print m; else print 0}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" && ${{mark}} == "mCG" ]]; then
                    zmini="0"
                    zmaxi="{params.cg_scale}"
                elif [[ "${{test}}" == "yes" && ${{mark}} == "mCHG" ]]; then
                    zmini="0"
                    zmaxi="{params.chg_scale}"
                elif [[ "${{test}}" == "yes" && ${{mark}} == "mCHH" ]]; then
                    zmini="0"
                    zmaxi="{params.chh_scale}"
                elif [[ "${{test}}" == "yes" ]]; then
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
                zmini=$(grep "${{sample}}" {output.temp_values} | awk '{{if ($5 != "nan") print $5; else print 0}}')
                zmaxi=$(grep "${{sample}}" {output.temp_values} | awk '{{if ($6 != "nan") print $6; else print 0}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" && "${{sample}}" =~ mCG ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.cg_scale}")
                elif [[ "${{test}}" == "yes" && "${{sample}}" =~ mCHG ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.chg_scale}")
                elif [[ "${{test}}" == "yes" && "${{sample}}" =~ mCHH ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.chh_scale}")
                elif [[ "${{test}}" == "yes" ]]; then
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
        mem_mb=config["resources"]["plotting_heatmap_on_targetfile"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_heatmap_on_targetfile"]["tmp_mb"],
        qos=config["resources"]["plotting_heatmap_on_targetfile"]["qos"]
    shell:
        """
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        reg="$(cat {input.params_regions})"
        heat="$(cat {input.params_heatmap})"
        prof="$(cat {input.params_profile})"
        printf "Plotting heatmap {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        plotHeatmap -m {input.matrix} -out {output.plot} {params.plot_params} {params.sort} ${{reg}} ${{heat}} ${{prof}} ${{add}} --outFileSortedRegions {output.sorted_regions}
        """

rule sort_heatmap:
    input: 
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__mC__{analysis_name}__{ref_genome}__{target_name}.gz",
        sorted_regions = "results/combined/matrix/Heatmap__{matrix_param}__most__{analysis_name}__{ref_genome}__{target_name}_sorted_regions.bed",
        params_regions = "results/combined/matrix/params_regions_final_matrix_{matrix_param}__most__{analysis_name}__{ref_genome}__{target_name}.txt"
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
        mem_mb=config["resources"]["sort_heatmap"]["mem_mb"],
        tmp_mb=config["resources"]["sort_heatmap"]["tmp_mb"],
        qos=config["resources"]["sort_heatmap"]["qos"]
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
        mem_mb=config["resources"]["plotting_sorted_heatmap_on_targetfile"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_sorted_heatmap_on_targetfile"]["tmp_mb"],
        qos=config["resources"]["plotting_sorted_heatmap_on_targetfile"]["qos"]
    shell:
        """
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        heat="$(cat {input.params_heatmap})"
        prof="$(cat {input.params_profile})"
        printf "Plotting heatmap {params.matrix} for mC {params.target_name} on {params.ref_genome}\n"
        plotHeatmap -m {input.matrix} -out {output.plot} {params.plot_params} --sortRegions 'keep' ${{heat}} ${{prof}} ${{add}}
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
        mem_mb=config["resources"]["plotting_profile_on_targetfile"]["mem_mb"],
        tmp_mb=config["resources"]["plotting_profile_on_targetfile"]["tmp_mb"],
        qos=config["resources"]["plotting_profile_on_targetfile"]["qos"]
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
        reg="$(cat {input.params_regions})"
        prof="$(cat {input.params_profile})"
        plotProfile -m {input.matrix} -out {output.plot1} {params.plot_params} ${{reg}} ${{prof}} ${{add}}
        
        printf "Plotting per group profile {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        ymin=$(cat {input.params_profile} | awk 'BEGIN {{y=99999}} {{for (i=1; i<=NF; i++) {{if ($i == "--yMin") {{for (j=i+1; j<=NF && $j !~ /^--/; j++) {{if ($j<y) y=$j}} break}} }} }} END {{print y}}' )
        ymax=$(cat {input.params_profile} | awk 'BEGIN {{y=-99999}} {{for (i=1; i<=NF; i++) {{if ($i == "--yMax") {{for (j=i+1; j<=NF && $j !~ /^--/; j++) {{if ($j>y) y=$j}} break}} }} }} END {{print y}}' )
        plotProfile -m {input.matrix} -out {output.plot2} {params.plot_params} ${{reg}} --yMin ${{ymin}} --yMax ${{ymax}} ${{add}} --perGroup
        }} 2>&1 | tee -a "{log}"
        """

###
# rules to plot browser shots
rule prep_chromosomes_for_browser:
    input: 
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes"
    output:
        bedfile = "results/combined/bedfiles/full_chromosomes__{ref_genome}.bed"
    params:
        chromosome_bs = config['chromosome_bs']
    log:
        temp(return_log_combined("bedfile", "{ref_genome}", "prep_chromosomes"))
    conda: CONDA_ENV
    threads: config["resources"]["prep_chromosomes_for_browser"]["threads"]
    resources:
        mem_mb=config["resources"]["prep_chromosomes_for_browser"]["mem_mb"],
        tmp_mb=config["resources"]["prep_chromosomes_for_browser"]["tmp_mb"],
        qos=config["resources"]["prep_chromosomes_for_browser"]["qos"]
    shell:
        """
        {{
        awk -v OFS="\t" -v c={params.chromosome_bs} 'NR <= 50 {{if ($2/c > 1) b=c; else b=int($2/500); print $1,"1",$2,$1,b}}' {input.chrom_sizes} > {output.bedfile}
        }} 2>&1 | tee -a "{log}" 
        """
        
rule prep_browser_on_region:
    input:
        bigwigs = lambda wildcards: define_key_for_plots(wildcards, "bigwigs"),
        target_file = lambda wildcards: define_combined_target_file(wildcards),
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        gff = lambda wildcards: f"genomes/{wildcards.ref_genome}/{wildcards.ref_genome}.gff",
        all_genes = lambda wildcards: f"results/combined/bedfiles/{wildcards.ref_genome}__all_genes.bed",
        TE_file = lambda wildcards: f"genomes/{wildcards.ref_genome}/{wildcards.ref_genome}__TE_file.bed" if config['browser_TE_file'] else []
    output:
        filenames = "results/combined/matrix/filenames__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        genes = "results/combined/matrix/genes_in_locus__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.gff",
        tes = "results/combined/matrix/tes_in_locus__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.bed",
        htstart = "results/combined/matrix/highlight_start__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        htwidth = "results/combined/matrix/highlight_width__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        name = "results/combined/matrix/name__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        tempgenes = temp("results/combined/matrix/genes_in_locus__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt"),
        templocus = temp("results/combined/matrix/locus_{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.bed"),
        temparray = temp("results/combined/matrix/array__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.npz"),
        tempvalues = temp("results/combined/matrix/values__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.tab")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        labels = lambda wildcards: define_key_for_plots(wildcards, "labels"),
        sample_table = lambda wildcards: define_key_for_plots(wildcards, "table"),
        trackfolder = lambda wildcards: f"results/combined/matrix/tracks_{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}",
        regionID = lambda wildcards: wildcards.regionID,
        browser_scales = config['browser_scales'],
        mc_scales = config['fixed_mc_scales'],
        cg_scale = config['fixed_mcg'],
        chg_scale = config['fixed_mchg'],
        chh_scale = config['fixed_mchh']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "prep_files_{target_name}_{regionID}"))
    conda: CONDA_ENV
    threads: config["resources"]["prep_browser_on_region"]["threads"]
    resources:
        mem_mb=config["resources"]["prep_browser_on_region"]["mem_mb"],
        tmp_mb=config["resources"]["prep_browser_on_region"]["tmp_mb"],
        qos=config["resources"]["prep_browser_on_region"]["qos"]
    shell:
        """
        {{
        printf "Extracting values for {params.regionID}\n"
        tmpregion="{params.regionID}"
        line_nb=$(echo "${{tmpregion}}" | sed 's/line//')
        chr=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $1}}')
        start=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $2}}')
        end=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $3}}')
        printf "${{chr}}\t${{start}}\t${{end}}\n" > {output.templocus}
        region="${{chr}}:${{start}}:${{end}}"

        regionID=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $4}}')
        printf "${{regionID}}\n" > {output.name}
        
        binsize=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $5}}')
        
        htstart=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $6}}')
        htwidth=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $7}}')
        if [[ ${{htstart}} != "" ]]; then
            printf "${{htstart}}\n" | awk -F"," '{{for (i=1;i<=NF;i++) print $i}}' > "{output.htstart}"
            printf "${{htwidth}}\n" | awk -F"," '{{for (i=1;i<=NF;i++) print $i}}' > "{output.htwidth}"
        else
            touch {output.htstart}
            touch {output.htwidth}
        fi
        
        ### To get genes in the region
        bedtools intersect -a {input.all_genes} -b {output.templocus} | awk '{{print $4}}' > {output.tempgenes}
        if [[ "{params.target_name}" == "full_chromosomes" ]]; then
            printf "Do not include genes in whole chromosomes\n"
            touch {output.genes}
        elif [[ -s "{output.tempgenes}" ]]; then
            printf "Getting gene track\n"
            bedtools intersect -wb -a {input.gff} -b {output.templocus} | awk -v OFS="\t" '{{if ($7!="+" && $7!="-") $7="*"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' > {output.genes}
        else
            printf "No genes in this region\n"
            touch {output.genes}
        fi
        ### To get the bed files of TEs. For now relying on a bed file of TEs (only one, needing to match the species).
        if [[ -n {input.TE_file} && -s {input.TE_file} && "{params.target_name}" != "full_chromosomes" ]]; then
            printf "Getting TE track\n"
            bedtools intersect -a {input.TE_file} -b {output.templocus} | awk -v OFS="\t" '{{if ($6!="+" && $6!="-") $6="*"; print $0}}' > {output.tes}
        elif [[ -n {input.TE_file} ]]; then
            printf "No TE file to be included in browser\n"
            touch {output.tes}
        fi
        
        printf "Testing if all bw have data in this region\n"
        filelist2=()
        rm -rf "{params.trackfolder}" && mkdir "{params.trackfolder}"
        while read bw lab back track plus minus mark
        do
            path="{params.trackfolder}/${{lab}}_empty"
            bigWigToBedGraph -chrom=${{chr}} -start=${{start}} -end=${{end}} ${{bw}} "${{path}}.bg"
            if [[ -s "${{path}}.bg" ]]; then
                printf "${{lab}} has data on ${{chr}}:${{start}}-${{end}}\n"
                filelist2+=("${{bw}}")
            else
                printf "${{lab}} is empty on ${{chr}}:${{start}}-${{end}}\n"
                grep "${{chr}}" {input.chrom_sizes} | awk -v OFS="\t" '{{print $1,"1",$2,"0"}}' | bedtools sort -i - > "${{path}}.bg"
                bedGraphToBigWig "${{path}}.bg" {input.chrom_sizes} "${{path}}.bw"
                filelist2+=("${{path}}.bw")
            fi
            rm -f "${{path}}.bg"            
        done < {params.sample_table}
        printf "Summarize bigwigs in binsize of ${{binsize}} bp on {params.regionID}\n"
        multiBigwigSummary bins -b ${{filelist2[@]}} -l {params.labels} -r ${{region}} -p {threads} -bs ${{binsize}} -out {output.temparray} --outRawCounts {output.tempvalues}
        
        while read bw lab back track plus minus mark
        do
            path="{params.trackfolder}/${{lab}}_${{mark}}"
            printf "Making bw for ${{lab}}\n"
            col=($(awk -v ORS=" " -v t=${{lab}} 'NR==1 {{for (i=1;i<=NF;i++) if ($i~t) print i}}' {output.tempvalues}))
            if [[ "${{lab}}" == *_minus ]]; then
                awk -v OFS="\t" -v a=${{col}} 'NR>1 {{if ($a == "nan") b=0; else b=-$a; print $1,$2,$3,b}}' {output.tempvalues} | bedtools sort -g {input.chrom_sizes} > "${{path}}.bedGraph"
            else
                awk -v OFS="\t" -v a=${{col}} 'NR>1 {{if ($a == "nan") b=0; else b=$a; print $1,$2,$3,b}}' {output.tempvalues} | bedtools sort -g {input.chrom_sizes} > "${{path}}.bedGraph"
            fi
            bedGraphToBigWig "${{path}}.bedGraph" {input.chrom_sizes} "${{path}}.bw"
        done < {params.sample_table}
        
        printf "Name\tPath\tBackcolor\tTrackcolor\tFillcolorplus\tFillcolorminus\tYmin\tYmax\n" > {output.filenames}
        while read bw lab back track plus minus mark
        do
            path="{params.trackfolder}/${{lab}}_${{mark}}"
            if [[ {params.browser_scales} == "sample" ]]; then
                ymin=$(cat "${{path}}.bedGraph" | awk 'BEGIN {{a=9999}} {{if ($4<a) a=$4;}} END {{if (a<0) b=a*1.2; else b=a*0.8; print b}}')
                ymax=$(cat "${{path}}.bedGraph" | awk 'BEGIN {{a=-9999}} {{if ($4>a) a=$4;}} END {{if (a>0) b=a*1.2; else b=a*0.8; print b}}')
            elif [[ {params.browser_scales} == "type" ]]; then
                ymin=$(cat {params.trackfolder}/*_${{mark}}.bedGraph | awk 'BEGIN {{a=9999}} {{if ($4<a) a=$4;}} END {{if (a<0) b=a*1.2; else b=a*0.8; print b}}')
                ymax=$(cat {params.trackfolder}/*_${{mark}}.bedGraph | awk 'BEGIN {{a=-9999}} {{if ($4>a) a=$4;}} END {{if (a>0) b=a*1.2; else b=a*0.8; print b}}')
            fi
            if [[ {params.mc_scales} == "True" ]]; then
                if [[ ${{mark}} == "mCG" ]]; then
                    ymin=0
                    ymax={params.cg_scale}
                elif [[ ${{mark}} == "mCHG" ]]; then
                    ymin=0
                    ymax={params.chg_scale}
                elif [[ ${{mark}} == "mCHH" ]]; then
                    ymin=0
                    ymax={params.chh_scale}
                fi
            fi
            printf "${{lab}}\t${{path}}.bw\t${{back}}\t${{track}}\t${{plus}}\t${{minus}}\t${{ymin}}\t${{ymax}}\n" >> {output.filenames}
        done < {params.sample_table}
        rm -f {params.trackfolder}/*.bedGraph
        
        }} 2>&1 | tee -a "{log}" 
        """

rule make_single_loci_browser_plot:
    input:
        filenames = "results/combined/matrix/filenames__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        genes = "results/combined/matrix/genes_in_locus__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.gff",
        tes = "results/combined/matrix/tes_in_locus__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.bed",
        name = "results/combined/matrix/name__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        htstart = "results/combined/matrix/highlight_start__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt",
        htwidth = "results/combined/matrix/highlight_width__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.txt"
    output:
        browser_plot = temp("results/combined/plots/single_browser__{target_name}__{regionID}__{env}__{analysis_name}__{ref_genome}.pdf")
    params:
        regionID = lambda wildcards: wildcards.regionID,
        title = lambda wildcards: f"results/combined/plots/single_browser__{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}.pdf",
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_browser_plot.R"),
        trackfolder = lambda wildcards: f"results/combined/matrix/tracks_{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}"
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "single_browser_{target_name}_{regionID}"))
    conda: CONDA_ENV
    threads: config["resources"]["make_single_loci_browser_plot"]["threads"]
    resources:
        mem_mb=config["resources"]["make_single_loci_browser_plot"]["mem_mb"],
        tmp_mb=config["resources"]["make_single_loci_browser_plot"]["tmp_mb"],
        qos=config["resources"]["make_single_loci_browser_plot"]["qos"]
    shell:
        """
        {{
        name=$(cat {input.name})
        if [[ -s {input.htstart} ]]; then
            printf "\nPlotting browser on {params.regionID} with higlights\n\n"
            Rscript "{params.script}" "{input.filenames}" "{input.genes}" "{input.tes}" "${{name}}" "{params.title}" "{input.htstart}" "{input.htwidth}"
        else
            printf "\nPlotting browser on {params.regionID} without higlights\n\n"
            Rscript "{params.script}" "{input.filenames}" "{input.genes}" "{input.tes}" "${{name}}" "{params.title}"
        fi
        rm -rf {params.trackfolder}
        }} 2>&1 | tee -a "{log}"
        """

rule merge_region_browser_plots:
    input: 
        plots = lambda wildcards: define_individual_browser_plots(wildcards)
    output:
        merged_plots = "results/combined/plots/Browser_{target_name}__{env}__{analysis_name}__{ref_genome}.pdf"
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "merging_browser_{target_name}"))
    conda: CONDA_ENV
    threads: config["resources"]["merge_region_browser_plots"]["threads"]
    resources:
        mem_mb=config["resources"]["merge_region_browser_plots"]["mem_mb"],
        tmp_mb=config["resources"]["merge_region_browser_plots"]["tmp_mb"],
        qos=config["resources"]["merge_region_browser_plots"]["qos"]
    shell:
        """
        pdfunite {input.plots} {output.merged_plots}
        """

############## To plot PCA

rule summarize_tracks_pca:
    input:
        tracks = lambda wildcards: define_input_for_pca(wildcards, "tracks")
    output:
        array = "results/combined/matrix/pca_matrix__{env}__{analysis_name}__{ref_genome}.npz"
    params:
        labels = lambda wildcards: define_input_for_pca(wildcards, "labels"),
        bs = config['pca_bs'],
        step = config['pca_step']
    log:
        temp(return_log_mc("{analysis_name}__{ref_genome}", "summarize_tracks_pca", "{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["summarize_tracks_pca"]["threads"]
    resources:
        mem_mb=config["resources"]["summarize_tracks_pca"]["mem_mb"],
        tmp_mb=config["resources"]["summarize_tracks_pca"]["tmp_mb"],
        qos=config["resources"]["summarize_tracks_pca"]["qos"]
    shell:
        """
        {{
        if [[ {wildcards.env} == "mCG" || {wildcards.env} == "mCHG" || {wildcards.env} == "mCHH" ]]; then
            printf "Summarizing bigwigs for {wildcards.analysis_name} {wildcards.ref_genome} for mC samples in {wildcards.env} sequence context\n"
            multiBigwigSummary bins -b {input.tracks} -o {output.array} -l {params.labels} -bs {params.bs} -n {params.step} -p {threads}
        else
            printf "Summarizing bams for {wildcards.analysis_name} {wildcards.ref_genome} for {wildcards.env} samples\n"
            multiBamSummary bins -b {input.tracks} -o {output.array} -l {params.labels} -bs {params.bs} -n {params.step} -p {threads}
        fi
        }} 2>&1 | tee -a "{log}"
        """    

rule plot_PCA_correlation:
    input:
        array = "results/combined/matrix/pca_matrix__{env}__{analysis_name}__{ref_genome}.npz"
    output:
        plot = "results/combined/plots/PCA__{env}__{analysis_name}__{ref_genome}.pdf"
    params:
        colors = lambda wildcards: define_input_for_pca(wildcards, "colors"),
        bs = config['pca_bs']
    log:
        temp(return_log_mc("{analysis_name}__{ref_genome}", "plot_PCA_correlation", "{env}"))
    conda: CONDA_ENV
    threads: config["resources"]["plot_PCA_correlation"]["threads"]
    resources:
        mem_mb=config["resources"]["plot_PCA_correlation"]["mem_mb"],
        tmp_mb=config["resources"]["plot_PCA_correlation"]["tmp_mb"],
        qos=config["resources"]["plot_PCA_correlation"]["qos"]
    shell:
        """
        {{
        printf "Plotting PCA for {wildcards.analysis_name} {wildcards.ref_genome} for {wildcards.env} samples\n"
        plotPCA -in {input.array} -T "PCA for {wildcards.env} in {params.bs}bp bins" -o {output.plot} --colors {params.colors:q} --transpose
        }} 2>&1 | tee -a "{log}"
        """

###
# final rule
rule all_combined:
    input:
        stats = define_final_stats_output(),
        final = lambda wildcards: define_final_combined_output(wildcards.ref_genome)
    output:
        touch = "results/combined/chkpts/combined_analysis__{analysis_name}__{ref_genome}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """        
