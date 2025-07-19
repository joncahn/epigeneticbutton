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
    elif target_name in ["all_genes","protein_coding_genes"]:
        file = f"results/combined/tracks/{ref_genome}__{target_name}.bed"
    else:
        raise ValueError(   
            f"{target_name} does not match possible files." 
            "It can be 'combined_peaks', or the value of "
            "'combined_target_file_name' in the config file"
        )
    
    return file

def has_header(bedfile):
    with open(bedfile) as f:
        first_line = f.readline().strip().split('\t')
        try:
            return not (int(first_line[1]) >=0 and int(first_line[2]) >=0)
        except (ValueError, IndexError):
            return true
           
def is_stranded(bedfile):
    strand_values = set()
    with open(bedfile) as f:
        if has_header(bedfile):
            next(f)
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 6:
                return False
            else:
                strand_values.add(cols[5])
    if strand_values.issubset({"+","-"}):
        return True
    else:
        return False

def define_matrix_per_target_name(wildcards):
    tname = config['combined_target_file_label']
    stranded_heatmaps = config['heatmaps']['stranded_heatmaps']
    matrix_param = wildcards.matrix_param
    env = wildcards.env
    analysis_name=config['analysis_name']
    ref_genome = wildcards.ref_genome
    target_name = wildcards.target_name
    prefix = f"results/combined/matrix/matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}"
    file = define_combined_target_file(wildcards)
    
    if stranded_heatmaps and is_stranded(file):
        return [f"{prefix}__plus.gz", f"{prefix}__minus.gz"]
    else:
        return [f"{prefix}__unstranded.gz"]

def define_sort_options(wildcards):
    sort_options = config['heatmaps']['sort_options']
    matrix_param = wildcards.matrix_param
    env = wildcards.env
    analysis_name=config['analysis_name']
    ref_genome = wildcards.ref_genome
    if sort_options == "no":
        return "--sortRegions keep"
    elif env == "mC":
        if target_name.endswith("sorted_regions"):
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

def define_labels_per_env_and_ref(wildcards):
    labels = []
    srna_sizes = config['srna_heatmap_sizes']
    strand = wildcards.strand
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    if globenv == "all":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] != "mC") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:    
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        if row.env == "TF":
            label=f"{row.line}_{row.tissue}_{row.extra_info}"
            labels.append(label)
        elif row.env == "ChIP":
            label=f"{row.line}_{row.tissue}_{row.sample_type}"
            labels.append(label)
        elif row.env == "mC":
            for context in ["CG","CHG","CHH"]:
                label=f"{row.line}_{row.tissue}_{context}"
                labels.append(label)
        elif row.env == "RNA":
            if strand == "unstranded":
                label1=f"{row.line}_{row.tissue}_{row.sample_type}_plus"
                label2=f"{row.line}_{row.tissue}_{row.sample_type}_minus"
                labels.append(label1)
                labels.append(label2)
            else:
                label=f"{row.line}_{row.tissue}_{row.sample_type}"
                labels.append(label)
        elif row.env == "sRNA":
            for size in srna_sizes:
                if strand == "unstranded":
                    label1=f"{row.line}_{row.tissue}_{row.sample_type}_{size}nt_plus"
                    label2=f"{row.line}_{row.tissue}_{row.sample_type}_{size}nt_minus"
                    labels.append(label1)
                    labels.append(label2)
                else:
                    label=f"{row.line}_{row.tissue}_{row.sample_type}_{size}nt"
                    labels.append(label)
            
    return labels

def define_bigwigs_per_env_and_ref(wildcards):
    bigwigs = []
    srna_sizes = config['srna_heatmap_sizes']
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    strand = wildcards.strand
    if globenv == "all":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] != "mC") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        if row.env in ["TF","ChIP"]:
            merged = f"FC__merged__{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}.bw"
            reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
            onerep = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{reps[0]}__{row.ref_genome}"
            bw = f"results/{row.env}/tracks/{merged}" if len(reps) >=2 else f"results/{row.env}/tracks/FC__final__{onerep}.bw"
            bigwigs.append(bw)
        elif row.env == "RNA":
            if strand == "unstranded":
                merged = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}"
                reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
                onerep = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{reps[0]}__{row.ref_genome}"
                bw1 = f"results/{row.env}/tracks/{merged}__plus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__plus.bw"
                bw2 = f"results/{row.env}/tracks/{merged}__minus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__minus.bw"
                bigwigs.append(bw1)
                bigwigs.append(bw2)
            else:
                merged = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}"
                reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
                onerep = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{reps[0]}__{row.ref_genome}"
                bw = f"results/{row.env}/tracks/{merged}__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{strand}.bw"
                bigwigs.append(bw)
        elif row.env == "sRNA":
            for size in srna_sizes:
                if strand == "unstranded":
                    merged = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}"
                    reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
                    onerep = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{reps[0]}__{row.ref_genome}"
                    bw1 = f"results/{row.env}/tracks/{merged}__{size}nt__plus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__plus.bw"
                    bw2 = f"results/{row.env}/tracks/{merged}__{size}nt__minus.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__minus.bw"
                    bigwigs.append(bw1)
                    bigwigs.append(bw2)
                else:
                    merged = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}"
                    reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
                    onerep = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{reps[0]}__{row.ref_genome}"
                    bw = f"results/{row.env}/tracks/{merged}__{size}nt__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{size}nt__{strand}.bw"
                    bigwigs.append(bw)
        elif row.env == "mC":
            merged = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__merged__{row.ref_genome}"
            reps = analysis_to_replicates.get((row.data_type, row.line, row.tissue, row.sample_type, row.ref_genome), [])
            onerep = f"{row.data_type}__{row.line}__{row.tissue}__{row.sample_type}__{reps[0]}__{row.ref_genome}"
            for context in ["CG","CHG","CHH"]:
                if strand == "unstranded":
                    bw = f"results/{row.env}/tracks/{merged}__{context}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{context}.bw"
                    bigwigs.append(bw)
                else:            
                    bw = f"results/{row.env}/tracks/{merged}__{context}__{strand}.bw" if len(reps) >=2 else f"results/{row.env}/tracks/{onerep}__{context}__{strand}.bw"
                    bigwigs.append(bw2)
    
    return bigwigs

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
    
    plot_files.append(f"results/combined/plots/Heatmap__regions__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    plot_files.append(f"results/combined/plots/Heatmap__tss__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    plot_files.append(f"results/combined/plots/Heatmap__tes__all__{analysis_name}__{ref_genome}__all_genes.pdf")
    if len(mc_analysis_samples) >=1: 
        if len(all_analysis_samples) > len(mc_analysis_samples):
            plot_files.append(f"results/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_genes_sorted_regions.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_genes_sorted_regions.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_genes_sorted_regions.pdf")
        else:
            plot_files.append(f"results/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"results/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
    
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
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes"
    output:
        temp_bedfile = temp("results/combined/bedfiles/temp__{target_name}__{ref_genome}.bed"),
        annotated_file = "results/combined/bedfiles/annotated__{target_name}__{ref_genome}.bed"
    params:
        target_name = lambda wildcards: wildcards.target_name
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
        # checking for presence of header
        read -r chrom start end _ < {input.bedfile}
        if [[ "${{start}}" =~ ^[0-9]+$ ]] && [[ "${{end}}" =~ ^[0-9]+$ ]]; then
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
        bigwigs = lambda wildcards: define_bigwigs_per_env_and_ref(wildcards),
        target_file = lambda wildcards: define_combined_target_file(wildcards)
    output:
        matrix = temp("results/combined/matrix/matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}__{strand}.gz"),
        temp = temp("results/combined/matrix/temp_region_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}__{strand}.bed")
    wildcard_constraints:
        strand = "plus|minus"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        labels = lambda wildcards: define_labels_per_env_and_ref(wildcards),
        matrix = lambda wildcards: wildcards.matrix_param,
        strand = lambda wildcards: wildcards.strand,
        header = lambda wildcards: "yes" if has_header(define_combined_target_file(wildcards)) else "no",
        params = lambda wildcards: config['heatmaps'][wildcards.matrix_param]['base'],
        bs = lambda wildcards: config['heatmaps'][wildcards.matrix_param]['bs'],
        before = lambda wildcards: config['heatmaps'][wildcards.matrix_param]['before'],
        after = lambda wildcards: config['heatmaps'][wildcards.matrix_param]['after'],
        middle = lambda wildcards: config['heatmaps'][wildcards.matrix_param]['middle']
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
        if [[ "{params.strand}" == "unstranded" ]]; then
            if [[ "{params.header}" == "no" ]]; then
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
        printf "Making {params.strand} strand {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
        computeMatrix {params.params} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.labels} -bs {params.bs} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        }} 2>&1 | tee -a "{log}"
        """
        
rule merging_matrix:
    input:
        matrix = lambda wildcards: define_matrix_per_target_name(wildcards)
    output:
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.gz"
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
        count=$(echo {input.matrix} | wc -w)
        if [[ ${{count}} -eq 2 ]]; then
            printf "\nMerging stranded matrices aligned by {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
			computeMatrixOperations rbind -m {input.matrix[0]} {input.matrix[1]} -o {output.matrix}
        else
            cp {input.matrix} {output.matrix}
        fi
        }} 2>&1 | tee -a "{log}"
        """
        
rule computing_matrix_scales:
    input:
        matrix = lambda wildcards: define_matrix_per_target_name(wildcards),
        target_file = lambda wildcards: define_combined_target_file(wildcards)
    output:
        params = "results/combined/matrix/params_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt",
        temp_values = temp("results/combined/matrix/temp_values_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"),
        temp_profile = temp("results/combined/matrix/temp_profile_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf"),
        temp_profile_values = temp("results/combined/matrix/temp_profile_values_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        scales = config['heatmaps']['scales'],
        header = lambda wildcards: "yes" if has_header(define_combined_target_file(wildcards)) else "no"
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
        count=$(wc -l {input.target_file} | cut -d' ' -f 1)
        if [[ "{params.header}" == "yes" ]]; then
            count=$((count-1))
        fi

        if [[ "{params.scales}" == "default" ]]; then
            awk -v ORS="" -v r=${{count}} -v n={params.target_name} 'BEGIN {{print "--regionsLabel "n"("r")"}}' > {output.params}
            touch {output.temp_values}
            touch {output.temp_profile}
            touch {output.temp_profile_values}
        elif [[ "{params.scales}" == "type" ]]; then
            printf "Getting scales for {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations dataRange -m {input.matrix} > {output.temp_values}
            plotProfile -m {input.matrix} -out {output.temp_profile} --averageType mean --outFileNameData {output.temp_profile}
            
            awk -v ORS="" -v r=${{count}} -v n={params.target_name} 'BEGIN {{print "--regionsLabel "n"("r")"}}' > {output.params}
            
        elif [[ "{params.scales}" == "sample" ]]; then
            
            printf "--regionsLabel ${{regionlabel}}" > {output.params}
            computeMatrixOperations dataRange -m {input.matrix} > {output.temp_values}
            plotProfile -m {input.matrix} -out {output.temp_profile} --averageType mean --outFileNameData {output.temp_profile}
            
            awk -v ORS="" -v r=${{count}} -v n={params.target_name} 'BEGIN {{print "--regionsLabel "n"("r")"}}' > {output.params}
        else
            printf "{params.scales} unknown. Returning default\n"
            awk -v ORS="" -v r=${{count}} -v n={params.target_name} 'BEGIN {{print "--regionsLabel "n"("r")"}}' > {output.params}
            touch {output.temp_values}
            touch {output.temp_profile}
            touch {output.temp_profile_values}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_heatmap_on_targetfile:
    input:
        matrix = "results/combined/matrix/final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.gz",
        params = "results/combined/matrix/params_final_matrix_{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.txt"
    output:
        plot = "results/combined/plots/Heatmap__{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}.pdf",
        sorted_regions = "results/combined/matrix/Heatmap__{matrix_param}__{env}__{analysis_name}__{ref_genome}__{target_name}_sorted_regions.bed"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        env = lambda wildcards: wildcards.env,
        plot_params = lambda wildcards: config['heatmaps']['plot_params'][wildcards.env],
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
        new_params="$(cat {input.params})"
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
