# function to access logs more easily
def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER,"results","combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_samplenames_per_env(wildcards, env):
    names = []
    ref_genome = wildcards.ref_genome
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == env) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        names.append(sample_name_str(row, 'analysis'))
    
    return names

def define_input_bedfile(wildcards):
    bedname = wildcards.bedname
    ref_genome = wildcards.ref_genome
    
    if bedname.startswith("combined_peakfiles"):
        file = f"results/combined/bedfiles/{bedname}__{ref_genome}.bed"
        
    return file

def define_final_combined_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    analysis_name = config['analysis_name']
    text_files = []
    plot_files = []
    
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
    
    if len(filtered_analysis_samples) >=2:
        text_files.append(f"results/combined/bedfiles/annotated__combined_peakfiles_{analysis_name}__{ref_genome}.bed")
        
    if analysis:
        results += plot_files + text_files
    
    return results

rule combined_peakfiles:
    input:
        chrom_sizes = lambda wildcards: f"genomes/{wildcards.ref_genome}/chrom.sizes",
        peakfiles = lambda wildcards: [ f"results/ChIP/peaks/selected_peaks__{names}.bedPeak" for names in define_peakfiles_for_combined(wildcards, "ChIP") ]
    output:
        temp1_file = "results/combined/bedfiles/temp1_combined_peakfiles_{analysis_name}__{ref_genome}.bed",
        temp2_file = "results/combined/bedfiles/temp2_combined_peakfiles_{analysis_name}__{ref_genome}.bed",
        merged_file = "results/combined/bedfiles/combined_peakfiles_{analysis_name}__{ref_genome}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        names = lambda wildcards: define_peakfiles_for_combined(wildcards, "ChIP"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_rna("{analysis_name}", "{ref_genome}", "annotate_bedfile"))
    threads: config["resources"]["combined_peakfiles"]["threads"]
    resources:
        mem=config["resources"]["combined_peakfiles"]["mem"],
        tmp=config["resources"]["combined_peakfiles"]["tmp"]
    shell:
        """
        {{
        for sample in {params.names}; do
            awk -v OFS="\t" -v s=${{sample}} '{{print $1,$2,$3,s}}' results/ChIP/peaks/selected_peaks__${{sample}}.bedPeak >> {output.temp_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" '{{print $1,$2,$3,"Peak_"NR,$4}}' > {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule get_annotations_for_bedfile:
    input:
        bedfile = lambda wildcard: define_input_bedfile(wildcards),
        region_file = lambda wildcard: f"results/combined/tracks/{wildcards.ref_genome}__all_genes.bed",
        chrom_sizes = lambda wildcard: f"genomes/{wildcards.ref_genome}/chrom.sizes"
    output:
        temp_bedfile = temp("results/combined/bedfiles/temp__{bedname}__{ref_genome}.bed")
        annotated_file = "results/combined/bedfiles/annotated__{bedname}__{ref_genome}.bed"
    params:
        bedname = lambda wildcards: wildcards.bedname
    log:
        temp(return_log_rna("{bedname}", "{ref_genome}", "annotate_bedfile"))
    threads: config["resources"]["get_annotations_for_bedfile"]["threads"]
    resources:
        mem=config["resources"]["get_annotations_for_bedfile"]["mem"],
        tmp=config["resources"]["get_annotations_for_bedfile"]["tmp"]
    shell:
        """
        {{
        awk -v OFS="\t" -v n={params.bedname} '{{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile)
        bedtools closest -a {output.temp_bedfile} -b {input.region_file} -g {input.chrom_sizes} -D ref | awk -v OFS="\t" '{{if ($11=="+") print $1,$2,$3,$4,$12,$11,$5,$9; else print $1,$2,$3,$4,-$12,$11,$5,$9}}' | awk -F"[=;]" -v OFS="\t" '{{print $1,$2}}' | sed 's/gene://' | awk -v OFS="\t" '{{print $1,$2,$3,$4,$5,$6,$7,$9}}' > {output.annotated_file}
        }} 2>&1 | tee -a "{log}"
        """

rule all_combined:
    input:
        final = lambda wildcards: define_final_combined_output(wildcards.ref_genome)
    output:
        touch = "results/combined/chkpts/combined_analysis__{analysis_name}__{ref_genome}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """        
