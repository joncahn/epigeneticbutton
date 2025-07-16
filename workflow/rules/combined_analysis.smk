# function to access logs more easily
def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER,"results","combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_input_bedfile(bedname):
    if bedname 
    selected_peaks__{spname}.bed

def define_refgenome_bedfile(bedname):
    if bedname 
    selected_peaks__{spname}.bed

def define_final_combined_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    plot_files = []
    
    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
    
    if len(filtered_analysis_samples) >=2:
            
    if analysis:
        results += plot_files
    
    return results

rule get_annotations_for_bedfile:
    input:
        bedfile = lambda wildcard: define_input_bedfile(wildcards.bedname),
        region_file = lambda wildcard: f"results/combined/tracks/{define_refgenome_bedfile(wildcards.bedname)}__all_genes.bed",
        chrom_sizes = lambda wildcard: f"genomes/{define_refgenome_bedfile(wildcards.bedname)}/chrom.sizes"
    output:
        temp_bedfile = temp("results/combined/bedfiles/temp__{bedname}.bed")
        annotated_file = "results/combined/bedfiles/annotated__{bedname}.bed"
    params:
        bedname = lambda wildcards: wildcards.bedname
    log:
        temp(return_log_rna("{bedname}", "annotate_bedfile", ""))
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
