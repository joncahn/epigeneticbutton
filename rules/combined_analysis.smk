# function to access logs more easily
def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER,"combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_input_bedfile(bedname):
    if bedname 

def define_refgenome_bedfile(bedname):
    if bedname 

rule get_annotations_for_bedfile:
    input:
        bedfile = lambda wildcard: define_input_bedfile(wildcards.bedname),
        gff_file = lambda wildcard: f"genomes/{define_refgenome_bedfile(wildcards.bedname)}/{define_refgenome_bedfile(wildcards.bedname)}.gff",
        chrom_sizes = lambda wildcard: f"genomes/{define_refgenome_bedfile(wildcards.bedname)}/chrom.sizes"
    output:
        annotated_file = "combined/bedfiles/annotated__{bedname}.bed"




rule prep_files_for_differential_enrichment:
    input: 
        lambda wildcards: define_input_for_differential_enrichment(wildcards.env, wildcards.ref_genome)
    output:
        rna_samples = "RNA/DEG/samples__{analysis_name}__{ref_genome}.txt",
        rna_counts = "RNA/DEG/counts__{analysis_name}__{ref_genome}.txt"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_rna("{ref_genome}", "prep_for_DEGs", "{analysis_name}"))
    threads: config["resources"]["prep_files_for_DEGs"]["threads"]
    resources:
        mem=config["resources"]["prep_files_for_DEGs"]["mem"],
        tmp=config["resources"]["prep_files_for_DEGs"]["tmp"]
    run:
        filtered_samples = samples[ (samples['data_type'] == 'RNAseq') & (samples['ref_genome'] == params.ref_genome) ].copy()
        filtered_samples['Sample'] = filtered_samples['line'] + "__" + filtered_samples['tissue']
        filtered_samples['Replicate'] = filtered_samples['Sample'] + "__" + filtered_samples['replicate'].astype(str)
        
        RNA_samples = filtered_samples[['Replicate','Sample']].drop_duplicates()    
        RNA_samples = RNA_samples.sort_values(by=['Sample', 'Replicate'],ascending=[True, True]).reset_index(drop=True)
        RNA_samples['Color'] = pd.factorize(RNA_samples['Sample'])[0] + 1

        RNA_samples.to_csv(output.rna_samples, sep="\t", index=False)
        
        RNA_counts = None
        replicates = filtered_samples[['sample_name', 'Replicate']].drop_duplicates()
        for sname, rep in replicates.values:
            file_path = f"RNA/DEG/counts__{sname}.tab"
            temp = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1])
            temp.columns = ['GID', rep]

            if RNA_counts is None:
                RNA_counts = temp
            else:
                RNA_counts = pd.merge(RNA_counts, temp, on='GID', how='outer')
            
        replicate_order = RNA_samples['Replicate'].tolist()
        column_order = ['GID'] + replicate_order
        RNA_counts = RNA_counts[column_order]
        RNA_counts.to_csv(output.rna_counts, sep="\t", index=False)