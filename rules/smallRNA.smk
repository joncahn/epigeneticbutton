# function to access logs more easily
def return_log_smallrna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER,"smallRNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def define_input_file_for_shortstack(sample_name):
    paired = get_sample_info_from_name(sample_name, samples, 'paired')
    if paired == "se":
        return "cleaned_{sample_name}__R0" if config['structural_rna_depletion'] else "trim_{sample_name}__R0"
    
rule shortstack_map:
    input:
        fastq = lambda wildcards: f"sRNA/fastq/{define_input_file_for_shortstack(wildcards.sample_name)}.fastq.gz",
        fasta = lambda wildcards: f"genomes/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa"
    output:
        count_file = "RNA/mapped/{sample_name}/ShortStack_All.gff3",
        bam_file = temp("RNA/mapped/star_se__{sample_name}_ReadsPerGene.out.tab"),
        metrics_map = "RNA/reports/star_se__{sample_name}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        srna_params = config['srna_mapping_params']
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "SE"))
    conda: CONDA_ENV
    threads: config["resources"]["shortstack_map"]["threads"]
    resources:
        mem=config["resources"]["shortstack_map"]["mem"],
        tmp=config["resources"]["shortstack_map"]["tmp"]
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with Shortstack version:\n"
        ShortStack --version
        ShortStack --readfile {input.fastq} --genomefile {input.fasta} --bowtie_cores $threads --sort_mem {resources.mem} {params.srna_params} --outdir sRNA/mapped/{params.sample_name}
        }} 2>&1 | tee -a "{log}"
        """