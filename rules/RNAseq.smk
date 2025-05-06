# function to access logs more easily
def return_log(ref, step):
    return os.path.join(REPO_FOLDER,"logs",f"tmp_chip_{step}_{ref}.log")

CONDA_ENV=os.path.join(REPO_FOLDER,"envs/reference.yaml")

rule stat_file:
    output:
        stat_file = "RNA/reports/summary_mapping_stats.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """

rule make_RNA_indices:
    input:
        fasta = "genomes/{ref}/temp_{ref}.fa",
        gtf = "genomes/{ref}/temp_{ref}.gtf"
    output:
        indices = "combined/genomes/{ref}/STAR_index"
    log:
        os.path.join(REPO_FOLDER,"logs",f"STAR_index_{ref}.log")
    threads: workflow.cores
    shell:
        """
        if [ ! -d {output.indices} ]; then
            printf "\nBuilding STAR index directory for ${ref}\n"
            mkdir {output.indices}
            STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.indices} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf}
        else
            printf "\nSTAR index already exists for {ref}\n"
        fi
        """
