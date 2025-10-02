# function to access logs more easily
def return_log_env(ref_genome, step):
    return os.path.join(REPO_FOLDER,"results","logs",f"tmp_{step}_{ref_genome}.log")

# Rule to summarize the preparation of the reference genome
rule prepare_reference:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gff = "genomes/{ref_genome}/{ref_genome}.gff",
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        region_files = ["results/combined/tracks/{ref_genome}__protein_coding_genes.bed", "results/combined/tracks/{ref_genome}__all_genes.bed"],
        logs = lambda wildcards: [ return_log_env(wildcards.ref_genome, step) for step in ["fasta", "gff", "gtf", "chrom_sizes", "region_file"] ]
    output:
        chkpt = "results/combined/chkpts/ref__{ref_genome}.done",
        log = os.path.join(REPO_FOLDER,"results","logs","ref_prep__{ref_genome}.log")
    localrule: True
    shell:
        """
        cat {input.logs} > {output.log}
        rm {input.logs}
        touch {output.chkpt}
        """

# Rule to make sure a fasta file is found, and unzipped it if needed
rule check_fasta:
    output:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    params:
        fasta = lambda wildcards: config[wildcards.ref_genome]['fasta_file']
    log:
        temp(return_log_env("{ref_genome}", "fasta"))
    conda: CONDA_ENV
    threads: config["resources"]["check_fasta"]["threads"]
    resources:
        mem=config["resources"]["check_fasta"]["mem"],
        tmp=config["resources"]["check_fasta"]["tmp"]
    shell:
        """
        {{
        # Looks if fasta is present and gzipped or not
        if [[ {params.fasta} == *.fa.gz || {params.fasta} == *.fasta.gz ]]; then
            printf "\nGzipped fasta file found: {params.fasta}\n"
            pigz -p {threads} -dc {params.fasta} > {output.fasta}
        elif [[ {params.fasta} == *.fa || {params.fasta} == *.fasta ]]; then
            printf "\nUnzipped fasta file found: {params.fasta}\n"
            cp {params.fasta} {output.fasta}
        else
            printf "\nNo fasta file (with .fasta(.gz) or .fa(.gz) extension) found at the location:\n {params.fasta}\n"
            exit 1
        fi
        }} >> 2>&1 | tee -a "{log}"
        """
        
rule check_gff:
    output:
        gff = "genomes/{ref_genome}/{ref_genome}.gff"
    params:
        gff = lambda wildcards: config[wildcards.ref_genome]['gff_file']
    log:
        temp(return_log_env("{ref_genome}", "gff"))
    conda: CONDA_ENV
    threads: config["resources"]["check_gff"]["threads"]
    resources:
        mem=config["resources"]["check_gff"]["mem"],
        tmp=config["resources"]["check_gff"]["tmp"]
    shell:
        """
        {{
        # Looks if gff is present and is gzipped or not
        if [[ {params.gff} == *.gff*.gz ]]; then
            printf "\nGzipped gff file found: {params.gff}\n"
            pigz -p {threads} -dc {params.gff} > {output.gff}
        elif [[ {params.gff} == *.gff* ]]; then
            printf "\nUnzipped gff file found: {params.gff}\n"
            cp {params.gff} {output.gff}
        else
            printf "\nNo gff file (with .gff*(.gz) extension) found at the location:\n {params.gff}\n"
            exit 1
        fi
        }} >> 2>&1 | tee -a "{log}"
        """

rule check_gtf:
    output:
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf"
    params:
        gtf = lambda wildcards: config[wildcards.ref_genome]['gtf_file']
    log:
        temp(return_log_env("{ref_genome}", "gtf"))
    conda: CONDA_ENV
    threads: config["resources"]["check_gtf"]["threads"]
    resources:
        mem=config["resources"]["check_gtf"]["mem"],
        tmp=config["resources"]["check_gtf"]["tmp"]
    shell:
        """
        {{
        # Looks if gtf is present and is gzipped or not
        if [[ {params.gtf} == *.gtf.gz ]]; then
            printf "\nGzipped gtf file found: {params.gtf}\n"
            pigz -p {threads} -dc {params.gtf} > {output.gtf}
        elif [[ {params.gtf} == *.gtf ]]; then
            printf "\nUnzipped gtf file found: {params.gtf}\n"
            cp {params.gtf} {output.gtf}
        else
            printf "\nNo gtf file (with .gtf(.gz) extension) found at the location:\n {params.gtf}\n"
            exit 1
        fi
        }} >> 2>&1 | tee -a "{log}"
        """
        
rule check_chrom_sizes:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    output:
        fasta_index = "genomes/{ref_genome}/{ref_genome}.fa.fai",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "chrom_sizes"))
    conda: CONDA_ENV
    threads: config["resources"]["check_chrom_sizes"]["threads"]
    resources:
        mem=config["resources"]["check_chrom_sizes"]["mem"],
        tmp=config["resources"]["check_chrom_sizes"]["tmp"]
    shell:
        """
        printf "\nMaking chrom.sizes file for {params.ref_genome}\n" >> {log} 2>&1
        samtools faidx {input.fasta}
        cut -f1,2 {output.fasta_index} > {output.chrom_sizes}
        """

rule prep_region_file:
    input:
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        gff = "genomes/{ref_genome}/{ref_genome}.gff"
    output:
        region_file1 = "results/combined/tracks/{ref_genome}__protein_coding_genes.bed",
        region_file2 = "results/combined/tracks/{ref_genome}__all_genes.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "region_file"))
    conda: CONDA_ENV
    threads: config["resources"]["prep_region_file"]["threads"]
    resources:
        mem=config["resources"]["prep_region_file"]["mem"],
        tmp=config["resources"]["prep_region_file"]["tmp"]
    shell:
        """
        printf "\nMaking a bed file with gene coordinates from {params.ref_genome}\n" >> {log} 2>&1
        awk -v OFS="\t" '$3=="gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        """
        
