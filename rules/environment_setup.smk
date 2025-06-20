# CONDA_ENV=os.path.join(REPO_FOLDER,"envs/reference.yaml")

# function to access logs more easily
def return_log_env(ref_genome, step):
    return os.path.join(REPO_FOLDER,"logs",f"tmp_{step}_{ref_genome}.log")

# Function to create directories
def create_directories(unique_envs, dirs):
    for env in unique_envs:
        for d in ["fastq", "mapped", "tracks", "reports", "logs", "chkpts", "plots"]:
            os.makedirs(f"{env}/{d}", exist_ok=True)
        if env in ["ChIP", "TF"]:
            os.makedirs(f"{env}/peaks", exist_ok=True)
        if env in ["mC"]:
            os.makedirs(f"{env}/methylcall", exist_ok=True)
            os.makedirs(f"{env}/DMRs", exist_ok=True)
        if env in ["RNA"]:
            os.makedirs(f"{env}/DEG", exist_ok=True)
    
    for key, value in dirs.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                os.makedirs(sub_value, exist_ok=True)
        else:
            os.makedirs(value, exist_ok=True)

# Rule to summarize the preparation of the reference genome
rule prepare_reference:
    input:
        setup = "chkpts/directories_setup.done",
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gff = "genomes/{ref_genome}/{ref_genome}.gff",
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        region_files = ["combined/tracks/{ref_genome}__protein_coding_genes.bed", "combined/tracks/{ref_genome}__all_genes.bed"],
        logs = lambda wildcards: [ return_log_env(wildcards.ref_genome, step) for step in ["fasta", "gff", "gtf", "chrom_sizes", "region_file"] ]
    output:
        chkpt = "chkpts/ref__{ref_genome}.done",
        log = os.path.join(REPO_FOLDER,"logs","ref_prep__{ref_genome}.log")
    threads: 1
    resources:
        mem=1,
        tmp=1
    shell:
        """
        cat {input.logs} > {output.log}
        rm {input.logs}
        touch {output.chkpt}
        """

# Call the function to create directories
rule setup_directories:
    output:
        touch = "chkpts/directories_setup.done"
    run:
        create_directories(UNIQUE_ENVS, DIRS)
        with open(output.touch, "w") as f:
            f.write("Setup complete\n")

# Rule to make sure a fasta file is found, and unzipped it if needed
rule check_fasta:
    output:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref_genome)
    log:
        temp(return_log_env("{ref_genome}", "fasta"))
    conda: CONDA_ENV
    threads: config["resources"]["use_pigz"]["threads"]
    resources:
        mem=config["resources"]["use_pigz"]["mem"],
        tmp=config["resources"]["use_pigz"]["tmp"]
    shell:
        """
        # Search for fasta file
        if [ -s {params.ref_dir}/*.fa.gz ]; then
            fa_file=$(ls {params.ref_dir}/*.fa.gz)
            fa_filename=${{fa_file##*/}}
            printf "\nGzipped fasta file found in {params.ref_dir}:\n ${{fa_filename}}\n" >> {log} 2>&1
            pigz -p {threads} -dc ${{fa_file}} > {output.fasta}
        elif [ -s {params.ref_dir}/*.fa ]; then
            fa_file=$(ls {params.ref_dir}/*.fa)
            fa_filename=${{fa_file##*/}}
            printf "\nUnzipped fasta file found in {params.ref_dir}:\n ${{fa_filename}}\n" >> {log} 2>&1
            cp ${{fa_file}} {output.fasta}
        elif [ -s {params.ref_dir}/*.fasta.gz ]; then
            fa_file=$(ls {params.ref_dir}/*.fasta.gz)
            fa_filename=${{fa_file##*/}}
            printf "\nGzipped fasta file found in {params.ref_dir}:\n ${{fa_filename}}\n" >> {log} 2>&1
            pigz -p {threads} -dc ${{fa_file}} > {output.fasta}
        elif [ -s {params.ref_dir}/*.fasta ]; then
            fa_file=$(ls {params.ref_dir}/*.fasta)
            fa_filename=${{fa_file##*/}}
            printf "\nUnzipped fasta file found in {params.ref_dir}:\n ${{fa_filename}}\n" >> {log} 2>&1
            cp ${{fa_file}} {output.fasta}
        else
            printf "\nNo fasta file found in reference directory:\n {params.ref_dir}\n" >> {log} 2>&1
            exit 1
        fi
        """
        
rule check_gff:
    output:
        gff = "genomes/{ref_genome}/{ref_genome}.gff"
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref_genome)
    log:
        temp(return_log_env("{ref_genome}", "gff"))
    conda: CONDA_ENV
    threads: config["resources"]["use_pigz"]["threads"]
    resources:
        mem=config["resources"]["use_pigz"]["mem"],
        tmp=config["resources"]["use_pigz"]["tmp"]
    shell:
        """
        if [ -s {params.ref_dir}/*.gff*.gz ]; then
            gff_file=$(ls {params.ref_dir}/*gff*.gz)
            gff_filename=${{gff_file##*/}}
            printf "\nGzipped GFF annotation file found in {params.ref_dir}:\n ${{gff_filename}}\n" >> {log} 2>&1
            pigz -p {threads} -dc ${{gff_file}} > {output.gff}	
        elif [ -s {params.ref_dir}/*.gff* ]; then
            gff_file=$(ls {params.ref_dir}/*.gff*)
            gff_filename=${{gff_file##*/}}
            printf "\nUnzipped GFF annotation file found in {params.ref_dir}:\n ${{gff_filename}}\n" >> {log} 2>&1
            cp ${{gff_file}} {output.gff}
        else
            printf "\nNo gff annotation file found in reference directory:\n {params.ref_dir}\n" >> {log} 2>&1
            exit 1
        fi
        """

rule check_gtf:
    output:
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf"
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref_genome)
    log:
        temp(return_log_env("{ref_genome}", "gtf"))
    conda: CONDA_ENV
    threads: config["resources"]["use_pigz"]["threads"]
    resources:
        mem=config["resources"]["use_pigz"]["mem"],
        tmp=config["resources"]["use_pigz"]["tmp"]
    shell:
        """
        if [ -s {params.ref_dir}/*.gtf.gz ]; then
            gtf_file=$(ls {params.ref_dir}/*gtf.gz)
            gtf_filename=${{gtf_file##*/}}
            printf "\nGzipped GTF annotation file found in {params.ref_dir}:\n ${{gtf_filename}}\n" >> {log} 2>&1
            sed 's/gene://' ${{gtf_file}} | pigz -p {threads} -dc > {output.gtf}	
        elif [ -s {params.ref_dir}/*.gtf ]; then
            gtf_file=$(ls {params.ref_dir}/*.gtf)
            gtf_filename=${{gtf_file##*/}}
            printf "\nUnzipped GTF annotation file found in {params.ref_dir}:\n ${{gtf_filename}}\n" >> {log} 2>&1
            sed 's/gene://' ${{gtf_file}} > {output.gtf}
        else
            printf "\nNo GTF annotation file found in reference directory:\n {params.ref_dir}\n" >> {log} 2>&1
            exit 1
        fi
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
    threads: config["resources"]["chrom_sizes"]["threads"]
    resources:
        mem=config["resources"]["chrom_sizes"]["mem"],
        tmp=config["resources"]["chrom_sizes"]["tmp"]
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
        region_file1 = "combined/tracks/{ref_genome}__protein_coding_genes.bed",
        region_file2 = "combined/tracks/{ref_genome}__all_genes.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "region_file"))
    conda: CONDA_ENV
    threads: config["resources"]["region_file"]["threads"]
    resources:
        mem=config["resources"]["region_file"]["mem"],
        tmp=config["resources"]["region_file"]["tmp"]
    shell:
        """
        printf "\nMaking a bed file with gene coordinates from {params.ref_genome}\n" >> {log} 2>&1
        awk -v OFS="\t" '$3=="gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        """
        
