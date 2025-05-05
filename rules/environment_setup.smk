# function to access logs more easily
def return_log(ref, step):
    return os.path.join(REPO_FOLDER,"logs",f"tmp_{step}_{ref}.log")

CONDA_ENV=os.path.join(REPO_FOLDER,"envs/reference.yaml")

# Rule to summarize the preparation of the reference genome
rule prepare_reference:
    input:
        fasta = "genomes/{ref}/temp_{ref}.fa",
        gff = "genomes/{ref}/temp_{ref}.gff",
        gtf = "genomes/{ref}/temp_{ref}.gtf",
        chrom_sizes = "genomes/{ref}/chrom.sizes",
        region_files = ["combined/tracks/{ref}_protein_coding_genes.bed", "combined/tracks/{ref}_all_genes.bed"],
        logs = lambda wildcards: [ return_log(wildcards.ref, step) for step in ["fasta", "gff", "gtf", "chrom_sizes", "region_file"] ]
    output:
        chkpt = "chkpts/ref__{ref}.done",
        log = os.path.join(REPO_FOLDER,"logs","ref_prep__{ref}.log")
    shell:
        """
        cat {input.logs} > {output.log}
        rm {input.logs}
        touch {output.chkpt}
        """

# Rule to make sure a fasta file is found, and unzipped it if needed
rule check_fasta:
    output:
        fasta = "genomes/{ref}/temp_{ref}.fa"
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref)
    log:
        return_log("{ref}", "fasta")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        # Search for fasta file
        if [ -s {params.ref_dir}/*.fa.gz ]; then
            fa_file=$(ls ${params.ref_dir}/*.fa.gz)
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
        gff = "genomes/{ref}/temp_{ref}.gff"
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref)
    log:
        return_log("{ref}", "gff")
    conda:
        CONDA_ENV
    threads: workflow.cores
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
        gtf = "genomes/{ref}/temp_{ref}.gtf"
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref)
    log:
        return_log("{ref}", "gtf")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        if [ -s {params.ref_dir}/*.gtf.gz ]; then
            gtf_file=$(ls {params.ref_dir}/*gtf.gz)
            gtf_filename=${{gtf_file##*/}}
            printf "\nGzipped GTF annotation file found in {params.ref_dir}:\n ${{gtf_filename}}\n" >> {log} 2>&1
            pigz -p {threads} -dc ${{gtf_file}} > {output.gtf}	
        elif [ -s {params.ref_dir}/*.gtf ]; then
            gtf_file=$(ls {params.ref_dir}/*.gtf)
            gtf_filename=${{gtf_file##*/}}
            printf "\nUnzipped GTF annotation file found in {params.ref_dir}:\n ${{gtf_filename}}\n" >> {log} 2>&1
            cp ${{gtf_file}} {output.gtf}
        else
            printf "\nNo GTF annotation file found in reference directory:\n {params.ref_dir}\n" >> {log} 2>&1
            exit 1
        fi
        """
        
rule check_chrom_sizes:
    input:
        fasta = "genomes/{ref}/temp_{ref}.fa"
    output:
        fasta_index = "genomes/{ref}/temp_{ref}.fa.fai",
        chrom_sizes = "genomes/{ref}/chrom.sizes"
    log:
        return_log("{ref}", "chrom_sizes")
    conda:
        CONDA_ENV
    shell:
        """
        printf "\nMaking chrom.sizes file for {ref}\n" >> {log} 2>&1
        samtools faidx {input.fasta}
        cut -f1,2 {output.fasta_index} > {output.chrom_sizes}
        """

rule prep_region_file:
    input:
        chrom_sizes = "genomes/{ref}/chrom.sizes",
        gff = "genomes/{ref}/temp_{ref}.gtf"
    output:
        region_file1 = "combined/tracks/{ref}_protein_coding_genes.bed",
        region_file2 = "combined/tracks/{ref}_all_genes.bed"
    log:
        return_log("{ref}", "region_file")
    conda:
        CONDA_ENV
    shell:
        """
        printf "\nMaking a bed file with gene coordinates from {ref}\n" >> {log} 2>&1
        awk -v OFS="\t" '$3=="gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        """
        
