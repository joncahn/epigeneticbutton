# function to access logs more easily
def return_log(ref, step):
    return os.path.join(REPO_FOLDER,"logs",f"tmp_{step}_{ref}.log")

CONDA_ENV=os.path.join(REPO_FOLDER,"envs/reference.yaml")

# Rule to summarize the preparation of the reference genome
rule prepare_reference:
    input:
        fasta = os.path.join(REF_PATH, "{ref}", "temp_{ref}.fa"),
        gff = os.path.join(REF_PATH, "{ref}", "temp_{ref}.gff"),
        gtf = os.path.join(REF_PATH, "{ref}", "temp_{ref}.gtf"),
        chrom_sizes = os.path.join(REF_PATH, "{ref}","chrom.sizes"),
        region_files = ["combined/tracks/{ref}_protein_coding_genes.bed", "combined/tracks/{ref}_all_genes.bed"],
        logs = [ret_log("{ref}", "fasta"), ret_log("{ref}", "gff"), ret_log("{ref}", "gtf"), ret_log("{ref}", "chrom_sizes"), ret_log("{ref}", "region_file")]
    output:
        chkpt = "chkpts/ref__{ref}.done",
        log = os.path.join(REPO_FOLDER,"logs",f"ref_prep__{ref}.log")
    shell:
        """
        cat {input.logs} > {output.log}
        rm {input.logs}
        touch {output.chkpt}
        """

# Rule to make sure a fasta file is found, and unzipped it if needed
rule check_fasta:
    output:
        touch = "chkpts/genome_prep_{ref}.done",
        fasta = os.path.join(REF_PATH, "{ref}", "temp_{ref}.fa")
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref)
    log:
        return_log("{ref}", "fasta")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        # Search for fasta file
        ref_dir={params.ref_dir}
        if [ -s ${ref_dir}/*.fa.gz ]; then
            fa_file=$(ls ${ref_dir}/*.fa.gz)
            fa_filename=${fa_file##*/}
            printf "\nGzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            pigz -p {threads} -dc ${fa_file} > {output.fasta}
        elif [ -s ${ref_dir}/*.fa ]; then
            fa_file=$(ls ${ref_dir}/*.fa)
            fa_filename=${fa_file##*/}
            printf "\nUnzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            cp ${fa_file} {output.fasta}
        elif [ -s ${ref_dir}/*.fasta.gz ]; then
            fa_file=$(ls ${ref_dir}/*.fasta.gz)
            fa_filename=${fa_file##*/}
            printf "\nGzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            pigz -p {threads} -dc ${fa_file} > {output.fasta}
        elif [ -s ${ref_dir}/*.fasta ]; then
            fa_file=$(ls ${ref_dir}/*.fasta)
            fa_filename=${fa_file##*/}
            printf "\nUnzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            cp ${fa_file} {output.fasta}
        else
            printf "\nNo fasta file found in reference directory:\n ${ref_dir}\n"
            exit 1
        fi
        }} &> {log}
        """
        
rule check_gff:
    output:
        gff = os.path.join(REF_PATH, "{ref}", "temp_{ref}.gff")
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref)
    log:
        return_log("{ref}", "gff")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        """
        {{
        ref_dir={params.ref_dir}
        if [ -s ${ref_dir}/*.gff*.gz ]; then
            gff_file=$(ls ${ref_dir}/*gff*.gz)
            gff_filename=${gff_file##*/}
            printf "\nGzipped GFF annotation file found in ${ref_dir}:\n ${gff_filename}\n"
            pigz -p {threads} -dc ${gff_file} > {output.gff}	
        elif [ -s ${ref_dir}/*.gff* ]; then
            gff_file=$(ls ${ref_dir}/*.gff*)
            gff_filename=${gff_file##*/}
            printf "\nUnzipped GFF annotation file found in ${ref_dir}:\n ${gff_filename}\n"
            cp ${gff_file} {output.gff}
        else
            printf "\nNo gff annotation file found in reference directory:\n ${ref_dir}\n"
            exit 1
        fi
        }} &> {log}
        """

rule check_gtf:
    output:
        gff = os.path.join(REF_PATH, "{ref}", "temp_{ref}.gtf")
    params:
        ref_dir = lambda wildcards: os.path.join(REF_PATH, wildcards.ref)
    log:
        return_log("{ref}", "gtf")
    conda:
        CONDA_ENV
    threads: workflow.cores
    shell:
        r"""
        {{
        ref_dir={params.ref_dir}
        if [ -s ${ref_dir}/*.gtf.gz ]; then
            gtf_file=$(ls ${ref_dir}/*gtf.gz)
            gtf_filename=${gtf_file##*/}
            printf "\nGzipped GTF annotation file found in ${ref_dir}:\n ${gtf_filename}\n" 
            pigz -p {threads} -dc ${gtf_file} > {ouptut.gtf}	
        elif [ -s ${ref_dir}/*.gtf ]; then
            gtf_file=$(ls ${ref_dir}/*.gtf)
            gtf_filename=${gtf_file##*/}
            printf "\nUnzipped GTF annotation file found in ${ref_dir}:\n ${gtf_filename}\n"
            cp ${gtf_file} {ouptut.gtf}
        else
            printf "\nNo GTF annotation file found in reference directory:\n ${ref_dir}\n"
            exit 1
        fi
        }} &> {log}
        """
        
rule check_chrom_sizes:
    input:
        fasta = os.path.join(REF_PATH, "{ref}", "temp_{ref}.fa")
    output:
        fasta_index = os.path.join(REF_PATH, "{ref}", "temp_{ref}.fa.fai"),
        chrom_sizes = os.path.join(REF_PATH, "{ref}", "chrom.sizes")
    log:
        return_log("{ref}", "chrom_sizes")
    conda:
        CONDA_ENV
    shell:
        r"""
        {{
        printf "\nMaking chrom.sizes file for {ref}\n" 
        samtools faidx {input.fasta}
        cut -f1,2 {output.fasta_index} > {output.chrom.sizes}
        }} &> {log}
        """

rule prep_region_file:
    input:
        chrom_sizes = os.path.join(REF_PATH, "{ref}", "chrom.sizes"),
        gff = os.path.join(REF_PATH, "{ref}", "temp_{ref}.gtf")
    output:
        region_file1 = "combined/tracks/{ref}_protein_coding_genes.bed",
        region_file2 = "combined/tracks/{ref}_all_genes.bed"
    log:
        return_log("{ref}", "region_file")
    conda:
        CONDA_ENV
    shell:
        r"""
        {{
        printf "\nMaking a bed file with gene coordinates from {ref}\n"
        awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {print $1,$4-1,$5,$9,".",$7}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        }} &> {log}
        """
        
