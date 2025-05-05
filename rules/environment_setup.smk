# Rule to prepare reference genome for each data type
rule prepare_reference:
    input:
        fasta = "{ref_dir}/temp_{data_type}_{ref}.fa"
        gff = "{ref_dir}/temp_{data_type}_{ref}.gff"
        gtf = "{ref_dir}/temp_{data_type}_{ref}.gtf"
        chrom_sizes = "{ref_dir}/chrom.sizes"
        region_files = ["${data_type}/tracks/${ref}_protein_coding_genes.bed", "${data_type}/tracks/${ref}_all_genes.bed"]
        stat_file = "${data_type}/reports/summary_mapping_stats.txt"
        logs = ["logs/tmp_fasta_{data_type}_{ref}.log", "logs/tmp_gff_{data_type}_{ref}.log", "logs/tmp_gtf_{data_type}_{ref}.log", "logs/tmp_chrom_size_{data_type}_{ref}.log", "logs/tmp_region_file_{data_type}_{ref}.log"]
     output:
        chkpt = "chkpts/ref__{ref_genome}__{env}.done"
        log = "logs/prepare_ref__{ref_genome}__{env}.log"
      shell:
        """
        cat {input.logs} > {output.log}
        rm {input.logs}
        touch {output.chkpt}
        """

# Rule to make sure a fasta file is found, and unzipped it if needed
rule check_fasta:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
    output:
        fasta = "{ref_dir}/temp_{data_type}_{ref}.fa"
    log:
        "logs/tmp_fasta_{data_type}_{ref}.log"
    conda:
        "envs/reference.yaml"
    threads: workflow.cores
    shell:
        """
        # Search for fasta file
        if [ -s {input.ref_dir}/*.fa.gz ]; then
            fa_file=$(ls {input.ref_dir}/*.fa.gz)
            fa_filename=${fa_file##*/}
            printf "\nGzipped fasta file found in {input.ref_dir}:\n ${fa_filename}\n"
            pigz -p {threads} -dc ${fa_file} > {output.fasta}
        elif [ -s {input.ref_dir}/*.fa ]; then
            fa_file=$(ls {input.ref_dir}/*.fa)
            fa_filename=${fa_file##*/}
            printf "\nUnzipped fasta file found in {input.ref_dir}:\n ${fa_filename}\n"
            cp ${fa_file} {output.fasta}
        elif [ -s {input.ref_dir}/*.fasta.gz ]; then
            fa_file=$(ls {input.ref_dir}/*.fasta.gz)
            fa_filename=${fa_file##*/}
            printf "\nGzipped fasta file found in {input.ref_dir}:\n ${fa_filename}\n"
            pigz -p {threads} -dc ${fa_file} > {output.fasta}
        elif [ -s {input.ref_dir}/*.fasta ]; then
            fa_file=$(ls {input.ref_dir}/*.fasta)
            fa_filename=${fa_file##*/}
            printf "\nUnzipped fasta file found in {input.ref_dir}:\n ${fa_filename}\n"
            cp ${fa_file} {output.fasta}
        else
            printf "\nNo fasta file found in reference directory:\n {input.ref_dir}\n"
            exit 1
        fi
        """
        
rule check_gff:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
    output:
        gff = "{ref_dir}/temp_{data_type}_{ref}.gff"
    log:
        "logs/tmp_gff_{data_type}_{ref}.log"
    conda:
        "envs/reference.yaml"
    threads: workflow.cores
    shell:
        """
        if [ -s {input.ref_dir}/*.gff*.gz ]; then
            gff_file=$(ls {input.ref_dir}/*gff*.gz)
            gff_filename=${gff_file##*/}
            printf "\nGzipped GFF annotation file found in {input.ref_dir}:\n ${gff_filename}\n"
            pigz -p {threads} -dc ${gff_file} > {output.gff}	
        elif [ -s {input.ref_dir}/*.gff* ]; then
            gff_file=$(ls {input.ref_dir}/*.gff*)
            gff_filename=${gff_file##*/}
            printf "\nUnzipped GFF annotation file found in {input.ref_dir}:\n ${gff_filename}\n"
            cp ${gff_file} {output.gff}
        else
            printf "\nNo gff annotation file found in reference directory:\n {input.ref_dir}\n"
            exit 1
        fi
        """

rule check_gtf:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
    output:
        gtf = "{ref_dir}/temp_{data_type}_{ref}.gtf"
    log:
        "logs/tmp_gtf_{data_type}_{ref}.log"
    conda:
        "envs/reference.yaml"
    threads: workflow.cores
    shell:
        """
        if [ -s {input.ref_dir}/*.gtf.gz ]; then
            gtf_file=$(ls {input.ref_dir}/*gtf.gz)
            gtf_filename=${gtf_file##*/}
            printf "\nGzipped GTF annotation file found in {input.ref_dir}:\n ${gtf_filename}\n"
            pigz -p {threads} -dc ${gtf_file} > {ouptut.gtf}	
            gtf={input.ref_dir}/temp_{wildcards.data_type}_{wildcards.ref}.gtf
        elif [ -s {input.ref_dir}/*.gtf ]; then
            gtf_file=$(ls {input.ref_dir}/*.gtf)
            gtf_filename=${gtf_file##*/}
            printf "\nUnzipped GTF annotation file found in {input.ref_dir}:\n ${gtf_filename}\n"
            cp ${gtf_file} {ouptut.gtf}
        else
            printf "\nNo GTF annotation file found in reference directory:\n {input.ref_dir}\n"
            exit 1
        fi
        """
        
rule check_chrom_sizes:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome),
        fasta = "{ref_dir}/temp_{data_type}_{ref}.fa"
    output:
        fasta_index = "{ref_dir}/temp_{data_type}_{ref}.fa.fai",
        chrom_sizes = "{ref_dir}/chrom.sizes"
    log:
        "logs/tmp_chrom_size_{data_type}_{ref}.log"
    conda:
        "envs/reference.yaml"
    shell:
        """
        printf "\nMaking chrom.sizes file for {wildcards.ref}\n"
        samtools faidx {input.fasta}
        cut -f1,2 {output.fasta_index} > {output.chrom.sizes}
        """

rule prep_region_file:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome),
        chrom_sizes = "{ref_dir}/chrom.sizes",
        gff = "{ref_dir}/temp_{data_type}_{ref}.gff"
    output:
        region_file1 = "${data_type}/tracks/${ref}_protein_coding_genes.bed",
        region_file2 = "${data_type}/tracks/${ref}_all_genes.bed"
    log:
        "logs/tmp_region_file_{data_type}_{ref}.log"
    conda:
        "envs/reference.yaml"
    shell:
        """
        printf "\nMaking a bed file with gene coordinates from ${ref}\n"
        awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {print $1,$4-1,$5,$9,".",$7}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        """

rule prep_stat_file:
    output:
        stat_file = "${data_type}/reports/summary_mapping_stats.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """
    
