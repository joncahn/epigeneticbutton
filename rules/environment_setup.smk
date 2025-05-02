# Rule to prepare reference genome for each data type
rule prepare_reference:
    input:
        fasta = "{ref_dir}/temp_{ref}.fa"
        gff = 
        gtf = 
        chrom_sizes =
        region_file = 
        logs = ["logs/tmp_fasta_{data_type}_{ref}.log", "logs/tmp_gff_{data_type}_{ref}.log", "logs/tmp_gtf_{data_type}_{ref}.log", "logs/tmp_index_{data_type}_{ref}.log"]
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
    shell:
        """
        # Search for fasta file
        ref_dir={input.ref_dir}
        if [ -s ${ref_dir}/*.fa.gz ]; then
            fa_file=$(ls ${ref_dir}/*.fa.gz)
            fa_filename=${fa_file##*/}
            printf "\nGzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${datatype}_${ref}.fa
            fasta=${ref_dir}/temp_${datatype}_${ref}.fa
        elif [ -s ${ref_dir}/*.fa ]; then
            fa_file=$(ls ${ref_dir}/*.fa)
            fa_filename=${fa_file##*/}
            printf "\nUnzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            fasta=${fa_file}
        elif [ -s ${ref_dir}/*.fasta.gz ]; then
            fa_file=$(ls ${ref_dir}/*.fasta.gz)
            fa_filename=${fa_file##*/}
            printf "\nGzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            pigz -p ${threads} -dc ${fa_file} > ${ref_dir}/temp_${datatype}_${ref}.fa
            fasta=${ref_dir}/temp_${datatype}_${ref}.fa
        elif [ -s ${ref_dir}/*.fasta ]; then
            fa_file=$(ls ${ref_dir}/*.fasta)
            fa_filename=${fa_file##*/}
            printf "\nUnzipped fasta file found in ${ref_dir}:\n ${fa_filename}\n"
            fasta=${fa_file}
        else
            printf "\nNo fasta file found in reference directory:\n ${ref_dir}\n"
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
    shell:
        """
        ref_dir={input.ref_dir}
        if [ -s ${ref_dir}/*.gff*.gz ]; then
            gff_file=$(ls ${ref_dir}/*gff*.gz)
            gff_filename=${gff_file##*/}
            printf "\nGzipped GFF annotation file found in ${ref_dir}:\n ${gff_filename}\n"
            pigz -p ${threads} -dc ${gff_file} > ${ref_dir}/temp_${datatype}_${ref}.gff	
            gff=${ref_dir}/temp_${datatype}_${ref}.gff
        elif [ -s ${ref_dir}/*.gff* ]; then
            gff_file=$(ls ${ref_dir}/*.gff*)
            gff_filename=${gff_file##*/}
            printf "\nUnzipped GFF annotation file found in ${ref_dir}:\n ${gff_filename}\n"
            gff=${gff_file}
        else
            printf "\nNo gff annotation file found in reference directory:\n ${ref_dir}\n"
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
    shell:
        """
        ref_dir={input.ref_dir}
        if [ -s ${ref_dir}/*.gtf.gz ]; then
            gtf_file=$(ls ${ref_dir}/*gtf.gz)
            gtf_filename=${gtf_file##*/}
            printf "\nGzipped GTF annotation file found in ${ref_dir}:\n ${gtf_filename}\n"
            pigz -p ${threads} -dc ${gtf_file} > ${ref_dir}/temp_${datatype}_${ref}.gtf	
            gtf=${ref_dir}/temp_${datatype}_${ref}.gtf
        elif [ -s ${ref_dir}/*.gtf ]; then
            gtf_file=$(ls ${ref_dir}/*.gtf)
            gtf_filename=${gtf_file##*/}
            printf "\nUnzipped GTF annotation file found in ${ref_dir}:\n ${gtf_filename}\n"
            gtf=${gtf_file}
        else
            printf "\nNo GTF annotation file found in reference directory:\n ${ref_dir}\n"
            exit 1
        fi
        """
        
rule check_chrom_sizes:
    input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome)
        fasta = "{ref_dir}/temp_{data_type}_{ref}.fa"
    output:
        chrom_sizes = "{ref_dir}/chrom.sizes"
    conda:
        "envs/reference.yaml"
    shell:
        """
        ref_dir={input.ref_dir}
        if [ ! -s ${ref_dir}/chrom.sizes ]; then
            printf "\nMaking chrom.sizes file for ${ref}\n"
            samtools faidx {input.fasta}
            cut -f1,2 ${fasta}.fai > ${ref_dir}/chrom.sizes
        fi
        """


