rule make_ChIP_indices:
      input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome),
        chrom_sizes = "{ref_dir}/chrom.sizes",
        gff = "{ref_dir}/temp_{data_type}_{ref}.gff"
      output:
        indices="{ref_dir}/*.bt2*"
      log:
        "logs/tmp_index_{data_type}_{ref}.log"
      threads: workflow.cores
      shell:
        """
        ### There could be an issue betwene overlapping indices being built between ChIP and TF, to be resolved
        if ls {input.ref_dir}/*.bt2* 1> /dev/null 2>&1; then
            printf "\nBowtie2 index already exists for {wildcards.ref} in {input.ref_dir}\n"
        else
            printf "\nBuilding Bowtie2 index for {wildcards.ref}\n"
            bowtie2-build --threads {threads} {input.fasta} {input.ref_dir}/{wildcards.ref}
        fi
        """
