rule stat_file:
    output:
        stat_file = "{data_type}/reports/summary_mapping_stats.txt"
    shell:
        """
        if [ ! -s {output.stat_file} ]; then
            printf "Line\tTissue\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        fi
        """

rule make_ChIP_indices:
      input:
        ref_dir = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome),
        chrom_sizes = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome, "chrom.sizes"),
        gff = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome, f"temp_{wildcards.ref}.gff")
      output:
        indices = lambda wildcards: os.path.join(config["ref_path"], wildcards.ref_genome, f"temp_{wildcards.ref}*.bt2*")
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
