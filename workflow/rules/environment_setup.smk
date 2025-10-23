def return_log_env(ref_genome, step):
    return os.path.join(REPO_FOLDER,"results","logs",f"tmp_{step}_{ref_genome}.log")

rule prepare_reference:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gff = "genomes/{ref_genome}/{ref_genome}.gff",
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        region_files = ["results/combined/bedfiles/{ref_genome}__protein_coding_genes.bed", "results/combined/bedfiles/{ref_genome}__all_genes.bed"],
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

rule check_fasta:
    output:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa"
    params:
        fasta = lambda wildcards: config[wildcards.ref_genome]['fasta_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "fasta"))
    conda: CONDA_ENV
    threads: config["resources"]["check_fasta"]["threads"]
    resources:
        mem_mb=config["resources"]["check_fasta"]["mem_mb"],
        tmp_mb=config["resources"]["check_fasta"]["tmp_mb"],
        qos=config["resources"]["check_fasta"]["qos"]
    shell:
        """
        {{
        if [[ ! -s {params.fasta} ]]; then
            printf "\nFasta file for {params.ref_genome} does not exist:\n{params.fasta}\n"
            exit 1
        elif [[ {params.fasta} == *.fa.gz || {params.fasta} == *.fasta.gz ]]; then
            printf "\nGzipped fasta file found: {params.fasta}\n"
            pigz -p {threads} -dc {params.fasta} > {output.fasta}
        elif [[ {params.fasta} == *.fa || {params.fasta} == *.fasta ]]; then
            printf "\nUnzipped fasta file found: {params.fasta}\n"
            cp {params.fasta} {output.fasta}
        else
            printf "\nExtension of fasta file unknown, should be .fasta(.gz) or .fa(.gz):\n {params.fasta}\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """
        
rule check_gff:
    output:
        gff = "genomes/{ref_genome}/{ref_genome}.gff"
    params:
        gff = lambda wildcards: config[wildcards.ref_genome]['gff_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "gff"))
    conda: CONDA_ENV
    threads: config["resources"]["check_gff"]["threads"]
    resources:
        mem_mb=config["resources"]["check_gff"]["mem_mb"],
        tmp_mb=config["resources"]["check_gff"]["tmp_mb"],
        qos=config["resources"]["check_gff"]["qos"]
    shell:
        """
        {{
        if [[ ! -s {params.gff} ]]; then
            printf "\nGFF file for {params.ref_genome} does not exist:\n{params.gff}\n"
            exit 1
        elif [[ {params.gff} == *.gff*.gz ]]; then
            printf "\nGzipped gff file found: {params.gff}\n"
            pigz -p {threads} -dc {params.gff} > {output.gff}
        elif [[ {params.gff} == *.gff* ]]; then
            printf "\nUnzipped gff file found: {params.gff}\n"
            cp {params.gff} {output.gff}
        else
            printf "\nExtension of gff file unknown, should be .gff*(.gz):\n {params.gff}\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule check_gtf:
    output:
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf"
    params:
        gtf = lambda wildcards: config[wildcards.ref_genome]['gtf_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "gtf"))
    conda: CONDA_ENV
    threads: config["resources"]["check_gtf"]["threads"]
    resources:
        mem_mb=config["resources"]["check_gtf"]["mem_mb"],
        tmp_mb=config["resources"]["check_gtf"]["tmp_mb"],
        qos=config["resources"]["check_gtf"]["qos"]
    shell:
        """
        {{
        if [[ ! -s {params.gtf} ]]; then
            printf "\nGTF file for {params.ref_genome} does not exist:\n{params.gtf}\n"
            exit 1
        elif [[ {params.gtf} == *.gtf.gz ]]; then
            printf "\nGzipped gtf file found: {params.gtf}\n"
            pigz -p {threads} -dc {params.gtf} > {output.gtf}
        elif [[ {params.gtf} == *.gtf ]]; then
            printf "\nUnzipped gtf file found: {params.gtf}\n"
            cp {params.gtf} {output.gtf}
        else
            printf "\nExtension of gtf file unknown, should be .gtf(.gz):\n {params.gtf}\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
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
        mem_mb=config["resources"]["check_chrom_sizes"]["mem_mb"],
        tmp_mb=config["resources"]["check_chrom_sizes"]["tmp_mb"],
        qos=config["resources"]["check_chrom_sizes"]["qos"]
    shell:
        """
        {{
        printf "\nMaking chrom.sizes file for {params.ref_genome}\n"
        samtools faidx {input.fasta}
        cut -f1,2 {output.fasta_index} > {output.chrom_sizes}
        }} 2>&1 | tee -a "{log}"
        """

rule prep_region_file:
    input:
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        gff = "genomes/{ref_genome}/{ref_genome}.gff"
    output:
        region_file1 = "results/combined/bedfiles/{ref_genome}__protein_coding_genes.bed",
        region_file2 = "results/combined/bedfiles/{ref_genome}__all_genes.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "region_file"))
    conda: CONDA_ENV
    threads: config["resources"]["prep_region_file"]["threads"]
    resources:
        mem_mb=config["resources"]["prep_region_file"]["mem_mb"],
        tmp_mb=config["resources"]["prep_region_file"]["tmp_mb"],
        qos=config["resources"]["prep_region_file"]["qos"]
    shell:
        """
        {{
        printf "\nMaking a bed file with gene coordinates from {params.ref_genome}\n" >> {log} 2>&1
        awk -v OFS="\t" '$3=="gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        }} 2>&1 | tee -a "{log}"
        """
        
rule check_te_file:
    output:
        te_file = "genomes/{ref_genome}/{ref_genome}__TE_file.bed"
    params:
        te_file = lambda wildcards: config[wildcards.ref_genome]['te_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "TEs"))
    conda: CONDA_ENV
    threads: config["resources"]["check_te_file"]["threads"]
    resources:
        mem_mb=config["resources"]["check_te_file"]["mem_mb"],
        tmp_mb=config["resources"]["check_te_file"]["tmp_mb"],
        qos=config["resources"]["check_te_file"]["qos"]
    shell:
        """
        {{
        if [[ ! -s {params.te_file} ]]; then
            printf "\nThe bed file of TEs for {wildcards.ref_genome} does not exist:\n {params.te_file}\n"
            exit 1
        elif [[ {params.te_file} == *.bed.gz ]]; then
            printf "\nGzipped TE file found: {params.te_file}\n"
            pigz -p {threads} -dc {params.te_file} > {output.te_file}
        elif [[ {params.te_file} == *.bed ]]; then
            printf "\nUnzipped TE file found: {params.te_file}\n"
            cp {params.te_file} {output.te_file}
        else
            printf "\nExtension of bed file of TEs unknown, should be .bed(.gz):\n {params.te_file}\n"
            exit 1
        fi
        tot=$(cat {output.te_file} | wc -l)
        unique=$(cat {output.te_file} | cut -f4 | sort -u | wc -l)
        if [[ ${{unique}} -ne ${{tot}} ]]; then
            printf "\nNot all the names of TEs are unique. This is required for follow-up. Remove redundant rows or add unique identifiers.\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """