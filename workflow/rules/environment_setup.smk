def return_log_env(ref_genome, step):
    return os.path.join(REPO_FOLDER,"results","logs",f"tmp_{step}_{ref_genome}.log")

rule prepare_reference:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        gff = "genomes/{ref_genome}/{ref_genome}.gff",
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf",
        chrom_sizes = "genomes/{ref_genome}/chrom.sizes",
        genome_stats = "genomes/{ref_genome}/genome_stats.json",
        taxid = "genomes/{ref_genome}/taxid.json",
        region_files = ["results/combined/bedfiles/{ref_genome}__protein_coding_genes.bed", "results/combined/bedfiles/{ref_genome}__all_genes.bed"],
        logs = lambda wildcards: [ return_log_env(wildcards.ref_genome, step) for step in ["fasta", "gff", "gtf", "chrom_sizes", "genome_stats", "taxid", "region_file"] ]
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
        fasta = lambda wildcards: config["genomes"][wildcards.ref_genome]['fasta_file'],
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
        gff = lambda wildcards: config["genomes"][wildcards.ref_genome]['gff_file'],
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
    input:
        gff = "genomes/{ref_genome}/{ref_genome}.gff"
    output:
        gtf = "genomes/{ref_genome}/{ref_genome}.gtf"
    params:
        gtf = lambda wildcards: config["genomes"][wildcards.ref_genome].get('gtf_file', '<auto>'),
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
        override="{params.gtf}"
        if [[ "$override" != "<auto>" && -n "$override" ]]; then
            if [[ ! -s "$override" ]]; then
                printf "\nGTF file for {params.ref_genome} does not exist:\n$override\n"
                exit 1
            elif [[ "$override" == *.gtf.gz ]]; then
                printf "\nGzipped gtf file found: $override\n"
                pigz -p {threads} -dc "$override" > {output.gtf}
            elif [[ "$override" == *.gtf ]]; then
                printf "\nUnzipped gtf file found: $override\n"
                cp "$override" {output.gtf}
            else
                printf "\nExtension of gtf file unknown, should be .gtf(.gz):\n $override\n"
                exit 1
            fi
        else
            printf "\nDeriving GTF from GFF for {params.ref_genome} using gffread\n"
            gffread -T {input.gff} -o {output.gtf}
            if [[ ! -s {output.gtf} ]]; then
                printf "\nERROR: gffread produced empty GTF from {input.gff}. Supply a GTF explicitly via gtf_file.\n"
                exit 1
            fi
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

rule compute_genome_stats:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        fai = "genomes/{ref_genome}/{ref_genome}.fa.fai"
    output:
        stats = "genomes/{ref_genome}/genome_stats.json"
    log:
        temp(return_log_env("{ref_genome}", "genome_stats"))
    conda: CONDA_ENV
    threads: config["resources"]["compute_genome_stats"]["threads"]
    resources:
        mem_mb=config["resources"]["compute_genome_stats"]["mem_mb"],
        tmp_mb=config["resources"]["compute_genome_stats"]["tmp_mb"],
        qos=config["resources"]["compute_genome_stats"]["qos"]
    shell:
        """
        {{
        printf "\nComputing genome stats for {wildcards.ref_genome}\n"
        total=$(awk '{{s+=$2}} END {{print s}}' {input.fai})
        n_bases=$(awk '!/^>/ {{gsub(/[^Nn]/,""); n+=length}} END {{print n+0}}' {input.fasta})
        effective=$((total - n_bases))
        star_nbases=$(python3 -c "import math; print(min(14, int(math.log2($total)/2 - 1)))")
        printf '{{"total_bases":%d,"n_bases":%d,"effective_size":%d,"star_sa_index_nbases":%d}}\n' \
            "$total" "$n_bases" "$effective" "$star_nbases" > {output.stats}
        printf "Genome stats: total=%d, N=%d, effective=%d, STAR_SAindexNbases=%d\n" \
            "$total" "$n_bases" "$effective" "$star_nbases"
        }} 2>&1 | tee -a "{log}"
        """

rule resolve_taxid:
    output:
        taxid_file = "genomes/{ref_genome}/taxid.json"
    params:
        genus = lambda wildcards: config["genomes"][wildcards.ref_genome].get("genus", ""),
        species = lambda wildcards: config["genomes"][wildcards.ref_genome].get("species", ""),
        override = lambda wildcards: str(config["genomes"][wildcards.ref_genome].get("ncbi_taxid", ""))
    log:
        temp(return_log_env("{ref_genome}", "taxid"))
    conda: CONDA_ENV
    threads: config["resources"]["resolve_taxid"]["threads"]
    resources:
        mem_mb=config["resources"]["resolve_taxid"]["mem_mb"],
        tmp_mb=config["resources"]["resolve_taxid"]["tmp_mb"],
        qos=config["resources"]["resolve_taxid"]["qos"]
    shell:
        """
        {{
        override="{params.override}"
        if [[ -n "$override" && "$override" != "<auto>" ]]; then
            printf "Using user-provided NCBI TaxId: %s\n" "$override"
            printf '{{"ncbi_taxid":"%s"}}\n' "$override" > {output.taxid_file}
        else
            printf "Looking up NCBI TaxId for {params.genus} {params.species}\n"
            taxid=$(datasets summary taxonomy taxon "{params.genus} {params.species}" \
                | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['reports'][0]['taxonomy']['tax_id'])" 2>/dev/null) || true
            if [[ -n "$taxid" ]]; then
                printf "Resolved NCBI TaxId: %s\n" "$taxid"
                printf '{{"ncbi_taxid":"%s"}}\n' "$taxid" > {output.taxid_file}
            else
                printf "WARNING: Could not resolve NCBI TaxId for {params.genus} {params.species}. GO analysis requiring TaxId may fail.\n" >&2
                printf '{{"ncbi_taxid":null}}\n' > {output.taxid_file}
            fi
        fi
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
        
rule download_rfam:
    output:
        cm = "genomes/rfam/Rfam.cm",
        clanin = "genomes/rfam/Rfam.clanin",
        pressed = touch("genomes/rfam/Rfam.cm.pressed")
    log:
        os.path.join(REPO_FOLDER,"results","logs","download_rfam.log")
    conda: CONDA_ENV
    threads: config["resources"]["download_rfam"]["threads"]
    resources:
        mem_mb=config["resources"]["download_rfam"]["mem_mb"],
        tmp_mb=config["resources"]["download_rfam"]["tmp_mb"],
        qos=config["resources"]["download_rfam"]["qos"]
    shell:
        """
        {{
        printf "Downloading Rfam covariance models and clan file\\n"
        mkdir -p genomes/rfam
        curl -fSL https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz \
            | pigz -dc > {output.cm}
        curl -fSL https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin \
            > {output.clanin}

        printf "Pressing covariance model database\\n"
        cmpress -F {output.cm}
        }} 2>&1 | tee -a "{log}"
        """

rule build_structural_rna_db:
    input:
        fasta = "genomes/{ref_genome}/{ref_genome}.fa",
        fai = "genomes/{ref_genome}/{ref_genome}.fa.fai",
        rfam_cm = "genomes/rfam/Rfam.cm",
        rfam_clanin = "genomes/rfam/Rfam.clanin",
        rfam_pressed = "genomes/rfam/Rfam.cm.pressed"
    output:
        structural_fa = "genomes/structural_RNAs/{ref_genome}/structural_rnas.fa",
        tblout = "genomes/structural_RNAs/{ref_genome}/infernal.tblout"
    params:
        threshold = config.get("infernal_threshold", "ga"),
        bin_size_bp = 1000000
    log:
        return_log_env("{ref_genome}", "structural_rna")
    conda: CONDA_ENV
    threads: config["resources"]["build_structural_rna_db"]["threads"]
    resources:
        mem_mb=config["resources"]["build_structural_rna_db"]["mem_mb"],
        tmp_mb=config["resources"]["build_structural_rna_db"]["tmp_mb"],
        qos=config["resources"]["build_structural_rna_db"]["qos"]
    shell:
        """
        {{
        printf "Building structural RNA database for {wildcards.ref_genome} via Infernal\\n"

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT
        chroms_dir="$workdir/chroms"
        mkdir -p "$chroms_dir"

        # Compute search space: total bases * 2 (both strands), in Mb
        search_space=$(awk '{{s+=$2}} END {{printf "%d", s*2/1000000}}' {input.fai})
        printf "Search space (-Z): %s Mb\\n" "$search_space"

        # Threshold flag
        threshold="{params.threshold}"
        if [[ "$threshold" == "ga" ]]; then
            threshold_flag="--cut_ga"
        else
            threshold_flag="-E $threshold --incE $threshold"
        fi
        printf "Infernal threshold: %s (flag: %s)\\n" "$threshold" "$threshold_flag"

        # Split genome: large contigs (>= bin_size) get their own file,
        # small contigs (< bin_size) are binned together
        bin_size={params.bin_size_bp}
        bin_idx=0
        bin_total=0
        bin_file="$chroms_dir/bin_${{bin_idx}}.fa"

        while IFS=$'\\t' read -r chrom size rest; do
            if [[ "$size" -ge "$bin_size" ]]; then
                samtools faidx {input.fasta} "$chrom" > "$chroms_dir/${{chrom}}.fa"
            else
                if [[ "$bin_total" -gt 0 && $(( bin_total + size )) -gt "$bin_size" ]]; then
                    bin_idx=$(( bin_idx + 1 ))
                    bin_total=0
                    bin_file="$chroms_dir/bin_${{bin_idx}}.fa"
                fi
                samtools faidx {input.fasta} "$chrom" >> "$bin_file"
                bin_total=$(( bin_total + size ))
            fi
        done < {input.fai}

        n_chunks=$(find "$chroms_dir" -name '*.fa' | wc -l)
        printf "Split genome into %d chunks (bin_size=%d bp)\\n" "$n_chunks" "$bin_size"

        # Determine parallelism: allocate threads across chunks
        if [[ "$n_chunks" -le 0 ]]; then
            printf "ERROR: No genome chunks produced\\n" >&2
            exit 1
        fi
        threads_per=$(( {threads} / n_chunks ))
        if [[ "$threads_per" -lt 1 ]]; then threads_per=1; fi
        n_parallel=$(( {threads} / threads_per ))
        if [[ "$n_parallel" -lt 1 ]]; then n_parallel=1; fi
        printf "Running %d parallel cmscan jobs with %d threads each\\n" "$n_parallel" "$threads_per"

        # Run cmscan per chunk
        find "$chroms_dir" -name '*.fa' | \
            xargs -P "$n_parallel" -I {{}} bash -c '
                cmscan --cpu '"$threads_per"' \
                    -Z '"$search_space"' \
                    '"$threshold_flag"' \
                    --rfam --nohmmonly --fmt 2 \
                    --tblout "{{}}.tblout" \
                    --clanin {input.rfam_clanin} \
                    {input.rfam_cm} \
                    {{}} > /dev/null
            '

        # Merge tblout results: keep header from first file, concatenate data lines,
        # remove lower-scoring clan overlaps (marked with " = " in olp column)
        printf "Merging and filtering Infernal results\\n"
        head -2 "$(find "$chroms_dir" -name '*.tblout' | head -1)" > {output.tblout}
        cat "$chroms_dir"/*.tblout | grep -v '^#' | grep -v ' = ' \
            | sort -k16,16g >> {output.tblout} || true

        # Count significant hits
        n_hits=$(grep -vc '^#' {output.tblout} || echo 0)
        printf "Found %d structural RNA loci\\n" "$n_hits"

        if [[ "$n_hits" -eq 0 ]]; then
            printf "WARNING: No structural RNAs found for {wildcards.ref_genome}. Creating empty FASTA.\\n" >&2
            touch {output.structural_fa}
        else
            # Convert tblout to BED (0-based coords)
            awk '$0 !~ /^#/ {{
                chrom = $4; start = $10; end = $11; name = $2; strand = $12
                if (strand == "+") {{ s = start - 1; e = end }}
                else {{ s = end - 1; e = start }}
                printf "%s\\t%d\\t%d\\t%s::%s:%d-%d(%s)\\t0\\t%s\\n",
                    chrom, s, e, name, chrom, s, e, strand, strand
            }}' {output.tblout} \
                | sort -k1,1 -k2,2n > "$workdir/structural.bed"

            # Extract sequences
            bedtools getfasta -fi {input.fasta} -bed "$workdir/structural.bed" \
                -s -name -fo {output.structural_fa}
        fi

        printf "Structural RNA database built: %s\\n" "{output.structural_fa}"
        }} 2>&1 | tee -a "{log}"
        """

rule check_te_file:
    output:
        te_file = "genomes/{ref_genome}/{ref_genome}__TE_file.bed"
    params:
        te_file = lambda wildcards: config["genomes"][wildcards.ref_genome]['te_file'],
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