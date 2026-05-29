#!/usr/bin/env bash
# Download and prepare hg38 chr21 reference genome and annotation
# One-time data prep for the human chr21 integration test
set -euo pipefail

OUTDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$OUTDIR"

echo "=== Preparing hg38 chr21 reference ==="

# ---------------------------------------------------------------------------
# 1. Download full hg38 FASTA (needed for alignment before subsetting)
# ---------------------------------------------------------------------------
if [[ ! -f hg38_full.fa ]]; then
    echo "Downloading full hg38 FASTA from UCSC ..."
    curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
        -o hg38_full.fa.gz
    echo "  Decompressing ..."
    gunzip hg38_full.fa.gz
    echo "  Indexing ..."
    samtools faidx hg38_full.fa
    echo "  Done: hg38_full.fa (~3.1 GB)"
else
    echo "  hg38_full.fa already exists"
    [[ -f hg38_full.fa.fai ]] || samtools faidx hg38_full.fa
fi

# ---------------------------------------------------------------------------
# 2. Extract chr21 FASTA (pipeline test reference)
# ---------------------------------------------------------------------------
if [[ ! -f hg38_chr21.fa.gz ]]; then
    echo "Extracting chr21 from full FASTA ..."
    samtools faidx hg38_full.fa chr21 | gzip > hg38_chr21.fa.gz
    echo "  Done: hg38_chr21.fa.gz"
else
    echo "  hg38_chr21.fa.gz already exists"
fi

# ---------------------------------------------------------------------------
# 3. Download GENCODE annotation and subset to chr21
# ---------------------------------------------------------------------------
GENCODE_VER="v46"
GENCODE_GFF="gencode.${GENCODE_VER}.annotation.gff3.gz"
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VER#v}/${GENCODE_GFF}"

if [[ ! -f hg38_chr21.gff3.gz ]]; then
    echo "Downloading GENCODE ${GENCODE_VER} GFF3 ..."
    if [[ ! -f "$GENCODE_GFF" ]]; then
        curl -sL "$GENCODE_URL" -o "$GENCODE_GFF"
    fi
    echo "  Subsetting to chr21 ..."
    zcat "$GENCODE_GFF" | awk '$1 == "chr21" || /^#/' | gzip > hg38_chr21.gff3.gz
    echo "  Done: hg38_chr21.gff3.gz"
    # Clean up full annotation
    rm -f "$GENCODE_GFF"
else
    echo "  hg38_chr21.gff3.gz already exists"
fi

# ---------------------------------------------------------------------------
# 3. Validate outputs
# ---------------------------------------------------------------------------
echo ""
echo "=== Validation ==="
echo -n "Full FASTA: "; du -h hg38_full.fa | cut -f1
echo -n "Full FASTA chromosomes: "; grep -c '^>' hg38_full.fa
echo -n "Chr21 FASTA: "; zcat hg38_chr21.fa.gz | grep -o '^>.*' || true
echo -n "Chr21 FASTA size: "; du -h hg38_chr21.fa.gz | cut -f1
echo -n "GFF3 features: "; zcat hg38_chr21.gff3.gz | grep -vc '^#' || true
echo -n "GFF3 size: "; du -h hg38_chr21.gff3.gz | cut -f1

echo ""
echo "=== Reference preparation complete ==="
echo "Full genome (for alignment): ${OUTDIR}/hg38_full.fa"
echo "Chr21 FASTA (pipeline test): ${OUTDIR}/hg38_chr21.fa.gz"
echo "Chr21 GFF3 (pipeline test):  ${OUTDIR}/hg38_chr21.gff3.gz"
