# A. thaliana ColCEN Test Case Design

## Objective

A multi-omics Arabidopsis test case covering 6 assay types: ChIP_broad
(CenH3, H3K9me2), ATAC, EMseq, PBAT, and dmC. Uses the full Col-CEN v1.2
reference genome (~133 Mbp T2T assembly). All inputs are publicly downloadable:
genome from GitHub, ChIP/ATAC from SRA, EMseq from SRA, PBAT from DDBJ,
dmC modBAMs from lemna.org.

## Reference Genome

**Assembly**: Col-CEN v1.2 (Naish et al. 2021, *Science* 374:abk3085)
**Source**: https://github.com/schatzlab/Col-CEN/tree/main/v1.2

| File | URL |
|------|-----|
| FASTA | `https://github.com/schatzlab/Col-CEN/raw/main/v1.2/Col-CEN_v1.2.fasta.gz` |
| GFF3 (genes) | `https://github.com/schatzlab/Col-CEN/raw/main/v1.2/Col-CEN_v1.2_genes.tair10.gff3.gz` |
| GFF3 (TEs) | `https://github.com/schatzlab/Col-CEN/raw/main/v1.2/t2t-col.20201227.fasta.mod.EDTA.TEanno.gff3.gz` |

The TE annotation is in EDTA GFF3 format (auto-converted to BED6 by the
pipeline using the `ID=` attribute as the name column).

## Dataset 1: ChIP-seq (CenH3, H3K9me2) — Shimada et al. 2024

**Source**: Shimada et al. (2024) "HIRA-dependent H3.3 deposition facilitates
PRC2-mediated H3K27me3 and is essential for RdDM-independent heterochromatin
maintenance." *Nature Plants* 10, 1453-1467. doi:10.1038/s41477-024-01779-7

**GEO**: [GSE132005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132005)

### Experimental design

- **Factor**: `genotype` — WT, rdr126ddm1, rdr126ddm1hp5
- **Tissue**: leaf (all samples)
- **Replicates**: 2 per genotype × mark × IP/Input
- **Layout**: PE, Illumina

### CenH3 ChIP (12 samples)

| Genotype | Rep | IP SRR | Input SRR |
|----------|-----|--------|-----------|
| WT | rep1 | SRR28453410 | SRR28453412 |
| WT | rep2 | SRR28453411 | SRR28453413 |
| rdr126ddm1 | rep1 | SRR28453414 | SRR28453416 |
| rdr126ddm1 | rep2 | SRR28453415 | SRR28453417 |
| rdr126ddm1hp5 | rep1 | SRR28453418 | SRR28453420 |
| rdr126ddm1hp5 | rep2 | SRR28453419 | SRR28453421 |

### H3K9me2 ChIP (8 samples)

No WT H3K9me2 exists in GSE132005. Input samples are distinct from CenH3
Inputs (different GEO entries, different library preparations).

| Genotype | Rep | IP SRR | Input SRR |
|----------|-----|--------|-----------|
| rdr126ddm1 | rep1 | SRR28453402 | SRR28453404 |
| rdr126ddm1 | rep2 | SRR28453403 | SRR28453405 |
| rdr126ddm1hp5 | rep1 | SRR28453406 | SRR28453408 |
| rdr126ddm1hp5 | rep2 | SRR28453407 | SRR28453409 |

## Dataset 2: dmC (ONT modBAM) — Shimada et al. 2024

**Source**: Same GEO series (GSE132005). ONT native methylation data deposited
as modBAMs. These are hosted on lemna.org because the SRA deposits were
mislabeled as bisulfite-seq.

- **Tissue**: leaf
- **Replicates**: 1 per genotype
- **Layout**: SE (long reads)

| Genotype | URL |
|----------|-----|
| rdr126ddm1 | `https://lemna.org/data/shimada2024/SRR28453434.bam` |
| rdr126ddm1het | `https://lemna.org/data/shimada2024/SRR28453435.bam` |
| rdr126ddm1hp5 | `https://lemna.org/data/shimada2024/SRR28453436.bam` |
| rdr126ddm1hethp5 | `https://lemna.org/data/shimada2024/SRR28453437.bam` |

## Dataset 3: EMseq — Trasser et al. 2024

**Source**: Trasser et al. (2024) "PTGS is dispensable for the initiation of
epigenetic silencing of an active transposon in Arabidopsis." *EMBO Reports*
25(12), 5780-5809. doi:10.1038/s44319-024-00304-5

**BioProject**: [PRJNA1111825](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1111825)

### Experimental design

- **Factor**: `genotype` — WT (Col-0) vs rdr6
- **Tissue**: flower (inflorescence)
- **Replicates**: 2 per genotype
- **Layout**: PE 150bp, Illumina NovaSeq 6000
- **Library**: NEBNext Enzymatic Methyl-seq (EM-seq) Kit

| Genotype | Rep | SRR |
|----------|-----|-----|
| WT (Col-0) | rep1 | SRR29036840 |
| WT (Col-0) | rep2 | SRR29036839 |
| rdr6 | rep1 | SRR29036836 |
| rdr6 | rep2 | SRR29036835 |

The full dataset contains 12 samples (including rdr6 F6 epigenetic variants
with active/silenced EVD transposon states). Only the parental genotypes
(Col-0 and rdr6) are used here.

## Dataset 4: PBAT — Takei et al. 2024

**Source**: Takei et al. (2024) *Plant Physiology* 195(2), 1333.
https://academic.oup.com/plphys/article/195/2/1333/7623186

**BioProject**: [PRJDB14218](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB14218)
(Watanabe Lab, University of Tokyo)

### Experimental design

- **Factor**: `genotype` — WT (Col-0) vs tarp1tarp2 double mutant
- **Tissue**: seedling (10-day-old)
- **Replicates**: 2 per genotype
- **Layout**: PE, Illumina NovaSeq 6000
- **Library**: PBAT (post-bisulfite adaptor tagging)

Uses DDBJ accessions (DRR prefix); fasterq-dump and ENA download handle these.

| Genotype | Rep | DRR |
|----------|-----|-----|
| WT (Col-0) | rep1 | DRR400324 |
| WT (Col-0) | rep2 | DRR400325 |
| tarp1tarp2 | rep1 | DRR400326 |
| tarp1tarp2 | rep2 | DRR400327 |

The full dataset also contains sRNA-seq samples (DRR400322-23, DRR528207-08)
which are not used here.

## Dataset 5: ATAC-seq — Crisp et al. 2020

**Source**: Crisp et al. (2020) "Stable unmethylated DNA demarcates expressed
genes and their cis-regulatory space in plant genomes." *PNAS* 117(38),
23991-24000. doi:10.1073/pnas.2010250117

**GEO**: [GSE155023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155023)

### Experimental design

- **Factor**: `genotype` — WT, ddm1, nrpd1
- **Tissue**: leaf
- **Replicates**: 2 per genotype (subset of 4 available)
- **Layout**: PE

| Genotype | Rep | SRR |
|----------|-----|-----|
| WT | rep1 | SRR12362020 |
| WT | rep2 | SRR12362021 |
| ddm1 | rep1 | SRR12362026 |
| ddm1 | rep2 | SRR12362027 |
| nrpd1 | rep1 | SRR12362044 |
| nrpd1 | rep2 | SRR12362045 |

## Assay Coverage Matrix

| Assay | Genotypes | Tissue | Layout | Reps | Samples | Source |
|-------|-----------|--------|--------|------|---------|--------|
| ChIP_broad (CenH3) | WT, rdr126ddm1, rdr126ddm1hp5 | leaf | PE | 2 | 12 | SRA (GSE132005) |
| ChIP_broad (H3K9me2) | rdr126ddm1, rdr126ddm1hp5 | leaf | PE | 2 | 8 | SRA (GSE132005) |
| dmC (ONT modBAM) | rdr126ddm1, rdr126ddm1het, rdr126ddm1hp5, rdr126ddm1hethp5 | leaf | SE | 1 | 4 | lemna.org URLs |
| EMseq | WT, rdr6 | flower | PE | 2 | 4 | SRA (PRJNA1111825) |
| PBAT | WT, tarp1tarp2 | seedling | PE | 2 | 4 | DDBJ (PRJDB14218) |
| ATAC | WT, ddm1, nrpd1 | leaf | PE | 2 | 6 | SRA (GSE155023) |

**Total**: 38 samples across 6 assay types.

## Known Design Considerations

### Multiple tissues across assays

The datasets come from different labs and use different tissues: leaf (ChIP,
ATAC, dmC), flower/inflorescence (EMseq), and seedling (PBAT). This means
cross-assay comparisons (e.g. combined heatmaps) will reflect tissue
differences in addition to genotype effects. This is acceptable for a pipeline
test case — the goal is to exercise all code paths, not to produce
biologically interpretable cross-assay comparisons.

### Multiple genotype sets

Each dataset uses a different set of genotypes. Only WT (Col-0) appears across
all datasets, enabling at least some cross-assay visualization at the WT level.
The mutant genotypes (rdr126ddm1, rdr6, tarp1tarp2, ddm1, nrpd1) are
dataset-specific.

### dmC samples have no replicates

Only one replicate per genotype is available for the dmC samples. This is
sufficient to exercise the dmC pipeline path (modBAM download, alignment,
modkit pileup, bedMethyl conversion, CX report) but DMR calling between
dmC genotypes will have limited statistical power.

### PBAT uses DDBJ accessions

The PBAT samples use DRR-prefix accessions from DDBJ rather than the more
common SRR (NCBI) or ERR (ENA) prefixes. The pipeline's SRA detection regex
was expanded to accept all three prefixes (`[SED]RR\d+`).

### Full genome (not subsetted)

Unlike the hg38_chr21 test case, this test case uses the full ColCEN genome
(~133 Mbp). Full-genome runs are slower but avoid the complexity of subsetting
plant data to a single chromosome (plants lack the clean chr21-like subsetting
target that works well for human data).
