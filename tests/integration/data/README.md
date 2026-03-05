# S. pombe Test Data

## Source

ChIP-seq and sRNA-seq samples are from **Kim, Roche, Bhattacharjee et al. (2024)**
"Clr4SUV39H1 ubiquitination and non-coding RNA mediate transcriptional silencing
of heterochromatin via Swi6 phase separation." *Nature Communications* 15, 9384.
doi:10.1038/s41467-024-53417-9 — GEO accession
[GSE156069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156069).

RNA-seq and H3K4me3 ChIP-seq samples are from separate public datasets (see
sample sheet for SRA accessions).

## Sample Sheet

`test_samples_pombe.tsv` — 18 samples across 4 assay types:

| Assay | Genotypes | Replicates | Notes |
|-------|-----------|------------|-------|
| ChIP_broad (H3K9me2) | WT, dcr1 | 2 (WT), 1 (dcr1) | Shallow MiSeq libraries |
| ChIP_broad (H3K9me3) | WT, dcr1 | 2 (WT), 1 (dcr1) | dcr1 extremely shallow (436K reads) |
| ChIP_narrow (H3K4me3) | WT | 2 | Good depth; exercises full IDR path |
| RNAseq | WT, dcr1 | 2 each | |
| sRNA | WT, dcr1 | 2 (WT), 1 (dcr1) | |

## Known Caveats

### Shallow H3K9me2/me3 ChIP-seq libraries

The H3K9me2 and H3K9me3 samples from GSE156069 were sequenced on an Illumina
MiSeq with low read depth:

| Sample | SRA | Read pairs |
|--------|-----|------------|
| WT H3K9me2 rep1 | SRR20678305 | ~2M |
| WT H3K9me2 rep2 | SRR20678333 | ~2M |
| WT H3K9me3 rep1 | SRR20678311 | ~4.5M |
| WT H3K9me3 rep2 | SRR20678336 | ~4.5M |
| dcr1 H3K9me2 rep1 | SRR20678308 | ~2M |
| dcr1 H3K9me3 rep1 | SRR20678314 | **~436K** |

This is well below the 10M+ reads typically needed for reliable broad-peak
ChIP-seq. As a result, MACS2 calls very few peaks (3-74 per replicate for
H3K9me marks), which causes:

- **IDR failure** for WT H3K9me2 (too few peaks to model reproducibility)
- **MAnorm skip** for H3K9me3 WT-vs-dcr1 (dcr1 has only 1 selected peak)

These outcomes are expected and exercise the pipeline's graceful-failure paths.

### Biology: few H3K9me peaks expected in S. pombe

Even with adequate sequencing depth, S. pombe has relatively few H3K9me2/me3
peaks. Heterochromatin is concentrated at three main loci: centromeres,
telomeres, and the mating-type region. In dcr1 deletion mutants, the RNAi
pathway is disrupted, but H3K9me2 is partially retained (~30-55% of WT) via
RNAi-independent mechanisms (SHREC, Seb1). The H3K9me2-to-me3 transition is
more severely impaired in dcr1 mutants.

### H3K4me3: the well-behaved mark

H3K4me3 (WT, 2 replicates) has adequate depth and produces ~4200-4500 peaks
per replicate, with **1681 IDR-passing peaks**. This mark fully exercises the
IDR, selected-peaks, and differential-analysis pipeline paths.

### WCE control shared across assay types

`WT_WCE_rep1` (SRR5445712) is used as the ChIP input control for all broad-peak
samples. It comes from a separate SRA submission but the same organism/strain.
`WT_Input_rep1` (SRR31073785) is used for narrow-peak (H3K4me3) samples.
