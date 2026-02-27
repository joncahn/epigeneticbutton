# ColCEN Chr5 Minimal Test Reference

This directory contains a minimal reference genome for integration testing,
extracted from the full Arabidopsis thaliana Col-CEN v1.2 reference.

## Contents

| File | Description | Size |
|------|-------------|------|
| `ColCEN_Chr5.fa.gz` | Chr5 sequence only | ~8 MB |
| `ColCEN_Chr5.gff3.gz` | Gene annotations for Chr5 | ~1 MB |
| `ColCEN_Chr5.gtf.gz` | Gene annotations (GTF format) | ~700 KB |
| `chrom.sizes` | Chromosome sizes file | 18 bytes |

## Source

Extracted from:
- `/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/ColCEN/Col-CEN_v1.2.fasta.gz`
- `/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/ColCEN/Col-CEN_v1.2_genes.tair10.gff3.gz`
- `/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/ColCEN/Col-CEN_v1.2_genes.tair10.gtf.gz`

## Chromosome Details

- **Chr5**: 29,480,885 bp (~29.5 Mb)

## Usage

Use with `tests/integration/data/test_config_chr5.yaml` for integration testing:

```bash
snakemake --configfile tests/integration/data/test_config_chr5.yaml \
    results/mC/tracks/rdr126ddm1_leaf_dmC_rep1__CG.bw
```

## Adding Test Samples

Update `tests/integration/data/test_samples_chr5.tsv` with actual dmC modBAM paths:

```
Sample_ID	Assay	Genome	Levels	Replicate_ID	Read_files	Read_layout	IP_target	Control
Col0_leaf_dmC_rep1	dmC	ColCEN_Chr5	genotype:Col0,tissue:leaf	rep1	/path/to/modbam.bam	SE
```

For best testing, the modBAM should be filtered to Chr5 reads only to keep file sizes small.
