# EPICC Pipeline Benchmarks

Execution profiles for the integration test cases, generated with `dev/profile_snakemake_log.py`.

**Interactive profile report**: [profile_test_cases.html](profile_test_cases.html) (tabbed, with Gantt timelines)

## H. sapiens hg38 chr21 — SLURM cluster

**Date**: 2026-03-12
**Environment**: Elzar SLURM cluster (bamdev2), up to 16 concurrent jobs
**Dataset**: 33 samples across 6 assay types (ChIP_broad, ChIP_narrow, ATAC, RNAseq/RAMPAGE, sRNA, WGBS, dmC)
**Reference**: hg38 chr21 (~46 Mb)

| Metric | Value |
|--------|-------|
| Wall time | 2h 0m |
| Total jobs | 428 |
| Total CPU time | 136h 34m |
| Avg parallelism | 68.2x |

### Bottlenecks

| Rule | Total CPU | Max single job | Notes |
|------|-----------|----------------|-------|
| bismark_map_pe | 4h 23m | 1h 27m (colon_WGBS_D2, 12.8M reads) | Alignment + dedup + methylation extraction |
| make_srna_stranded_bigwigs | 1h 4m | 11m 26s | 24 jobs (per-sample × per-size) |
| get_fastq_pe + get_fastq_se | 1h 52m | 5m 33s | Network-bound SRA downloads |
| STAR_map_pe | 38m 32s | 7m 29s | RNA-seq alignment |
| process_fastq_pe | 44m 26s | 4m 3s | fastp trimming |

### Phase breakdown

| Phase | Total CPU | % of wall |
|-------|-----------|-----------|
| Other (misc analysis) | 9h 33m | 477% |
| Download | 3h 1m | 151% |
| Tracks (bigwigs) | 2h 5m | 104% |
| Alignment | 1h 27m | 73% |
| QC | 1h 9m | 57% |
| sRNA analysis | 54m | 45% |
| Peak calling | 47m | 39% |
| Differential | 35m | 29% |

### Key observations

- **Bismark dominates the critical path.** The longest single job (colon_WGBS_D2) runs for 1h 27m. With `--multicore 8` (16 threads, `threads // 2`), bismark uses ~8 cores effectively. The methylation extraction step (`bismark_methylation_extractor --cytosine_report --CX`) is the slowest sub-phase within each bismark job.
- **sRNA bigwig generation is surprisingly heavy** — 24 jobs totaling over 1 hour. Each sample produces bigwigs for 3 size classes (21/22/24 nt) × 2 strands, and merged replicates add more.
- **Downloads parallelize well** but are network-bound. ENA-first strategy helps; fasterq-dump fallback is slower.
- **All non-bismark jobs complete in under 17 minutes** individually.

## S. pombe — SLURM cluster (Elzar)

**Date**: 2026-03-12
**Environment**: Elzar SLURM cluster (bamdev2), up to 16 concurrent jobs
**Dataset**: 17 samples across 4 assay types (ChIP_broad, ChIP_narrow, RNAseq, sRNA)
**Reference**: S. pombe genome (~12 Mb)

| Metric | Value |
|--------|-------|
| Wall time | 54m 0s |
| Total jobs | 274 |
| Total CPU time | 49h 30m |
| Avg parallelism | 55.0x |

### Bottlenecks

| Rule | Total CPU | Max single job | Notes |
|------|-----------|----------------|-------|
| calling_peaks_macs2_pe | 43m | 3m 23s | 16 jobs — SLURM overhead dominates |
| make_bigwig_chip | 41m | 4m 43s | 11 jobs |
| make_srna_stranded_bigwigs | 37m | 3m 22s | 16 jobs |
| shortstack_map | 33m | 12m 51s | Longest single job |
| get_fastq_pe + get_fastq_se | 55m | 5m 11s | SRA downloads (network-bound) |

### Phase breakdown

| Phase | Total CPU | % of wall |
|-------|-----------|-----------|
| Other (misc analysis) | 3h 42m | 411% |
| Tracks (bigwigs) | 1h 40m | 186% |
| Download | 1h 37m | 180% |
| Peak calling | 1h 11m | 132% |
| QC | 1h 9m | 128% |
| Alignment | 1h 6m | 123% |
| sRNA analysis | 44m | 81% |

### Key observations

- No WGBS/mC samples (S. pombe lacks DNA methylation), so bismark bottleneck does not apply.
- **SLURM job overhead is significant for small jobs** — the cluster run takes 54 min vs 24 min local, despite higher CPU parallelism (55x vs 28x). Each job incurs ~1-2 min of scheduling + conda activation overhead, which dominates for jobs that run in seconds locally.
- **ShortStack is the longest single job** (12m 51s), followed by filter_bam_pe (8m 18s for WT_WCE_rep1).
- The pombe dataset is too small to benefit from cluster execution; local runs are ~2x faster.

## S. pombe — local execution (gemmule)

**Date**: 2026-03-12
**Environment**: gemmule (56 threads, local execution)
**Dataset**: 17 samples across 4 assay types (ChIP_broad, ChIP_narrow, RNAseq, sRNA)
**Reference**: S. pombe genome (~12 Mb)

| Metric | Value |
|--------|-------|
| Wall time | 24m 24s |
| Total jobs | 274 |
| Total CPU time | 11h 29m |
| Avg parallelism | 28.2x |

### Bottlenecks

| Rule | Total CPU | Max single job | Notes |
|------|-----------|----------------|-------|
| make_bigwig_chip | 13m 46s | 1m 23s | 11 jobs |
| making_stranded_matrix_on_targetfile | 13m 17s | 2m 21s | 6 deeptools matrix computations |
| get_fastq_pe + get_fastq_se | 22m 48s | 10m 22s | SRA downloads (network-bound) |
| filter_bam_pe | 9m 21s | 6m 3s (WT_WCE_rep1) | bowtie2 + samtools pipe |
| make_rna_stranded_bigwigs | 7m 2s | 3m 48s | Stranded RNA-seq bigwigs |

### Phase breakdown

| Phase | Total CPU | % of wall |
|-------|-----------|-----------|
| Other (misc analysis) | 35m | 142% |
| Download | 32m | 132% |
| Alignment | 26m | 106% |
| Tracks (bigwigs) | 22m | 89% |
| Peak calling | 8m | 34% |
| QC | 6m | 23% |
| sRNA analysis | 3m | 14% |

### Key observations

- No WGBS/mC samples (S. pombe lacks DNA methylation), so bismark bottleneck does not apply.
- **SRA download is the single slowest job** (10m 22s for WT_WCE_rep1), despite the small genome.
- **Parallelism (28.2x) is well below the 56 available threads** — the DAG is the limiting factor, not CPU. Many jobs are fast (<30s) and have sequential dependencies.
- All individual jobs complete in under 10 minutes. No single rule dominates the critical path.

## Generating profiles

```bash
# Profile the latest run (markdown to stdout)
python dev/profile_snakemake_log.py --latest

# Generate HTML report with Gantt timeline
python dev/profile_snakemake_log.py --latest --html dev/docs/profile_run.html

# Profile a specific log file
python dev/profile_snakemake_log.py .snakemake/log/2026-03-12T202344.snakemake.log

# Generate consolidated multi-section HTML from multiple runs
python dev/profile_snakemake_log.py --html dev/docs/profile_test_cases.html --multi \
  "H. sapiens hg38 chr21 (cluster)=.snakemake/log/2026-03-12T202344.534344.snakemake.log" \
  "S. pombe (cluster)=.snakemake/log/2026-03-12T230550.072095.snakemake.log" \
  "S. pombe (local)=.snakemake/log/2026-03-12T223119.439216.snakemake.log"
```

## Resource tiers

Current resource presets defined in `config/epicc-options.yaml`:

| Tier | Threads | Memory | Tmp | QOS | Used by |
|------|---------|--------|-----|-----|---------|
| `low` | 1 | 1 GB | 1 GB | default | Stats, dispatchers |
| `standard` | 4 | 2 GB | 2 GB | default | Index checks, light processing |
| `download` | 8 | 16 GB | 48 GB | slow_nice | SRA downloads, fastp |
| `heavy` | 4 | 16 GB | 48 GB | slow_nice | Peak calling, filtering, STAR index |
| `heavier` | 8 | 32 GB | 64 GB | slow_nice | BAM filtering, STAR mapping |
| `max` | 16 | 32 GB | 96 GB | slow_nice | Bismark, DMR calling, modBAM alignment |
| `single` | 1 | 32 GB | 48 GB | slow_nice | High-memory single-thread (genome stats) |

Note: Bismark rules compute `tmp_mb` dynamically from input FASTQ sizes (overriding the tier default).
