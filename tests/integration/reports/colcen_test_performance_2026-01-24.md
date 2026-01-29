# Pipeline Performance Report - ColCEN Integration Test

**Date:** January 24, 2026
**Config:** `tests/integration/data/test_config_colcen.yaml`
**Reference Genome:** ColCEN (full Arabidopsis genome with centromere assemblies)

## Overall Metrics

| Metric | Value |
|--------|-------|
| **Total Duration** | 5 hours 1 minute |
| **Total Jobs** | 237 |
| **Start Time** | 15:21:56 |
| **End Time** | 20:22:52 |
| **Status** | Completed successfully |

## Test Data

- **ChIP samples:** 12 (CenH3 IP + Input, 2 replicates each, 3 genotypes: WT, rdr126ddm1, rdr126ddm1hp5)
- **ONT samples:** 4 (pre-aligned modBAMs for rdr126ddm1, rdr126ddm1het, rdr126ddm1hethp5, rdr126ddm1hp5)
- **Reference genome:** ColCEN (~135 Mb)

## Timing by Pipeline Stage

| Stage | Jobs | Duration | Notes |
|-------|------|----------|-------|
| **ChIP-seq Processing** | 93 | 4h 42m | FASTQ → trimming → alignment → filtering → bigWig |
| **ChIP Peak Calling** | 31 | 3h 8m | MACS2 peak calling + IDR reproducibility analysis |
| **ONT Processing** | 28 | 1h 1m | modBAM validation → alignment → pileup → bigWig |
| **DMR Analysis** | 6 | 14m | DMRcaller pairwise comparisons (ONT vs ONT) |
| **Combined Analysis** | 54 | 3h 50m | Heatmaps, metaprofiles, genome browser tracks |
| **Reports** | 7 | 4h 3m | Mapping stats, peak stats, UpSet plots |

## Detailed Rule Timings

### ChIP-seq Rules

| Rule | Jobs | First Completed | Last Completed | Span |
|------|------|-----------------|----------------|------|
| `get_fastq_pe` | 12 | 15:21:56 | 15:24:07 | 2m 11s |
| `process_fastq_pe` | 12 | 15:24:47 | 15:27:28 | 2m 41s |
| `bowtie2_map_pe` | 12 | 16:23:15 | 19:43:14 | 3h 20m |
| `filter_chip_pe` | 12 | 16:28:18 | 19:47:46 | 3h 19m |
| `make_chip_stats_pe` | 12 | 16:28:57 | 19:49:07 | 3h 20m |
| `merging_chip_replicates` | 6 | 16:52:25 | 19:53:19 | 3h 1m |
| `make_bigwig_chip` | 9 | 18:06:22 | 20:03:42 | 1h 57m |
| `make_fingerprint_plot` | 6 | 18:07:03 | 19:57:50 | 1h 51m |

### ChIP Peak Calling Rules

| Rule | Jobs | First Completed | Last Completed | Span |
|------|------|-----------------|----------------|------|
| `making_pseudo_replicates` | 6 | 16:57:37 | 19:57:10 | 3h 0m |
| `calling_peaks_macs2_pe` | 15 | 18:05:42 | 20:02:12 | 1h 57m |
| `idr_analysis_replicates` | 3 | 19:37:42 | 19:56:30 | 19m |
| `best_peaks_pseudoreps` | 3 | 19:41:54 | 20:03:44 | 22m |
| `make_peak_stats` | 3 | 19:42:34 | 20:05:13 | 23m |
| `combine_peakfiles` | 1 | 20:05:13 | - | - |

### ONT Processing Rules

| Rule | Jobs | First Completed | Last Completed | Span |
|------|------|-----------------|----------------|------|
| `get_modbam` | 4 | 15:27:28 | 15:28:10 | 42s |
| `align_modbam` | 4 | 15:28:10 | 15:33:20 | 5m 10s |
| `modkit_summary` | 4 | 15:30:50 | 15:34:01 | 3m 11s |
| `modkit_pileup` | 4 | 15:48:14 | 15:58:07 | 9m 53s |
| `convert_bedmethyl_to_cx_report` | 4 | 15:49:44 | 16:00:18 | 10m 34s |
| `make_mc_stats_ont` | 4 | 15:53:26 | 16:03:29 | 10m 3s |
| `make_ont_bigwig_files` | 4 | 16:09:31 | 16:28:17 | 18m 46s |

### DMR Analysis Rules

| Rule | Jobs | First Completed | Last Completed | Span |
|------|------|-----------------|----------------|------|
| `call_DMRs_pairwise` | 6 | 16:02:49 | 16:17:03 | 14m 14s |

### Combined Analysis Rules

| Rule | Jobs | First Completed | Last Completed | Span |
|------|------|-----------------|----------------|------|
| `making_stranded_matrix_on_targetfile` | 12 | 16:33:18 | 20:10:47 | 3h 37m |
| `merging_matrix` | 6 | 16:37:10 | 20:12:48 | 3h 36m |
| `computing_matrix_scales` | 6 | 16:38:30 | 20:14:09 | 3h 36m |
| `plotting_profile_on_targetfile` | 6 | 16:39:51 | 20:15:40 | 3h 36m |
| `prep_browser_on_region` | 7 | 20:07:55 | 20:13:29 | 5m 34s |
| `make_single_loci_browser_plot` | 7 | 20:10:07 | 20:15:40 | 5m 33s |
| `plotting_heatmap_on_targetfile` | 3 | 20:15:40 | 20:15:40 | 0s |
| `sort_heatmap` | 3 | 20:20:51 | 20:21:31 | 40s |
| `plotting_sorted_heatmap_on_targetfile` | 3 | 20:22:11 | 20:22:52 | 41s |
| `merge_region_browser_plots` | 1 | 20:17:10 | - | - |

## Resource Usage (SLURM)

Based on a sample of jobs from sacct:

| Metric | Min | Max | Average |
|--------|-----|-----|---------|
| **Job Duration** | 1m 15s | 3m 2s | 2m 23s |
| **Memory (MaxRSS)** | 126 MB | 210 MB | 142 MB |

Memory distribution:
- 0-150 MB: 9 jobs (75%)
- 150-200 MB: 2 jobs (17%)
- 200-300 MB: 1 job (8%)

All jobs completed successfully with exit code 0.

## Bottlenecks Identified

1. **ChIP bowtie2 mapping** (3h 20m span)
   - 12 samples processed with some parallelism
   - Alignment to full ColCEN genome is computationally intensive

2. **Combined matrix generation** (3h 35m span)
   - Depends on ChIP processing completion
   - Large matrices for genome-wide analysis

3. **Peak calling with MACS2** (2h span)
   - 15 peak calling jobs (IP, pooled, pseudoreps)
   - IDR analysis adds additional time

## ONT Performance Notes

The ONT direct methylation processing completed relatively quickly (~1 hour total):
- modBAM validation and alignment: ~6 minutes
- modkit pileup: ~10 minutes per sample
- Context conversion to CX report: ~10 minutes per sample
- BigWig generation: ~19 minutes total

The ONT workflow was not a bottleneck in this test run. The pre-aligned modBAMs (19GB each) were processed efficiently using the updated validation script with streaming reads.

## Recommendations

1. Consider increasing parallelism for bowtie2 alignment if cluster resources allow
2. The ONT workflow could potentially run independently of ChIP in a production setting
3. Memory allocation could be reduced for most jobs (currently using ~150 MB average)

---

*Report generated from Snakemake log analysis*
*Log file: `.snakemake/log/2026-01-24T150502.884001.snakemake.log`*
