# S. pombe Pipeline Profiling

**Date**: 2026-03-04
**Branch**: `big-refactor` @ `6cc1180`
**Options**: `tests/integration/data/test_options_pombe.yaml` (18 samples, 4 assays)
**Machine**: 112 cores available, run with `--cores 56`
**Wall clock**: 42 minutes (00:01 → 00:43)
**Steps completed**: 267 / 282 (95%)

## Errors (all data/network, not code bugs)

| Error | Count | Cause |
|-------|-------|-------|
| `get_fastq_pe` (WT_WCE_rep1) | 1 | Truncated ENA download (retried successfully) |
| `idr_analysis_replicates` | 2 | <20 peaks post-merge (H3K9me2/me3 in small genome) |
| `perform_pairwise_diff_peaks` | 2 | Cascading from IDR failures |

## Per-rule profiling (sorted by total CPU time)

| Rule | Jobs | Total(s) | Mean(s) | Max(s) | % total |
|------|------|----------|---------|--------|---------|
| **make_fingerprint_plot** | 8 | 6928 | 866 | 1005 | **45.1%** |
| get_fastq_pe | 8 | 1213 | 152 | 835 | 7.9% |
| making_stranded_matrix_on_targetfile | 6 | 878 | 146 | 154 | 5.7% |
| **filter_chip_pe** | 8 | 865 | 108 | 341 | 5.6% |
| **make_bigwig_chip** | 11 | 849 | 77 | 86 | 5.5% |
| get_fastq_se | 10 | 573 | 57 | 175 | 3.7% |
| process_fastq_pe | 8 | 359 | 45 | 109 | 2.3% |
| filter_rna_se | 4 | 335 | 84 | 99 | 2.2% |
| calling_peaks_macs2_pe | 16 | 315 | 20 | 27 | 2.1% |
| process_fastq_se | 10 | 311 | 31 | 61 | 2.0% |
| merging_rna_replicates | 2 | 303 | 152 | 168 | 2.0% |
| filter_chip_se | 3 | 289 | 96 | 161 | 1.9% |
| make_rna_stranded_bigwigs | 6 | 261 | 44 | 80 | 1.7% |
| STAR_map_se | 4 | 210 | 53 | 58 | 1.4% |
| making_pseudo_replicates | 8 | 191 | 24 | 64 | 1.2% |
| analyze_all_srna_samples_on_target_file | 2 | 172 | 86 | 98 | 1.1% |
| calling_peaks_macs2_se | 5 | 172 | 34 | 46 | 1.1% |
| shortstack_map | 3 | 170 | 57 | 79 | 1.1% |
| make_single_loci_browser_plot | 6 | 157 | 26 | 27 | 1.0% |
| merging_matrix | 3 | 153 | 51 | 55 | 1.0% |
| (41 other rules) | — | 716 | — | — | 4.7% |
| **Total** | | **15344** | | | |

## Key observations

### 1. plotFingerprint is the dominant bottleneck (45% of CPU time)

deeptools `plotFingerprint` reads entire BAM files to compute genome-wide coverage
distributions. Each job averaged **14.4 minutes** on the small S. pombe genome
(12 Mb). On larger genomes this would be far worse.

All 8 jobs ran concurrently in the final phase of the pipeline and were the last
thing to finish — they extended wall clock time by ~15 minutes after all other
rules had completed.

### 2. filter_chip_pe shows bimodal timing

Two WCE (whole-cell extract) samples took ~340s each while the six IP samples
averaged ~33s. WCE files are larger (more reads, no enrichment-based filtering
by MACS2 upstream), which explains the disparity. The chromap alignment step
itself is fast; the time is dominated by samtools sort + markdup on the larger
BAM files.

### 3. Downloads are variable but not a code bottleneck

One PE download took 835s (likely ENA server load or network congestion), but
most completed in <20s. The `gzip -t` validation added in commit `6cc1180`
catches truncated downloads; the `retries: 3` directive handles transient
failures.

### 4. deeptools matrix/bigwig rules are moderately expensive

`making_stranded_matrix_on_targetfile` (878s total) and `make_bigwig_chip`
(849s total) together account for 11.2% of CPU time. These are inherent to the
deeptools compute-then-plot workflow.

### 5. RNA and sRNA pipelines are fast

STAR mapping (210s total for 4 samples) and ShortStack (170s for 3 samples)
are well-optimized. RNA merging is slow (303s for 2 jobs) due to large BAM
sorting after merge.

## Optimization opportunities

### High impact

- **plotFingerprint**: Consider making it optional via an options flag
  (`chip_fingerprint_plots: true/false`), or reducing its `--numberOfSamples`
  parameter (default 500,000) to speed up computation at the cost of minor
  precision loss. Could also lower its thread allocation to allow more
  concurrent jobs.

- **Duplicate SRA downloads**: `WT_WCE_rep1` and `dcr1_WCE_rep1` share the
  same accession (SRR5445712). The pipeline downloads it twice independently.
  A shared download rule or symlink mechanism would save time and bandwidth.

### Medium impact

- **IDR minimum peak threshold**: The S. pombe test case consistently fails IDR
  with <20 peaks post-merge. The IDR rule could catch this gracefully (warn +
  create empty output) instead of erroring, avoiding cascading diff_peaks
  failures.

- **samtools sort memory**: The sort step in filter_chip_pe/se uses default
  memory. Adding `-m` (memory per thread) could reduce I/O pressure on large
  BAMs.

### Low impact

- **making_stranded_matrix_on_targetfile**: Already parallelized via deeptools.
  Limited room for improvement without changing the analysis.

- **process_fastq_pe/se (fastp)**: Already fast (31-45s mean). fastp is
  well-optimized.
