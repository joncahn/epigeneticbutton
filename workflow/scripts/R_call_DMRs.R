#!/usr/bin/env Rscript

library(dplyr)
library(DMRcaller)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)

threads<-as.numeric(args[1])
chromsizes<-read.table(args[2], col.names = c("chr", "length"))

# DMRcaller parallelizes its inner loop across the regions GRanges
# (one task per chromosome), so spawning more workers than chromosomes
# wastes memory without buying speed. Cap accordingly.
n_chrs <- nrow(chromsizes)
workers <- max(1L, min(threads, n_chrs))

# Use SnowParam (PSOCK) rather than MulticoreParam (fork). The latter
# triggered the reducer crash in
# github.com/joncahn/epigeneticbutton/issues/23: a forked worker hit
# OOM (COW pages of the methylation pool inflate per-worker RSS),
# parallel::mccollect returned no result for that slot, and
# BiocParallel's .reducer_add died on an empty index. PSOCK workers
# are independent R processes with a known baseline footprint and
# survive sibling deaths cleanly.
# OMP_NUM_THREADS=1 inside workers prevents BLAS/OpenMP from
# oversubscribing the slot (each worker inherits this on spawn).
if (workers > 1) {
    Sys.setenv(OMP_NUM_THREADS="1")
    register(SnowParam(workers=workers, type="SOCK"))
}
cat("R_call_DMRs.R parallelism setup: threads=", threads,
    " | n_chrs=", n_chrs, " | workers=", workers,
    " | BiocParallel backend=", class(bpparam())[1],
    "\n", sep="")
# args[3] is a comma-separated list of methylation contexts to call DMRs in,
# e.g. "CG,CHG,CHH" for plants or "CG" for animal genomes that lack
# substantive non-CpG methylation. Replaces the legacy "all"/"CG-only"
# string the script used to receive.
contexts<-strsplit(args[3], ",")[[1]]
sample1<-args[4]
sample2<-args[5]
nb_sample1<-as.numeric(args[6])
nb_sample2<-as.numeric(args[7])
output_dir<-args[8]
list_sample1<-args[9:(8+nb_sample1)]
list_sample2<-args[(9+nb_sample1):(8+nb_sample1+nb_sample2)]

chrs<-GRanges(seqnames = chromsizes$chr, ranges = IRanges(start = 1, end = chromsizes$length))

methylationDatasample1pool<-readBismarkPool(list_sample1)
methylationDatasample2pool<-readBismarkPool(list_sample2)

# Per-context DMRcaller parameters. CHG and CHH need different settings
# than CG (lower proportion-difference thresholds for the more diffuse
# asymmetric contexts; CHH uses bin-based instead of noise-filter mode).
context_params <- list(
  CG  = list(method="noise_filter", minProportionDifference=0.3),
  CHG = list(method="noise_filter", minProportionDifference=0.2),
  CHH = list(method="bins",          minProportionDifference=0.1)
)

call_dmrs_for_context <- function(ctx) {
  p <- context_params[[ctx]]
  base_args <- list(
    methylationDatasample1pool, methylationDatasample2pool,
    regions=chrs, context=ctx,
    method=p$method, binSize=200, test="score",
    pValueThreshold=0.01, minCytosinesCount=5,
    minProportionDifference=p$minProportionDifference,
    minGap=200, minSize=50, minReadsPerCytosine=3
  )
  # DMRcaller >=1.42 introduced a new `parallel=FALSE` default that
  # hardcodes BPPARAM to SerialParam regardless of cores= or any
  # globally-registered backend, so the legacy cores= parameter alone
  # is effectively ignored on those versions and the job runs serial
  # despite a multi-CPU allocation. Pass the explicit parallel/BPPARAM
  # pair when DMRcaller exposes them; fall back to cores= for older
  # versions (1.38 and earlier, which used parallel::mclapply directly).
  fnames <- names(formals(computeDMRs))
  if (workers > 1 && all(c("parallel", "BPPARAM") %in% fnames)) {
    base_args$parallel <- TRUE
    base_args$BPPARAM <- SnowParam(workers=workers, type="SOCK")
  }
  if ("cores" %in% fnames) {
    base_args$cores <- workers
  }
  do.call(computeDMRs, base_args)
}

summarize_dmrs <- function(dmrs, ctx) {
  count_col <- paste0(ctx, "_DMRs")
  if ( length(dmrs) > 0 ) {
    df <- data.frame(
      Chr=seqnames(dmrs), Start=start(dmrs)-1, End=end(dmrs),
      firstsample=mcols(dmrs)$proportion1, secondsample=mcols(dmrs)$proportion2,
      Pvalue=mcols(dmrs)$pValue
    ) %>% mutate(Delta=firstsample-secondsample)
    write.table(df,
      paste0(output_dir, "/mC/DMRs/", sample1, "__vs__", sample2, "__", ctx, "_DMRs.txt"),
      sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    summ <- mutate(df, Type=ifelse(Delta>0, "hyper", "hypo")) %>%
            group_by(Type) %>%
            summarize(!!count_col := n(), .groups = "drop")
  } else {
    summ <- tibble::tibble(Type=c("hyper", "hypo"), !!count_col := c(0, 0))
  }
  summ
}

# Iterate the requested contexts and merge per-context summaries on Type
# ("hyper"/"hypo"). Skipping CHG/CHH for animal genomes leaves the merged
# summary table with only the CG_DMRs column, which is exactly what we want.
summary_file <- NULL
for (ctx in contexts) {
  if (!ctx %in% names(context_params)) {
    warning(paste0("R_call_DMRs.R: unknown context '", ctx, "' — skipping"))
    next
  }
  dmrs <- call_dmrs_for_context(ctx)
  summ <- summarize_dmrs(dmrs, ctx)
  if (is.null(summary_file)) {
    summary_file <- summ
  } else {
    summary_file <- merge(summary_file, summ, by=c("Type"))
  }
}

# If no contexts were requested (shouldn't happen — Snakefile validates),
# emit an empty summary so the rule's output file still exists.
if (is.null(summary_file)) {
  summary_file <- tibble::tibble(Type=c("hyper", "hypo"))
}

summary_file<-mutate(summary_file, Sample=paste0(sample1,"_vs_",sample2)) %>% select(Sample, everything())
write.table(summary_file,paste0(output_dir,"/mC/DMRs/summary__",sample1,"__vs__",sample2,"__DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
