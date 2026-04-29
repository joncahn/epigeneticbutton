#!/usr/bin/env Rscript

library(dplyr)
library(DMRcaller)
library(BiocParallel)

args = commandArgs(trailingOnly=TRUE)

threads<-as.numeric(args[1])
chromsizes<-read.table(args[2], col.names = c("chr", "length"))

# Wire the allocated threads through every parallel mechanism we might
# hit. DMRcaller 1.38 uses parallel::mclapply with an explicit
# mc.cores=cores arg, so passing cores=threads to computeDMRs() below
# already covers the primary path. The setup below adds belt-and-
# suspenders for the remaining vectors:
#   - register(MulticoreParam(workers=threads)) makes bpparam() return
#     MulticoreParam instead of SerialParam — fixes the explicit
#     SerialParam observation from
#     github.com/joncahn/epigeneticbutton/issues/23 and future-proofs
#     against DMRcaller versions that switch to BiocParallel.
#   - options(mc.cores) is the default for any mclapply() that doesn't
#     set mc.cores explicitly.
#   - OMP_NUM_THREADS prevents any C/Fortran OpenMP region (e.g. in
#     BLAS-backed calls) from oversubscribing the slot.
if (threads > 1) {
    register(MulticoreParam(workers=threads))
    options(mc.cores=threads)
    Sys.setenv(OMP_NUM_THREADS=as.character(threads))
}
cat("R_call_DMRs.R parallelism setup: threads=", threads,
    " | BiocParallel backend=", class(bpparam())[1],
    " | options(mc.cores)=", getOption("mc.cores", 1L),
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
  computeDMRs(
    methylationDatasample1pool, methylationDatasample2pool,
    regions=chrs, context=ctx,
    method=p$method, binSize=200, test="score",
    pValueThreshold=0.01, minCytosinesCount=5,
    minProportionDifference=p$minProportionDifference,
    minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads
  )
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
