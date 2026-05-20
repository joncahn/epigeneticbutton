#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DMRcaller)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
threads         <- as.numeric(args[1])
context         <- args[2]
chromsizes_path <- args[3]
dmrs_out        <- args[4]
counts_out      <- args[5]
n1              <- as.integer(args[6])
reps1_paths     <- args[7:(6 + n1)]
reps2_paths     <- args[(7 + n1):length(args)]
n2 <- length(reps2_paths)

chromsizes <- read.table(chromsizes_path, col.names = c("chr", "length"))
chrs <- GRanges(seqnames = chromsizes$chr,
                ranges   = IRanges(start = 1, end = chromsizes$length))

# DMRcaller parallelizes its inner loop one task per chromosome; capping
# workers at min(threads, n_chrs) avoids forking empty workers. SnowParam
# (PSOCK) over MulticoreParam (fork) survives worker OOM cleanly — see
# github.com/joncahn/epigeneticbutton/issues/23.
n_chrs <- nrow(chromsizes)
workers <- max(1L, min(threads, n_chrs))
if (workers > 1) {
    Sys.setenv(OMP_NUM_THREADS = "1")
    register(SnowParam(workers = workers, type = "SOCK"))
}

reps1 <- lapply(reps1_paths, readRDS)
reps2 <- lapply(reps2_paths, readRDS)
cat(sprintf("R_call_DMRs_pair.R: n1=%d n2=%d context=%s threads=%d workers=%d backend=%s\n",
            n1, n2, context, threads, workers, class(bpparam())[1]))

# Replicate-aware path uses computeDMRsReplicates with betareg test, which
# only accepts methods {neighbourhood, bins}. We use bins for all three
# contexts: neighbourhood + betareg is computationally prohibitive on
# whole-plant-genome 2-vs-2 replicate data (ColCEN integration test:
# >24 h per CG/CHG pair vs ~hours per CHH/bins pair). Bins gives a
# fixed-window segmentation but preserves the replicate-variance
# modeling that's the whole point of computeDMRsReplicates. Pooled
# fallback for N<2 uses the original noise_filter (CG/CHG) / bins (CHH)
# pair with the score test that the pre-refactor pipeline was tuned for.
context_params_replicates <- list(
  CG  = list(method = "bins", minProportionDifference = 0.3),
  CHG = list(method = "bins", minProportionDifference = 0.2),
  CHH = list(method = "bins", minProportionDifference = 0.1)
)
context_params_pooled <- list(
  CG  = list(method = "noise_filter", minProportionDifference = 0.3),
  CHG = list(method = "noise_filter", minProportionDifference = 0.2),
  CHH = list(method = "bins",          minProportionDifference = 0.1)
)

use_replicates <- (n1 >= 2 && n2 >= 2)

if (use_replicates) {
    p <- context_params_replicates[[context]]
    joined <- Reduce(joinReplicates, c(reps1, reps2))
    condition <- c(rep("sample1", n1), rep("sample2", n2))

    base_args <- list(
      methylationData = joined,
      condition       = condition,
      regions         = chrs,
      context         = context,
      method          = p$method,
      binSize         = 200,
      test            = "betareg",
      pValueThreshold = 0.01,
      minCytosinesCount = 5,
      minProportionDifference = p$minProportionDifference,
      minGap          = 200,
      minSize         = 50,
      minReadsPerCytosine = 3
    )
    # parallel=FALSE: computeDMRsReplicates's SnowParam path has a
    # reducer bug (DMRcaller 1.42) that masks worker R errors as
    # "wrong args for environment subassignment" and crashes the job
    # even when the underlying betareg fits would have succeeded. A
    # serial run of the same input (verified on the
    # PBAT__Col0_seedling vs WGBS__Col0_seedling CHG pair, 2026-05-18)
    # completes in ~60 min and returns the expected DMRs. Trade
    # parallelism over chromosomes for correctness here; the bins
    # method is fast enough sequentially that the wall-time cost is
    # acceptable for ColCEN-scale inputs.
    fnames <- names(formals(computeDMRsReplicates))
    if ("parallel" %in% fnames) base_args$parallel <- FALSE
    if ("cores"    %in% fnames) base_args$cores    <- 1L
    dmrs <- do.call(computeDMRsReplicates, base_args)
} else {
    warning(sprintf("Pair has insufficient replicates (n1=%d, n2=%d) for the betareg replicate model; falling back to pooled computeDMRs",
                    n1, n2))
    p <- context_params_pooled[[context]]
    pool1 <- if (length(reps1) > 1) Reduce(poolTwoMethylationDatasets, reps1) else reps1[[1]]
    pool2 <- if (length(reps2) > 1) Reduce(poolTwoMethylationDatasets, reps2) else reps2[[1]]

    base_args <- list(
      pool1, pool2,
      regions = chrs, context = context,
      method = p$method, binSize = 200, test = "score",
      pValueThreshold = 0.01, minCytosinesCount = 5,
      minProportionDifference = p$minProportionDifference,
      minGap = 200, minSize = 50, minReadsPerCytosine = 3
    )
    fnames <- names(formals(computeDMRs))
    if (workers > 1 && all(c("parallel", "BPPARAM") %in% fnames)) {
        base_args$parallel <- TRUE
        base_args$BPPARAM  <- SnowParam(workers = workers, type = "SOCK")
    }
    if ("cores" %in% fnames) base_args$cores <- workers
    dmrs <- do.call(computeDMRs, base_args)
}

# Both functions return GRanges with proportion1, proportion2, pValue —
# the columns we use downstream. computeDMRsReplicates also adds
# direction/regionType, which we ignore (Delta sign derives Type below).
if (length(dmrs) > 0) {
    df <- data.frame(
      Chr          = seqnames(dmrs),
      Start        = start(dmrs) - 1,
      End          = end(dmrs),
      firstsample  = mcols(dmrs)$proportion1,
      secondsample = mcols(dmrs)$proportion2,
      Pvalue       = mcols(dmrs)$pValue
    )
    df$Delta <- df$firstsample - df$secondsample
    write.table(df, dmrs_out, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    n_hyper <- sum(df$Delta > 0)
    n_hypo  <- sum(df$Delta <= 0)
} else {
    write.table(data.frame(Chr = character(), Start = integer(), End = integer(),
                           firstsample = numeric(), secondsample = numeric(),
                           Pvalue = numeric(), Delta = numeric()),
                dmrs_out, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    n_hyper <- 0L; n_hypo <- 0L
}

counts <- data.frame(Type = c("hyper", "hypo"))
counts[[paste0(context, "_DMRs")]] <- c(n_hyper, n_hypo)
write.table(counts, counts_out, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
