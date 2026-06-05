#!/usr/bin/env Rscript
# Estimate the noise-floor sigma for auto-calibrated DMR min_diff threshold.
#
# For each position with data in both groups, computes the between-group mean
# methylation difference.  The vast majority of positions carry no real DMR,
# so the SD of those differences is dominated by measurement noise, which
# naturally scales with the context's baseline methylation:
#   CG  (p ~ 0.8) → p*(1-p)/N large → σ large → higher threshold
#   CHH (p ~ 0.02) → p*(1-p)/N small → σ small → lower threshold
# This avoids hard-coded per-context or per-organism values.
#
# Threshold = sigma_n * sigma (default sigma_n=3, ≈ 3σ above noise floor).
# A floor of 0.05 prevents pathologically low thresholds when coverage is
# very high or methylation is near zero everywhere.
#
# Output: one-line JSON consumed by the Snakemake params lambda in
# call_DMRs_for_pair_context when dmr_thresholds.min_diff: "auto".

suppressPackageStartupMessages(library(GenomicRanges))

args        <- commandArgs(trailingOnly = TRUE)
context     <- args[1]
json_out    <- args[2]
sigma_n     <- as.numeric(args[3])
n1          <- as.integer(args[4])
reps1_paths <- args[5:(4 + n1)]
reps2_paths <- args[(5 + n1):length(args)]
n2          <- length(reps2_paths)

cat(sprintf("R_calibrate_dmr_min_diff.R: context=%s n1=%d n2=%d sigma_n=%.1f\n",
            context, n1, n2, sigma_n))

reps1 <- lapply(reps1_paths, readRDS)
reps2 <- lapply(reps2_paths, readRDS)

# Compute per-position mean methylation rate for one group.
# Uses the same packed-key union as R_call_DMRs_pair_metilene.R so dmC
# inputs (where replicates may cover different position sets) work correctly.
group_mean_rates <- function(reps) {
    all_chr <- unlist(lapply(reps, function(gr) as.character(seqnames(gr))), use.names = FALSE)
    all_pos <- unlist(lapply(reps, function(gr) start(gr)), use.names = FALSE)
    chr_levels <- unique(all_chr)
    keys <- match(all_chr, chr_levels) * 1e10 + all_pos
    keep <- !duplicated(keys)
    union_keys <- keys[keep]
    ord <- order(match(all_chr[keep], chr_levels), all_pos[keep])
    union_keys <- union_keys[ord]
    rates <- lapply(reps, function(gr) {
        gkeys <- match(as.character(seqnames(gr)), chr_levels) * 1e10 + start(gr)
        m <- mcols(gr)
        r <- m$readsM / m$readsN
        r[m$readsN == 0] <- NA
        r[match(union_keys, gkeys)]
    })
    list(keys = union_keys, mean = rowMeans(do.call(cbind, rates), na.rm = TRUE))
}

g1 <- group_mean_rates(reps1)
g2 <- group_mean_rates(reps2)

# Differences at positions present in both groups.
idx_in_g1 <- match(g2$keys, g1$keys)
shared    <- !is.na(idx_in_g1)
diffs     <- g1$mean[idx_in_g1[shared]] - g2$mean[shared]
diffs     <- diffs[!is.na(diffs)]
n_pos     <- length(diffs)

# Context defaults used as fallback when there are too few shared positions.
context_defaults <- c(CG = 0.3, CHG = 0.2, CHH = 0.1)

if (n_pos < 100) {
    warning(sprintf(
        "Only %d shared positions for sigma estimation in context %s; falling back to default %.1f",
        n_pos, context, context_defaults[[context]]
    ))
    threshold <- unname(context_defaults[[context]])
    sigma     <- NA_real_
} else {
    sigma     <- sd(diffs)
    threshold <- max(0.05, sigma_n * sigma)
    cat(sprintf("  n_positions=%d  sigma=%.4f  threshold=%.4f (%.1f * sigma)\n",
                n_pos, sigma, threshold, sigma_n))
}

sigma_str <- if (is.na(sigma)) "null" else sprintf("%.6f", sigma)
writeLines(
    sprintf('{"context":"%s","min_diff":%.6f,"sigma":%s,"sigma_n":%.1f,"n_positions":%d}',
            context, threshold, sigma_str, sigma_n, n_pos),
    json_out
)
