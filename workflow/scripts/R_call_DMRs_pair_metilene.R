#!/usr/bin/env Rscript
# Metilene-backed driver for the call_DMRs_for_pair_context rule.
#
# Reads the per-(replicate, context) RDS caches produced by
# cache_mc_replicate_for_context, writes a metilene-format TSV (chr,
# pos, per-replicate methylation rate columns), calls the metilene
# binary, parses its BED-like output, and writes our standard DMR
# table + per-context hyper/hypo counts file.
#
# Per-context parameters are plant-tuned: -d (mean methylation
# difference threshold) matches DMRcaller's minProportionDifference
# (CG=0.3, CHG=0.2, CHH=0.1); -v (per-position valley filter) is
# lowered for CHH because CHH per-site signal is intrinsically noisy
# at the low-methylation baseline plants typically have there.

suppressPackageStartupMessages({
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
threads         <- as.integer(args[1])
context         <- args[2]
chromsizes_path <- args[3]  # unused (kept for arg-order parity with DMRcaller path)
dmrs_out        <- args[4]
counts_out      <- args[5]
n1              <- as.integer(args[6])
reps1_paths     <- args[7:(6 + n1)]
reps2_paths     <- args[(7 + n1):length(args)]
n2 <- length(reps2_paths)

# ---------------------------------------------------------------------------
# Per-context metilene parameters
# ---------------------------------------------------------------------------
metilene_params <- list(
  CG  = list(d = 0.3, m = 5L, M = 300L, v = 0.7),
  CHG = list(d = 0.2, m = 5L, M = 300L, v = 0.7),
  CHH = list(d = 0.1, m = 5L, M = 300L, v = 0.3)
)
p <- metilene_params[[context]]
if (is.null(p)) stop(sprintf("Unknown context '%s'", context))

# ---------------------------------------------------------------------------
# Load per-replicate RDS caches and compute methylation rate per position
# ---------------------------------------------------------------------------
reps1 <- lapply(reps1_paths, readRDS)
reps2 <- lapply(reps2_paths, readRDS)
all_reps <- c(reps1, reps2)

cat(sprintf("R_call_DMRs_pair_metilene.R: n1=%d n2=%d context=%s threads=%d\n",
            n1, n2, context, threads))

rate_at_position <- function(gr) {
  m <- mcols(gr)
  rate <- m$readsM / m$readsN
  rate[m$readsN == 0] <- NA  # NA when no coverage; metilene treats "." as missing
  rate
}

rates <- lapply(all_reps, rate_at_position)

# CX reports for a given genome are deterministic on positions, so all
# replicates have the same chr/pos order. Use rep 1 as the template.
template <- all_reps[[1]]
n_pos <- length(template)
stopifnot(all(sapply(all_reps, length) == n_pos))

# ---------------------------------------------------------------------------
# Write the metilene input TSV
# ---------------------------------------------------------------------------
sample_names <- c(paste0("g1_v", seq_len(n1)), paste0("g2_v", seq_len(n2)))

tmpdir <- Sys.getenv("TMPDIR", "/tmp")
pid <- Sys.getpid()
input_tsv  <- file.path(tmpdir, sprintf("metilene_input_%s_%d.tsv",  context, pid))
output_tsv <- file.path(tmpdir, sprintf("metilene_output_%s_%d.tsv", context, pid))

# Build a data.frame keyed (chr, pos) with one column per replicate.
# Sort by chr then pos because metilene requires sorted input.
ord <- order(as.character(seqnames(template)), start(template))
input_df <- data.frame(
  chr = as.character(seqnames(template))[ord],
  pos = start(template)[ord],
  stringsAsFactors = FALSE
)
for (i in seq_along(rates)) {
  rate_col <- rates[[i]][ord]
  rate_col_chr <- sprintf("%.4f", rate_col)
  rate_col_chr[is.na(rate_col)] <- "."
  input_df[[sample_names[i]]] <- rate_col_chr
}

# Header line (metilene -h: tab-separated; "chr\tpos\t<sample1>\t<sample2>\t...")
con <- file(input_tsv, "w")
writeLines(paste(c("chr", "pos", sample_names), collapse = "\t"), con)
write.table(input_df, con, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
close(con)
cat(sprintf("Wrote metilene input: %s (%d positions, %d samples)\n",
            input_tsv, n_pos, length(sample_names)))

# ---------------------------------------------------------------------------
# Run metilene
# ---------------------------------------------------------------------------
metilene_args <- c(
  "-a", "g1",
  "-b", "g2",
  "-d", as.character(p$d),
  "-m", as.character(p$m),
  "-M", as.character(p$M),
  "-v", as.character(p$v),
  "-t", as.character(max(1L, threads)),
  "-X", "1",
  "-Y", "1",
  input_tsv
)
cat("Running: metilene", paste(metilene_args, collapse = " "), "\n")
rc <- system2("metilene", args = metilene_args, stdout = output_tsv)
cat(sprintf("metilene exit code: %d\n", rc))
if (rc != 0) stop(sprintf("metilene failed (exit %d)", rc))

# ---------------------------------------------------------------------------
# Parse metilene output and translate to our standard DMR table
# ---------------------------------------------------------------------------
# metilene de-novo (-f 1) output columns:
#   chr  start  end  q-value  meth-diff  #cpgs  pMWU  p2D-KS  mean-g1  mean-g2
# No header. start is 0-based, end is BED-style exclusive.
out_cols <- c("Chr", "Start", "End", "qValue", "meth_diff",
              "n_cpgs", "pMWU", "p2D_KS", "mean_g1", "mean_g2")
res <- tryCatch(
  read.table(output_tsv, sep = "\t", header = FALSE,
             col.names = out_cols, fill = TRUE,
             stringsAsFactors = FALSE),
  error = function(e) {
    cat(sprintf("Note: metilene output empty or unparseable (%s)\n",
                conditionMessage(e)))
    data.frame()
  }
)

if (nrow(res) > 0) {
  df_out <- data.frame(
    Chr          = res$Chr,
    Start        = res$Start,
    End          = res$End,
    firstsample  = res$mean_g1,
    secondsample = res$mean_g2,
    Pvalue       = res$qValue,
    Delta        = res$mean_g1 - res$mean_g2
  )
  write.table(df_out, dmrs_out, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  n_hyper <- sum(df_out$Delta > 0)
  n_hypo  <- sum(df_out$Delta <= 0)
  cat(sprintf("Wrote %d DMRs (hyper=%d, hypo=%d) to %s\n",
              nrow(df_out), n_hyper, n_hypo, dmrs_out))
} else {
  write.table(data.frame(Chr = character(), Start = integer(), End = integer(),
                         firstsample = numeric(), secondsample = numeric(),
                         Pvalue = numeric(), Delta = numeric()),
              dmrs_out, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  n_hyper <- 0L; n_hypo <- 0L
  cat(sprintf("No DMRs called; wrote empty header to %s\n", dmrs_out))
}

counts <- data.frame(Type = c("hyper", "hypo"))
counts[[paste0(context, "_DMRs")]] <- c(n_hyper, n_hypo)
write.table(counts, counts_out, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Cleanup TMPDIR scratch
suppressWarnings(file.remove(input_tsv, output_tsv))
