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
# G is metilene's --maxseg flag: caps the maximum continuous segment
# length to bound per-segment memory. Default -1 (no cap) works for CG
# and CHG. CHH on dense-coverage long chromosomes (e.g. ColCEN Chr1)
# produces single segments of ~8 million CpGs that drive metilene RSS to
# 400-700 GB and OOM-kill the job, so we cap CHH. The downside of
# capping is that DMRs straddling a chunk boundary may get split or
# (if the pieces fall below -m) missed. With -M 300, the per-DMR split
# risk is roughly M/G — 3% at G=10000. Larger G is also slightly faster
# (per-chunk overhead dominates per EpiDiverse), so the only cost of
# going larger is RSS. -G 10000 measured ~33 GB / 50 min on ColCEN CHH
# dmC pairs. -G 50000 fit those dmC pairs but OOM-killed PBAT/WGBS CHH
# pairs at 200-300 GB (denser CX input than dmC), so we stay at 10000.
metilene_params <- list(
  CG  = list(d = 0.3, m = 5L, M = 300L, v = 0.7, G = -1L),
  CHG = list(d = 0.2, m = 5L, M = 300L, v = 0.7, G = -1L),
  CHH = list(d = 0.1, m = 5L, M = 300L, v = 0.3, G = 10000L)
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

# Prominent warning when either group has fewer than 2 replicates.
# Metilene's MWU and 2D-KS tests assume biological replication; with
# n=1 in either group there is no variance estimate and the reported
# q-values are NOT meaningful significance values (effectively a sign
# test on the single observation).
if (n1 < 2 || n2 < 2) {
  bar <- strrep("=", 78)
  msg <- paste(
    "", bar,
    sprintf("WARNING: metilene called with insufficient replicates (n1=%d, n2=%d)", n1, n2),
    bar,
    "  Metilene's MWU and 2D-KS tests assume biological replication.",
    "  With n=1 in either group there is no variance estimate, so the",
    "  reported q-values are NOT meaningful significance values --",
    "  treat the output as candidate DMR coordinates only, not as",
    "  statistically supported calls.",
    "",
    "  For meaningful q-values, rerun with >=2 biological replicates",
    "  in each group.",
    bar, "",
    sep = "\n"
  )
  cat(msg, file = stderr())
}

# Compute the union of (chr, pos) across all replicates.
#
# Bismark-aligned samples have CX reports covering every genomic position
# of the context, so all per-replicate caches have identical positions.
# dmC samples come from `modkit pileup → bedMethyl → CX`, where bedMethyl
# only contains positions with coverage — so different replicates have
# different position sets. Take the union and pad missing positions with
# NA, which metilene reads as "." (missing).
#
# Keep chr and pos as parallel vectors and dedupe via a packed numeric
# key (chr_idx * 1e10 + pos). The earlier paste/strsplit/do.call(rbind)
# pattern aborted R with SIGBUS on ~35M-key CHH inputs because rbind on
# a list of millions of length-2 char vectors blows up internally. 1e10
# is well above any plant chromosome length (<2e8) and well under 2^53,
# so the packed keys remain exact in double precision.
chr_list <- lapply(all_reps, function(gr) as.character(seqnames(gr)))
pos_list <- lapply(all_reps, function(gr) start(gr))
all_chr <- unlist(chr_list, use.names = FALSE)
all_pos <- unlist(pos_list, use.names = FALSE)
chr_levels <- unique(all_chr)
all_keys_num <- match(all_chr, chr_levels) * 1e10 + all_pos
keep <- !duplicated(all_keys_num)
union_chr <- all_chr[keep]
union_pos <- all_pos[keep]
union_keys_num <- all_keys_num[keep]
ord <- order(union_chr, union_pos)
union_chr <- union_chr[ord]
union_pos <- union_pos[ord]
union_keys_num <- union_keys_num[ord]
n_pos <- length(union_chr)

# Per-replicate methylation rates aligned to the union of positions.
rates <- lapply(all_reps, function(gr) {
  rep_keys_num <- match(as.character(seqnames(gr)), chr_levels) * 1e10 + start(gr)
  m <- mcols(gr)
  rep_rate <- m$readsM / m$readsN
  rep_rate[m$readsN == 0] <- NA
  rep_rate[match(union_keys_num, rep_keys_num)]
})

# ---------------------------------------------------------------------------
# Write the metilene input TSV
# ---------------------------------------------------------------------------
sample_names <- c(paste0("g1_v", seq_len(n1)), paste0("g2_v", seq_len(n2)))

tmpdir <- Sys.getenv("TMPDIR", "/tmp")
pid <- Sys.getpid()
input_tsv  <- file.path(tmpdir, sprintf("metilene_input_%s_%d.tsv",  context, pid))
output_tsv <- file.path(tmpdir, sprintf("metilene_output_%s_%d.tsv", context, pid))

# Build a data.frame keyed (chr, pos) with one column per replicate.
input_df <- data.frame(
  chr = union_chr,
  pos = union_pos,
  stringsAsFactors = FALSE
)
for (i in seq_along(rates)) {
  rate_col <- rates[[i]]
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
  # -G only when capping segment length. Passing `-G -1` as a "no cap"
  # sentinel breaks metilene's argparser (it reads `-1` as another flag);
  # metilene's internal default for --maxseg is already -1 (no cap), so
  # omitting -G is the correct way to express "no cap".
  if (p$G > 0L) c("-G", as.character(p$G)) else character(0),
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
