#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DMRcaller))

args <- commandArgs(trailingOnly = TRUE)
context     <- args[1]
output_path <- args[2]
cx_path     <- args[3]

mdata <- readBismark(cx_path)
mdata_ctx <- mdata[mdata$context == context]
cat(sprintf("Cached %d positions for context %s from %s\n",
            length(mdata_ctx), context, basename(cx_path)))
saveRDS(mdata_ctx, output_path)
