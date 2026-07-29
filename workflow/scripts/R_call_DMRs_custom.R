#!/usr/bin/env Rscript

# Custom parameter-sweep DMR caller (opt-in via custom_script_dmrs: true).
# Sweeps DMRcaller's method (noise_filter, bins) x binSize (100, 200, 500)
# for each requested methylation context, using the analytical thresholds
# from the options file (dmr_thresholds) rather than hardcoded values.

library(dplyr)
library(DMRcaller)

args = commandArgs(trailingOnly=TRUE)

threads<-as.numeric(args[1])
chromsizes<-read.table(args[2], col.names = c("chr", "length"))
contexts<-strsplit(args[3], ",")[[1]]            # e.g. "CG,CHG,CHH" or "CG"
sample1<-args[4]
sample2<-args[5]
nb_sample1<-as.numeric(args[6])
nb_sample2<-as.numeric(args[7])
output_dir<-args[8]
min_cytosines<-as.numeric(args[9])
q_value<-as.numeric(args[10])   # DMRcaller pValueThreshold (returned pValue is BH-adjusted)
min_gap<-as.numeric(args[11])
min_size<-as.numeric(args[12])
min_reads<-as.numeric(args[13])
min_diff_spec<-args[14]                           # "CG=0.3,CHG=0.2,CHH=0.1"
list_sample1<-args[15:(14+nb_sample1)]
list_sample2<-args[(15+nb_sample1):(14+nb_sample1+nb_sample2)]

# Parse the per-context minProportionDifference thresholds.
min_diff<-list()
for (kv in strsplit(min_diff_spec, ",")[[1]]) {
	parts<-strsplit(kv, "=")[[1]]
	min_diff[[parts[1]]]<-as.numeric(parts[2])
}

chrs<-GRanges(seqnames = chromsizes$chr, ranges = IRanges(start = 1, end = chromsizes$length))

methylationDatasample1pool<-readBismarkPool(list_sample1)
methylationDatasample2pool<-readBismarkPool(list_sample2)

tot_file <- NULL
for ( meth in c("noise_filter", "bins") ) {

	for ( bs in c(100, 200, 500) ) {

		summary_combo <- NULL
		for ( ctx in contexts ) {

			DMRs<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context=ctx, method=meth, binSize=bs, test="score", pValueThreshold=q_value, minCytosinesCount=min_cytosines, minProportionDifference=min_diff[[ctx]], minGap=min_gap, minSize=min_size, minReadsPerCytosine=min_reads, cores=threads)

			colname<-paste0(ctx, "_DMRs")
			if ( length(DMRs) > 0 ) {
				df<-data.frame(Chr=seqnames(DMRs),Start=start(DMRs)-1,End=end(DMRs),firstsample=mcols(DMRs)$proportion1,secondsample=mcols(DMRs)$proportion2, qValue=mcols(DMRs)$pValue) %>%
					mutate(Delta=firstsample-secondsample)

				write.table(df,paste0(output_dir,"/mC/DMRs/",sample1,"__vs__",sample2,"__",ctx,"_DMRs_",meth,"_",bs,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

				s<-mutate(df, Type=ifelse(Delta>0, "hyper", "hypo"), Method=meth, Binsize=bs) %>%
					group_by(Type, Method, Binsize) %>%
					summarize(n=n(), .groups = "drop")
				names(s)[names(s) == "n"]<-colname
			} else {
				s<-tibble::tibble(Type=c("hyper", "hypo"), Method=meth, Binsize=bs)
				s[[colname]]<-c(0, 0)
			}

			if (is.null(summary_combo)) {
				summary_combo<-s
			} else {
				summary_combo<-merge(summary_combo, s, by=c("Type", "Method", "Binsize"))
			}
		}
		tot_file <- bind_rows(tot_file, summary_combo)
	}
}

tot_file<-mutate(tot_file, Sample=paste0(sample1,"_vs_",sample2)) %>% select(Sample, everything())
write.table(tot_file,paste0(output_dir,"/mC/DMRs/summary__",sample1,"__vs__",sample2,"__DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
