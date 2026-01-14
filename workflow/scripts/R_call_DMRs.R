#!/usr/bin/env Rscript

library(dplyr)
library(DMRcaller)

args = commandArgs(trailingOnly=TRUE)

threads<-as.numeric(args[1])
chromsizes<-read.table(args[2], col.names = c("chr", "length"))
context<-args[3]
sample1<-args[4]
sample2<-args[5]
nb_sample1<-as.numeric(args[6])
nb_sample2<-as.numeric(args[7])
list_sample1<-args[8:(7+nb_sample1)]
list_sample2<-args[(8+nb_sample1):(7+nb_sample1+nb_sample2)]

chrs<-GRanges(seqnames = chromsizes$chr, ranges = IRanges(start = 1, end = chromsizes$length))

methylationDatasample1pool<-readBismarkPool(list_sample1)
methylationDatasample2pool<-readBismarkPool(list_sample2)

DMRsCGpool<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context="CG", method="noise_filter", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads)
if ( length(DMRsCGpool) > 0 ) {
	CGpool<-data.frame(Chr=seqnames(DMRsCGpool),Start=start(DMRsCGpool)-1,End=end(DMRsCGpool),firstsample=mcols(DMRsCGpool)$proportion1,secondsample=mcols(DMRsCGpool)$proportion2, Pvalue=mcols(DMRsCGpool)$pValue) %>%
			mutate(Delta=firstsample-secondsample)

	write.table(CGpool,paste0("results/mC/DMRs/",sample1,"__vs__",sample2,"__CG_DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

	summary_file<-mutate(CGpool, Type=ifelse(Delta>0, "hyper", "hypo")) %>%
				group_by(Type) %>%
				summarize(CG_DMRs=n(), .groups = "drop")
} else {
	summary_file<-tibble::tibble(Type=c("hyper", "hypo"), CG_DMRs=c(0, 0))
}

if (context == "all") {
	DMRsCHGpool<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context="CHG", method="noise_filter", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.2, minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads)
	if ( length(DMRsCHGpool) > 0 ) {
		CHGpool<-data.frame(Chr=seqnames(DMRsCHGpool),Start=start(DMRsCHGpool)-1,End=end(DMRsCHGpool),firstsample=mcols(DMRsCHGpool)$proportion1,secondsample=mcols(DMRsCHGpool)$proportion2, Pvalue=mcols(DMRsCHGpool)$pValue) %>%
				mutate(Delta=firstsample-secondsample)
		write.table(CHGpool,paste0("results/mC/DMRs/",sample1,"__vs__",sample2,"__CHG_DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	
		summary_fileCHG<-mutate(CHGpool, Type=ifelse(Delta>0, "hyper", "hypo")) %>%
						group_by(Type) %>%
						summarize(CHG_DMRs=n(), .groups = "drop")
	} else {
		summary_fileCHG<-tibble::tibble(Type=c("hyper", "hypo"), CHG_DMRs=c(0, 0))
	}
	
	DMRsCHHpool<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context="CHH", method="bins", binSize=200, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.1, minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads)
	if ( length(DMRsCHHpool) > 0 ) {
		CHHpool<-data.frame(Chr=seqnames(DMRsCHHpool),Start=start(DMRsCHHpool)-1,End=end(DMRsCHHpool),firstsample=mcols(DMRsCHHpool)$proportion1,secondsample=mcols(DMRsCHHpool)$proportion2, Pvalue=mcols(DMRsCHHpool)$pValue) %>%
				mutate(Delta=firstsample-secondsample)
		write.table(CHHpool,paste0("results/mC/DMRs/",sample1,"__vs__",sample2,"__CHH_DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	
		summary_fileCHH<-mutate(CHHpool, Type=ifelse(Delta>0, "hyper", "hypo")) %>%
						group_by(Type) %>%
						summarize(CHH_DMRs=n(), .groups = "drop")
	} else {
		summary_fileCHH<-tibble::tibble(Type=c("hyper", "hypo"), CHH_DMRs=c(0, 0))
	}					
	
	summary_file<-merge(summary_file, summary_fileCHG, by=c("Type"))
	summary_file<-merge(summary_file, summary_fileCHH, by=c("Type"))
}

summary_file<-mutate(summary_file, Sample=paste0(sample1,"_vs_",sample2)) %>% select(Sample, everything())
write.table(summary_file,paste0("results/mC/DMRs/summary__",sample1,"__vs__",sample2,"__DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
