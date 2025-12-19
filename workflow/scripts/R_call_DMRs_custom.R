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

tot_file <- NULL
for ( meth in c("noise_filter", "bins")) {
	
	for ( bs in c(100, 200, 500) ) {
		
		DMRsCGpool<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context="CG", method=meth, binSize=bs, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads)

		if ( length(DMRsCGpool) > 0 ) {
			CGpool<-data.frame(Chr=seqnames(DMRsCGpool),Start=start(DMRsCGpool)-1,End=end(DMRsCGpool),firstsample=mcols(DMRsCGpool)$proportion1,secondsample=mcols(DMRsCGpool)$proportion2, Pvalue=mcols(DMRsCGpool)$pValue) %>%
				mutate(Delta=firstsample-secondsample)

			write.table(CGpool,paste0("results/mC/DMRs/",sample1,"__vs__",sample2,"__CG_DMRs_",meth,"_",bs,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

			summary_file<-mutate(CGpool, Type=ifelse(Delta>0, "hyper", "hypo"), Method=meth, Binsize=bs) %>%
					group_by(Type, Method, Binsize) %>%
					summarize(CG_DMRs=n(), .groups = "drop")
		} else {
			summary_file<-tibble::tibble(Type=c("hyper", "hypo"), Method=meth, Binsize=bs, CG_DMRs=c(0, 0))
		}
		
		if (context == "all") {
		
			DMRsCHGpool<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context="CHG", method=meth, binSize=bs, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.2, minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads)
			if ( length(DMRsCHGpool) > 0 ) {
				CHGpool<-data.frame(Chr=seqnames(DMRsCHGpool),Start=start(DMRsCHGpool)-1,End=end(DMRsCHGpool),firstsample=mcols(DMRsCHGpool)$proportion1,secondsample=mcols(DMRsCHGpool)$proportion2, Pvalue=mcols(DMRsCHGpool)$pValue) %>%
				mutate(Delta=firstsample-secondsample)
				write.table(CHGpool,paste0("results/mC/DMRs/",sample1,"__vs__",sample2,"__CHG_DMRs_",meth,"_",bs,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
		
				summary_fileCHG<-mutate(CHGpool, Type=ifelse(Delta>0, "hyper", "hypo"), Method=meth, Binsize=bs) %>%
					group_by(Type, Method, Binsize) %>%
					summarize(CHG_DMRs=n(), .groups = "drop")
			} else {
				summary_fileCHG<-tibble::tibble(Type=c("hyper", "hypo"), Method=meth, Binsize=bs, CHG_DMRs=c(0, 0))
			}
			
			DMRsCHHpool<-computeDMRs(methylationDatasample1pool, methylationDatasample2pool, regions=chrs, context="CHH", method=meth, binSize=bs, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.1, minGap=200, minSize=50, minReadsPerCytosine=3, cores=threads)
			if ( length(DMRsCHHpool) > 0 ) {
				CHHpool<-data.frame(Chr=seqnames(DMRsCHHpool),Start=start(DMRsCHHpool)-1,End=end(DMRsCHHpool),firstsample=mcols(DMRsCHHpool)$proportion1,secondsample=mcols(DMRsCHHpool)$proportion2, Pvalue=mcols(DMRsCHHpool)$pValue) %>%
				mutate(Delta=firstsample-secondsample)
				write.table(CHHpool,paste0("results/mC/DMRs/",sample1,"__vs__",sample2,"__CHH_DMRs_",meth,"_",bs,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	
				summary_fileCHH<-mutate(CHHpool, Type=ifelse(Delta>0, "hyper", "hypo"), Method=meth, Binsize=bs) %>%
					group_by(Type, Method, Binsize) %>%
					summarize(CHH_DMRs=n(), .groups = "drop")
			} else {
				summary_fileCHH<-tibble::tibble(Type=c("hyper", "hypo"), Method=c(meth,meth), Binsize=c(bs,bs), CHH_DMRs=c(0, 0))
			}
					
		summary_file<-merge(summary_file, summary_fileCHG, by=c("Type", "Method", "Binsize"))
		summary_file<-merge(summary_file, summary_fileCHH, by=c("Type", "Method", "Binsize"))
		}
		tot_file <- bind_rows(tot_file, summary_file)
	}
}

tot_file<-mutate(tot_file, Sample=paste0(sample1,"_vs_",sample2)) %>% select(Sample, everything())
write.table(tot_file,paste0("results/mC/DMRs/summary__",sample1,"__vs__",sample2,"__DMRs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
