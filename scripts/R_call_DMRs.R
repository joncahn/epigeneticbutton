#!/usr/bin/env Rscript

library(DMRcaller)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

chromsizes<-read.table(args[1], col.names = c("chr", "length"))
sample1<-args[2]
sample2<-args[3]

chrs<-GRanges(seqnames = chromsizes$chr, ranges = IRanges(start = 1, end = chromsizes$length))

# methylationDataddm1hirapool<-readBismarkPool(c("CX_reports/ddm1hira_A_Rep1.deduplicated.CX_report.txt.gz","CX_reports/ddm1hira_A_Rep2.deduplicated.CX_report.txt.gz"))
# methylationDataddm1pool<-readBismarkPool(c("CX_reports/ddm1_A_Rep1.deduplicated.CX_report.txt.gz","CX_reports/ddm1_A_Rep2.deduplicated.CX_report.txt.gz"))
# DMRsCGpool<-computeDMRs(methylationDataddm1hirapool, methylationDataddm1pool, regions=tair10_chrs, context="CG", method="bins", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3, cores=6)

# CGpool<-data.frame(Chr=seqnames(DMRsCGpool),Start=start(DMRsCGpool)-1,End=end(DMRsCGpool),ddm1hira=elementMetadata(DMRsCGpool)[,3],ddm1=elementMetadata(DMRsCGpool)[,6], Pvalue=elementMetadata(DMRsCGpool)[,10]) %>%
	# mutate(Delta=ddm1hira-ddm1)

# methylationDataddm1hirarep1<-readBismarkPool("CX_reports/ddm1hira_A_Rep1.deduplicated.CX_report.txt.gz")
# methylationDataddm1rep1<-readBismark("CX_reports/ddm1_A_Rep1.deduplicated.CX_report.txt.gz")
# DMRsCGrep1<-computeDMRs(methylationDataddm1hirarep1, methylationDataddm1rep1, regions=tair10_chrs, context="CG", method="bins", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3, cores=6)

# CGrep1<-data.frame(Chr=seqnames(DMRsCGrep1),Start=start(DMRsCGrep1)-1,End=end(DMRsCGrep1),ddm1hirarep1=elementMetadata(DMRsCGrep1)[,3],ddm1rep1=elementMetadata(DMRsCGrep1)[,6], Pvaluerep1=elementMetadata(DMRsCGrep1)[,10]) %>%
	# mutate(Delta=ddm1hirarep1-ddm1rep1)
	
# methylationDataddm1hirarep2<-readBismarkPool("CX_reports/ddm1hira_A_Rep2.deduplicated.CX_report.txt.gz")
# methylationDataddm1rep2<-readBismark("CX_reports/ddm1_A_Rep2.deduplicated.CX_report.txt.gz")
# DMRsCGrep2<-computeDMRs(methylationDataddm1hirarep2, methylationDataddm1rep2, regions=tair10_chrs, context="CG", method="bins", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=0.3, minGap=200, minSize=50, minReadsPerCytosine=3, cores=6)

# CGrep2<-data.frame(Chr=seqnames(DMRsCGrep2),Start=start(DMRsCGrep2)-1,End=end(DMRsCGrep2),ddm1hirarep2=elementMetadata(DMRsCGrep2)[,3],ddm1rep2=elementMetadata(DMRsCGrep2)[,6], Pvaluerep2=elementMetadata(DMRsCGrep2)[,10]) %>%
	# mutate(Deltarep2=ddm1hirarep2-ddm1rep2)

# write.table(CGrep1,"DMRs_CG_rep1.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(CGrep2,"DMRs_CG_rep2.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(CGpool,"DMRs_CG_pool.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
