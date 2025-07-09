#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
			
args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]
refgenome<-args[2]
targetfile<-args[3]
name<-str(args[4])

load(paste0("results/RNA/DEG/ReadyToPlot__",analysisname,"__",refgenome,".RData"))

genelist<-read.delim(targetfile, header = FALSE)

if (name == "unique_DEGs") {
	for 
	pdf(paste0("results/combined/plots/plot_expression_",analysisname,"_",refgenome,"_",name,".pdf"), height=8, width=8)
	for (i in 1:(min(nrow(genelist),100)) {
		gene<-genelist[i,1]
		if ( gene != "GID" ) {
			label <- paste0("unique ",genelist[i,3]," in ",genelist[i,2])
			print(plot.Expression(gene, label))
		}
	}
	dev.off()
} else {
	pdf(paste0("results/combined/plots/plot_expression_",analysisname,"_",refgenome,"_",name,".pdf"), height=8, width=8)
	for (i in 1:(nrow(genelist)) {
		gene<-genelist[i,1]
		if ( gene != "GID" ) {
			label <- if (ncol(genelist) >= 2) genelist[i,2] else "NoLabel"
			print(plot.Expression(gene, label))
		}
	}
	dev.off()
}