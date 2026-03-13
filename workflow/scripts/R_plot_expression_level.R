#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
			
args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]
refgenome<-args[2]
targetfile<-args[3]
name<-args[4]
output_dir<-args[5]

load(paste0(output_dir,"/RNA/DEG/ReadyToPlot__",analysisname,"__",refgenome,".RData"))

genelist<-read.delim(targetfile, header = FALSE)

if (name == "unique_DEGs") {
	pdf(paste0(output_dir,"/RNA/plots/plot_expression__",analysisname,"__",refgenome,"__",name,".pdf"), height=8, width=8)
	for (i in 1:(min(nrow(genelist),100))) {
		gene<-genelist[i,1]
		if ( gene != "GID" ) {
			label <- "NoLabel"
			print(plot.Expression(gene, label))
		}
	}
	dev.off()
} else {
	pdf(paste0(output_dir,"/RNA/plots/plot_expression__",analysisname,"__",refgenome,"__",name,".pdf"), height=8, width=8)
	for (i in 1:(nrow(genelist))) {
		gene<-genelist[i,1]
		if ( gene != "GID" ) {
			label <- if (ncol(genelist) >= 2) genelist[i,2] else "NoLabel"
			print(plot.Expression(gene, label))
		}
	}
	dev.off()
}