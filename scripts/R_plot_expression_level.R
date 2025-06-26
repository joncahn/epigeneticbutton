#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
			
args = commandArgs(trailingOnly=TRUE)

analysisname<-args[1]
refgenome<-args[2]
targetfile<-args[3]
filename<-args[4]

if (filename == "") {
	name<-tools::file_path_sans_ext(basename(targetfile))
} else {
	name<-c(filename)
}

load(paste0("RNA/DEG/ReadyToPlot__",analysisname,"__",refgenome,".RData"))

genelist<-read.delim(targetfile, header = FALSE)

pdf(paste0("combined/plots/plot_expression_",analysisname,"_",refgenome,"_",name,".pdf"), height=8, width=8)
for (i in 1:(nrow(genelist))) {
	gene<-genelist[i,1]
	label<-genelist[i,2]
	print(plot.Expression(gene, label))
}
dev.off()