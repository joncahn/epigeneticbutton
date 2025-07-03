#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

genecount<-read.delim(args[1], header = TRUE)
genecount<-filter(genecount, ! grepl("^N_", GID))

targets<-read.delim(args[2], header = TRUE)
samples<-as.factor(targets$Sample)
reps<-as.factor(targets$Replicate)
genotypes<-unique(samples)

analysisname<-args[3]

refgenome<-args[4]

ref_genes<-read.delim(args[5], header = FALSE, 
                      col.names = c("Chr","Start","Stop","Name","Value","Strand"))
ref_genes<-mutate(ref_genes, GID=str_replace(ref_genes$Name, pattern = ".*ID=(gene:)?([^;]+).*", replacement = "\\2")) %>%
  select(-Name, -Value)
ref_genes$GID<-str_remove_all(ref_genes$GID, pattern = "_.")

all_rpkm<-data.frame()
for (sample1 in genotypes) {
	nbreps<-sum(grepl(sample1, reps))
	tmp<-genecount[, grepl(sample1, colnames(genecount)) | colnames(genecount) == "GID"]
	tmp<-mutate(tmp, avg = rowMeans(select(tmp, -GID)))
	tmp2<-merge(tmp, ref_genes) %>%
		mutate(RPKM=(avg*1000)/(Stop-Start), Sample=sample1) %>%
		select(GID,Sample,RPKM)
	all_rpkm<-rbind(all_rpkm,tmp2)
}

write.table(all_rpkm,paste0("results/RNA/DEG/genes_rpkm__",analysisname,"__",refgenome,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
