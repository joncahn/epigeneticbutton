#!/usr/bin/env Rscript

library(AnnotationForge)
library(rrvgo)
library(topGO)
library(dplyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

gaf<-read.delim(args[1], header=FALSE)
genes<-read.delim(args[2], header=TRUE) %>%
 rowwise() %>%
 mutate(desc=ifelse(Description=="protein_coding",Type,Description),
        typ=ifelse(Description=="protein_coding",Description,Type)) %>%
 select(-Description, -Type) %>%
 rename(Description=desc, Type=typ)

refgenome<-args[3]
genus<-args[4]
species<-args[5]
ncbiID<-args[6]
dbname<-paste0("org.",substr(genus,1,1),species,".eg.db")

fGO<-unique(gaf[,c(1,6,10)])
colnames(fGO)<-c("GID","GO","EVIDENCE")

fSym<-unique(select(genes, GID, Type, Description))
fSym$ENTREZID<-paste0("ent",fSym$GID)

fChr<-unique(select(genes, GID, Chr))

makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
              version="0.1",
              maintainer="user <user@epicbutton>",
              author="user <user@epicbutton>",
              outputDir = paste0("./genomes/",refgenome,"/GO/"),
              tax_id = ncbiID,
              genus = genus,
              species = species,
              goTable="go")

db<-paste0("./genomes/",refgenome,"/GO/")
setwd(db)
install.packages(dbname, repos=NULL, type="source")
setwd("../../..")