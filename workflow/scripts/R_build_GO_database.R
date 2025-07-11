#!/usr/bin/env Rscript

library(AnnotationForge)
library(rrvgo)
library(topGO)
library(dplyr)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

info<-read.delim(args[1], header=FALSE)
genes<-read.delim(args[2], header=TRUE) %>%
 rowwise() %>%
 mutate(desc=ifelse(Description=="protein_coding",Type,Description),
        typ=ifelse(Description=="protein_coding",Description,Type)) %>%
 select(-Description, -Type) %>%
 rename(Description=desc, Type=typ)

ref_genome<-args[3]
genus<-args[4]
species<-args[5]
ncbiID<-args[6]
dbname<-paste0("org.",substr(genus,1,1),species,".eg.db")

fGOzm<-unique(info[,c(2,5,7)])
colnames(fGOzm)<-c("GID","GO","EVIDENCE")

fSymzm<-select(genes, GID, Type, Description)
fSymzm$ENTREZID <- paste0("ent",fSymzm$GID)

fChrzm<-select(genes, GID, Chr)

makeOrgPackage(gene_info=fSymzm, chromosome=fChrzm, go=fGOzm,
              version="0.1",
              maintainer="@epicbutton",
              author="@epicbutton",
              outputDir = paste0("./genomes/",refgenome,"/GO/"),
              tax_id = ncbiID,
              genus = genus,
              species = species,
              goTable="go")

db<-paste0("./genomes/",refgenome,"/GO/")
setwd(db)
install.packages(dbname, repos=NULL, type="source")
setwd("../../..")