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
ncbi_taxid<-args[6]

dbname<-paste0("org.",substr(genus,1,1),species,".eg.db")

fGO<-unique(gaf[,c(1,6,10)])
colnames(fGO)<-c("GID","GO","EVIDENCE")

fSym<-unique(select(genes, GID, Type, Description))
fSym$ENTREZID<-paste0("ent",fSym$GID)

fChr<-unique(select(genes, GID, Chr))

godir <- normalizePath(paste0("./genomes/",refgenome,"/GO"), mustWork=FALSE)
dir.create(godir, showWarnings=FALSE, recursive=TRUE)

# makeOrgPackage writes a source package directory; install.packages with
# repos=NULL,type="source" refuses to install when lib equals the source
# package's parent dir ("cannot install to srcdir"). Avoid the conflict by
# writing the source to a sibling temp dir, then installing from there.
src_tmp <- tempfile(pattern="go_src_", tmpdir=dirname(godir))
dir.create(src_tmp)
on.exit(unlink(src_tmp, recursive=TRUE), add=TRUE)

makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
              version="0.1",
              maintainer="user <user@epicc>",
              author="user <user@epicc>",
              outputDir = src_tmp,
              tax_id = ncbi_taxid,
              genus = genus,
              species = species,
              goTable="go")

# Install into the per-genome GO directory and prepend it to .libPaths()
# so this run loads its own org.<G><species>.eg.db, isolated from any
# same-named package built for a different reference genome with the
# same binomial. AnnotationForge fixes the package name to
# org.<G><species>.eg.db, so the only way to keep multiple genomes
# coexisting is to scope each install to its own library directory.
install.packages(file.path(src_tmp, dbname), repos=NULL, type="source", lib=godir)
.libPaths(c(godir, .libPaths()))
