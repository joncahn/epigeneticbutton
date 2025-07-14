#!/usr/bin/env Rscript

library(edgeR)
library(AnnotationForge)
library(rrvgo)
library(dplyr)
library(topGO)
library(purrr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

dbname<-args[1]
analysisname<-args[2]
refgenome<-args[3]
targetfile<-args[4]
backgroundfile<-args[5]

db<-paste0("./genomes/",refgenome,"/GO/")
setwd(db)
library(dbname, character.only = TRUE)
setwd("../../..")

getGO<-function(genelist, target, ont, name) {
	GOdata<-new("topGOdata", 
				ontology = ont, 
				allGenes = genelist,
				annot = annFUN.gene2GO, 
				gene2GO = gene2GO)
	resultFisher<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
	resultFisherSummary <- summary(attributes(resultFisher)$score <= 0.01)
	nSigTerms<-0
	if (length(resultFisherSummary) == 3) {
		nSigTerms<-as.integer(resultFisherSummary[[3]])
	}
	summary<-GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nSigTerms, numChar=1000)
	tab<-summary %>%
		rename_with(.cols = starts_with("apply"), .fn = ~ { if (length(.) > 0) { paste0("classicFisher") } else { . } }) %>%
		mutate(classicFisher = classicFisher %>% str_replace(pattern= "< *1e-30", replacement = "1e-30") %>% as.numeric())
	sigTerms<-tab$GO.ID
	genesInTerms<-genesInTerm(GOdata, sigTerms)
	genesInTerms2<-map(genesInTerms, ~ intersect(.x, myInterestedGenes) %>% paste(collapse = ","))
	tab2<-tab %>% 
		left_join(tibble(GO.ID = names(genesInTerms2), genes = genesInTerms2) %>% 
		tidyr::unnest(genes), by = "GO.ID")
	tab3<-tab %>%
		rename(GO=GO.ID) %>%
		merge(target, by="GID") %>%
		arrange(GO) %>%
		unique()
	if (nrow(tab2) > 1) {
		write.table(tab2,paste0("results/RNA/GO/topGO_",name,"_",ont,"_GOs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	if (nrow(tab3) > 0) {
		write.table(tab3,paste0("results/RNA/GO/topGO_",name,"_",ont,"_GIDs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}	
  
	scores<-setNames(-log10(as.numeric(tab$classicFisher)), tab$GO.ID)
	reducedTerms <- tab2 %>%
			rename("go" = "GO.ID", "term" = "Term") %>%
			mutate(parentTerm = term, score = scores)
  
	if (nrow(tab) > 1 ) {
		simMatrix<-calculateSimMatrix(tab$GO.ID,
									orgdb=dbname,
									ont=ont,
									method="Rel")
		if (! is.null(nrow(simMatrix))) {
			reducedTerms<-reduceSimMatrix(simMatrix,
										scores,
										threshold = 0.7,
										orgdb=dbname)
		}
	}
	pdf(paste0("results/RNA/plots/topGO_",name,"_",ont,"_treemap.pdf"), width=8, height=8)
	treemapPlot(reducedTerms, size = "score")
	dev.off()
}

if (startsWith(backgroundfile, "results/RNA/DEG/counts__")) {
	genecount<-read.delim(backgroundfile, header = TRUE, row.names = "GID")
	target<-read.delim(targetfile, header = TRUE)

	keep.exprs<-rowSums(cpm(genecount)>1)>=2
	filtered<-genecount[keep.exprs,]
	filtered$GID<-row.names(filtered)

	allGenes<-unique(unlist(filtered$GID))
	
	for ( samp in unique(sampletable$Sample) ) {
		for deg in c("UP","DOWN") {
			samplename<-paste0(deg,"_in_",samp)
			for ( ont in c("BP","MF") ) {
				print(paste0("Getting ",ont," for ",samplename))
				sampletable<-filter(target, Sample==samp, DEG==deg)
				myInterestedGenes<-unique(unlist(sampletable$GID))
				geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
				names(geneList)<-allGenes
				getGO(geneList, target, ont, samplename)
			}
		}
	}

} else if (startsWith(backgroundfile, "results/combined/tracks/")) {
	
} else {

}




