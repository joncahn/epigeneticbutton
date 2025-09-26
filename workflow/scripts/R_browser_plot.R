#!/usr/bin/env Rscript

library(Gviz)
library(GenomicFeatures)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)

filenames<-read.delim(args[1], header=TRUE)
if ( file.exists(args[2]) ) {

	gff <- import(args[2], format = "gff3")
	gene_types <- c("gene")
	tx_types   <- c("mRNA", "transcript")
	genes_gff  <- gff[gff$type %in% gene_types]
	trans_gff  <- gff[gff$type %in% tx_types]
	exons_gff  <- gff[gff$type == "exon"]

	stopifnot(!is.null(mcols(trans_gff)$ID) || length(trans_gff) == 0)
	stopifnot(!is.null(mcols(trans_gff)$Parent) || length(trans_gff) == 0)
	stopifnot(!is.null(mcols(exons_gff)$Parent))

	parents_list <- strsplit(as.character(mcols(exons_gff)$Parent), ",")
	rep_idx      <- rep(seq_along(exons_gff), lengths(parents_list))
	trans_ids    <- unlist(parents_list, use.names = FALSE)
	exons_expanded <- exons_gff[rep_idx]
	mcols(exons_expanded)$transcript <- trans_ids

	tx_to_gene <- character(0)
	if (length(trans_gff) > 0 && !is.null(mcols(trans_gff)$ID) && !is.null(mcols(trans_gff)$Parent)) {
		tx_ids    <- as.character(mcols(trans_gff)$ID)
		tx_pars   <- as.character(mcols(trans_gff)$Parent)  # usually gene ID
		tx_to_gene <- setNames(tx_pars, tx_ids)
	} else if (length(genes_gff) > 0 && !is.null(mcols(genes_gff)$ID)) {
		overlaps <- findOverlaps(exons_expanded, genes_gff)
		gene_ids <- as.character(mcols(genes_gff)$ID)
		mapped <- rep(NA_character_, length(exons_expanded))
		mapped[queryHits(overlaps)] <- gene_ids[subjectHits(overlaps)]
		mcols(exons_expanded)$symbol <- mapped
	}

	if (length(tx_to_gene) > 0) {
		gene_labels <- tx_to_gene[ mcols(exons_expanded)$transcript ]
		na_idx <- which(is.na(gene_labels) | gene_labels == "")
		if (length(na_idx)) {
			gene_labels[na_idx] <- as.character(mcols(exons_expanded)$transcript[na_idx])
		}
		mcols(exons_expanded)$symbol <- as.character(gene_labels)
	}

	if (any(is.na(mcols(exons_expanded)$symbol))) {
		warning("Some symbol labels are NA; they will be replaced by transcript IDs.")
		na_idx <- which(is.na(mcols(exons_expanded)$symbol))
		mcols(exons_expanded)$symbol[na_idx] <- as.character(mcols(exons_expanded)$transcript[na_idx])
	}
	
	genetrack <- GeneRegionTrack(exons_expanded, name="Genes", shape="smallArrow", col="black", fill="grey60", rotation.title=0, cex.title=0.5, lwd=0.1, collapseTranscripts="meta", showId=TRUE, stacking="dense",transcriptAnnotation="symbol")

#	genes<-txdbmaker::makeTxDbFromGFF(args[2], format="gff")
# 	genetrack <- GeneRegionTrack(genes, name="Genes", shape="smallArrow", col="black", fill="grey60", rotation.title=0, cex.title=0.5, lwd=0.1, collapseTranscripts="meta", showId=TRUE, stacking="dense", transcriptAnnotation="symbol")
} else {
	genes<-c()
}

tes<-import(args[3], format="bed")
tetrack<-AnnotationTrack(tes, name="TEs", stacking = "dense", fill = "lightgreen", shape="box", rotation.title=0, cex.title=0.5, lwd=0.1)

plotname<-args[4]
pdfmame<-args[5]

tot<-nrow(filenames)

htcol<-c()
htcol2<-c()
if ( length(args) == 7 ) {
	htstarttable<-read.delim(args[6], header=FALSE)
	htwidthtable<-read.delim(args[7], header=FALSE)
	colors<-c("#B7E2FD","#fac0c7","#fac0c7","#fac0c7","#fac0c7","#fac0c7")
	colors2<-c("#F6FBFE","#fffafa","#fffafa","#fffafa","#fffafa","#fffafa")
	htstart<-c()
	htwidth<-c()
	for ( i in c(1:nrow(htstarttable)) ) {
		htstart<-c(htstart,htstarttable[i,])
		htwidth<-c(htwidth,htwidthtable[i,])
		htcol<-c(htcol,colors[i])
		htcol2<-c(htcol2,colors2[i])
	}
}

options(ucscChromosomeNames=FALSE)

tracksize<-c(1,1,0.5)
tracklist<-list()
for ( i in c(1:tot) ) {
	label<-filenames$Name[i]
	path<-filenames$Path[i]
	backcolor<-filenames$Backcolor[i]
	trackcolor<-filenames$Trackcolor[i]
	fillcolorplus<-filenames$Fillcolorplus[i]
	fillcolorminus<-filenames$Fillcolorminus[i]
	ymin<-filenames$Ymin[i]
	ymax<-filenames$Ymax[i]
	ymintick<-sign(ymin)*((floor(abs(ymin)*100)/100))
	ymaxtick<-sign(ymax)*((floor(abs(ymax)*100)/100))
	tracksize<-c(tracksize,1)
	print(paste0("Importing bw for ",label))
	bw<-import(path)
	print(paste0("Creating track for ",label))
	track<-DataTrack(bw, type="polygon", baseline=0, name=label, background.title = backcolor, col=trackcolor, fill.mountain=c(fillcolorminus,fillcolorplus), col.baseline="grey50", ylim=c(ymin,ymax), yTicksAt=c(ymintick,ymaxtick), rotation.title=0, cex.title=0.5, lwd=0.01, hjust.title=1)
	tracklist<-append(tracklist, track)
}

axistrack<-GenomeAxisTrack(scale=0.1, labelPos="above")

if ( length(htcol) > 0 ) {
	httrack <- HighlightTrack(trackList=tracklist, start=htstart, width=htwidth, col=htcol, fill=htcol2)
	pdf(pdfmame, width = 12, height = tot)
	plotTracks(list(axistrack, genetrack, tetrack, httrack), sizes=tracksize, main=plotname)
	dev.off()
} else {
	tracks<-append(list(axistrack, genetrack, tetrack), tracklist)
	pdf(pdfmame,paper="a4", width = 12, height = tot)
	plotTracks(tracks, sizes=tracksize, main=plotname)
	dev.off()
}





