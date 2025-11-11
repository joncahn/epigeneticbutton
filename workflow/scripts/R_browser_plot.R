#!/usr/bin/env Rscript

library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(txdbmaker)

args = commandArgs(trailingOnly=TRUE)

filenames<-read.delim(args[1], header=TRUE)
if ( file.info(args[2])$size > 0 ) {
	genes<-txdbmaker::makeTxDbFromGFF(args[2], format="gff")
} else {
	genes<-c()
}

tes<-import(args[3], format="bed")

plotname<-args[4]
pdfname<-args[5]

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
	lab<-gsub("_", "\n", label)
	track<-DataTrack(bw, type="polygon", baseline=0, name=lab, background.title = backcolor, col=trackcolor, fill.mountain=c(fillcolorminus,fillcolorplus), col.baseline="grey50", ylim=c(ymin,ymax), yTicksAt=c(ymintick,ymaxtick), cex.title=0.5, lwd=0.01, fontsize=12, fontcolor.title="black", col.axis="black", cex.axis=0.5)
	tracklist<-append(tracklist, track)
}

axistrack<-GenomeAxisTrack(scale=0.1, labelPos="above")
genetrack<-GeneRegionTrack(genes, name="Genes", shape="smallArrow", col="black", fill="grey60", rotation.title=0, cex.title=0.5, lwd=0.1, fontsize=12, fontcolor.title="black", collapseTranscripts=FALSE, showId=TRUE, stacking="squish")
tetrack<-AnnotationTrack(tes, name="TEs", stacking = "dense", fill = "lightgreen", shape="box", rotation.title=0, cex.title=0.5, lwd=0.1, fontsize=12, fontcolor.title="black")

if ( length(htcol) > 0 ) {
	httrack <- HighlightTrack(trackList=tracklist, start=htstart, width=htwidth, col=htcol, fill=htcol2)
	pdf(pdfname, width = 12, height = tot)
	plotTracks(list(axistrack, genetrack, tetrack, httrack), sizes=tracksize, main=plotname, cex.main = 1)
	dev.off()
} else {
	tracks<-append(list(axistrack, genetrack, tetrack), tracklist)
	pdf(pdfname, paper="a4", width = 12, height = tot)
	plotTracks(tracks, sizes=tracksize, main=plotname, cex.main = 1)
	dev.off()
}





