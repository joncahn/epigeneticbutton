#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

statfile<-args[1]
analysisname<-args[2]
zoommin<-args[3]
zoommax<-args[4]

summary_stats<-read.delim(statfile, header = TRUE)
summary_stats$Count<-as.numeric(summary_stats$Count)
minsize<-min(summary_stats$Size)
maxsize<-max(summary_stats$Size)

tot<-length(unique(summary_stats$Sample))

plot.sRNA.sizes<-function(stattable, sizemin, sizemax) {
	
	count<-filter(stattable, Size>=sizemin & Size<=sizemax)
	count$Count<-as.numeric(count$Count)
	count<-pivot_wider(count, names_from = Type, values_from = Count) 

	if (! "deduplicated" %in% colnames(count)) {
		count<-mutate(count, deduplicated=trimmed)
	}
	
	if (! "filtered" %in% colnames(count)) {
		count<-mutate(count, filtered=deduplicated)
	}

	count<-mutate(count, dedup=trimmed-deduplicated, struc=deduplicated-filtered, unmap=filtered-mapped) %>%
		select(-trimmed, -deduplicated, -filtered) %>%
		rename(deduplicated=dedup, structural_variants=struc, unmapped=unmap) %>%
		pivot_longer(cols = c(deduplicated, structural_variants, unmapped, mapped), names_to = "Type", values_to = "Count")
	count$Type<-factor(count$Type, levels=c("deduplicated","structural_variants","unmapped","mapped"))
	count$Sample <- factor(count$Sample, levels = sort(unique(count$Sample)))

	if ( sizemin == 20 && sizemax == 25) {
		breaksarray<-c(20, 21, 22, 23, 24, 25)
	} else {	
		a<-seq(sizemin, sizemax, by = 10)
		breaksarray<-sort(unique(c(a, 21, 24)))
	}

	plot <- ggplot(count, aes(Size, Count, fill=Type)) +
				geom_bar(stat="identity", position="stack", color="black", linewidth=0.01) +
				facet_wrap(~Sample, ncol=1, scales="free_y") +
				scale_fill_manual(values = c("deduplicated"="grey","structural_variants"="purple","unmapped"="blue","mapped"="darkgreen")) +
				labs(y="Counts", x="Sizes", fill="") +
				scale_x_continuous(breaks = breaksarray) +
				theme(axis.title.y=element_text(size=15), 
					axis.title.x=element_text(size=15),
					axis.text.x=element_text(size=10),
					panel.grid.major.y = element_line(colour="lightgrey"), 
					panel.grid.minor.y = element_blank(),
					panel.grid.major.x = element_line(colour="lightgrey",linewidth=0.1),
					panel.background = element_rect(fill = "white", colour = "black"),
					strip.background = element_rect(fill = 'white', colour = 'black'),
					legend.key=element_blank())
	plot
}  

pdf(paste0("results/combined/plots/srna_sizes_stats_",analysisname,"_sRNA.pdf"), height=tot*2, width=12)
plot.sRNA.sizes(summary_stats, minsize, maxsize)
dev.off()

pdf(paste0("results/combined/plots/srna_sizes_stats_zoom_",analysisname,"_sRNA.pdf"), height=tot*2, width=12)
plot.sRNA.sizes(summary_stats, zoommin, zoommax)
dev.off()
