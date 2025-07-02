#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

statfile<-args[1]
analysisname<-args[2]

summary_stats<-read.delim(statfile, header = TRUE)
summary_stats$Count<-as.numeric(summary_stats$Count)
minsize<-min(summary_stats$Size)
maxsize<-max(summary_stats$Size)

plot.sRNA.sizes<-function(stattable, sizemin, sizemax) {
	
	count<-filter(stattable, Size>=minsize & Size<=maxsize)
	count$Count<-as.numeric(count$Count)
	count<-pivot_wider(count, names_from = Type, values_from = Count) 

	if (! "filtered" %in% colnames(count)) {
		count<-mutate(count, filtered=mapped)
	}

	count<-mutate(count, trim=trimmed-filtered, filt=filtered-mapped) %>%
		select(-trimmed, -filtered) %>%
		rename(trimmed=trim, filtered=filt) %>%
		pivot_longer(cols = c(filtered, trimmed, mapped), names_to = "Type", values_to = "Count")
	count$Type<-factor(count$Type, levels=c("trimmed","filtered","mapped"))

	a<-seq(sizemin, sizemax, by = 10)
	breaksarray<-sort(unique(c(a, 21, 24)))

	plot <- ggplot(count, aes(Size, Count, fill=Type)) +
				geom_bar(stat="identity", position="stack", color="black", linewidth=0.01) +
				facet_wrap(~Sample, nrow = length(unique(count$Sample)), scales="free_y") +
				scale_fill_manual(labels=c("trimmed"="post-trimming","filtered"="post-filtering","mapped"="mapped"), 
						values = c("trimmed"="grey","filtered"="blue","mapped"="green")) +
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

pdf(paste0("combined/plots/srna_sizes_stats_",analysisname,"_sRNA.pdf"), height=10, width=12)
plot.sRNA.sizes(summary_stats, minsize, maxsize)
dev.off()

pdf(paste0("combined/plots/srna_sizes_stats_zoom_",analysisname,"_sRNA.pdf"), height=10, width=12)
plot.sRNA.sizes(summary_stats, 20, 25)
dev.off()
