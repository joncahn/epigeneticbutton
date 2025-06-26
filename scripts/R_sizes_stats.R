#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

summary_stats<-args[1]
analysisname<-args[2]
minsize<-args[3]
maxsize<-args[4]

plot.sRNA.sizes<-function(stattable, sizemin, sizemax) {
	
	count<-read.delim(stattable, header = TRUE) %>%
		filter(Count>=sizemin & Count<=sizemax) 
	
	rdvalue<-count$Size[1]
	rdsample<-count$Sample[1]
	if (! "filtered" %in% Count$Type) {
		count<-rbind(rdsample,"filtered",rdvalue,0)
	}
	count<-spread(count, key = Type, value=Count) %>%
		mutate(trim=trimmed-filtered, filt=filtered-mapped) %>%
		select(-trimmed, -filtered) %>%
		rename(trimmed=trim, filtered=filt) %>%
		gather(Type, Count, filtered, trimmed, mapped)
	count$Count<-replace_na(count$Count, 0)
	count$Type<-factor(count$Type, levels=c("trimmed","filtered","mapped"))

	a<-seq(sizemin, sizemax, by = 10)
	breaksarray<-sort(unique(c(a, 21, 24)))
	
	plot <- ggplot(count, aes(Size, Count, fill=Type)) +
			geom_bar(stat="identity", position="stack", color="black", size=0.01) +
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
				panel.grid.major.x = element_line(colour="lightgrey",size=0.1),
				panel.background = element_rect(fill = "white", colour = "black"),
				strip.background = element_rect(fill = 'white', colour = 'black'),
				legend.key=element_blank())
	plot
}  

pdf(paste0("combined/plots/srna_sizes_stats_",analysisname,"_sRNA.pdf"), height=10, width=12)
plot.sRNA.sizes(summary_stats, minsize, maxsize)
dev.off()