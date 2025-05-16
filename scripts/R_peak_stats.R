#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

summary_stats<-args[1]
analysisname<-args[2]
outputfile<-args[3]
env<-args[4]

plot.peak.stats<-function(stattable, name) {
  table<-read.delim(stattable, header = TRUE, sep="\t") %>%
    mutate(Label=paste(Line,Tissue,Sample)) %>%
    separate(Selected_peaks, into=c("Selected","temp1"), sep=" ") %>%
    rename(Rep1=Peaks_in_Rep1,Rep2=Peaks_in_Rep2,Merged=Peaks_in_merged, Pseudoreps=Peaks_in_pseudo_reps) %>%
    select(-temp1)
  table$Line<-as.factor(table$Line)
  table$Tissue<-as.factor(table$Tissue)
  table$Sample<-as.factor(table$Sample)
  table$Selected<-as.numeric(table$Selected)
  table<-gather(table, key="Peak_type",value="Number",Rep1,Rep2,Merged,Pseudoreps,Selected)
  table$Peak_type<-factor(table$Peak_type, levels=c("Rep1","Rep2","Merged","Pseudoreps","Selected"))
  table<-arrange(table, desc(Line), desc(Sample), desc(Tissue))
  
  plot<-ggplot(table, aes(Sample,Number,fill=Peak_type)) +
    geom_bar(stat="identity", position="dodge", color="black", show.legend = T) +
    labs(title=paste("Number of peaks in each Histone ChIPseq sample of",analysisname), 
         x="",y="Number of peaks", fill="Peaks in:") +
    scale_fill_manual(values = brewer.pal(5,"Paired"), guide = "legend") +
    facet_grid(~Line+Tissue) +
    theme(axis.text.x=element_text(color="black",size=10, angle=90, vjust=0.5, hjust = 1),
          title = element_text(size=15),
          axis.title.y = element_text(size=12),
          legend.title = element_text(size=12),
          panel.grid=element_blank(),
          panel.grid.major.y = element_line(color="grey"),
          panel.grid.minor.y = element_line(color="grey"),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          strip.text.x = element_text(size=12),
          strip.background = element_blank())

  plot
}

pdf(outputfile, height=10, width=12)
plot.peak.stats(summary_stats, analysisname)
dev.off()

