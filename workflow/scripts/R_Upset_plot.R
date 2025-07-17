#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)
library(ComplexUpset)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

merged<-read.delim(args[1], header = TRUE)
annotated<-read.delim(args[2], header = TRUE) %>%
	mutate(distance=abs(Distance)+1) %>%
	select("PeakID","Category",Gap=distance) %>%
	rename(Distance=Gap)
env<-args[3]
types<-args[4]
output<-args[5]

sampleslist<-unique(unlist(strsplit(merged$Samples, ",")))

mat<-separate_rows(merged, Samples, sep = ",") %>%
	mutate(value=1) %>%
	pivot_wider(names_from = Samples, values_from = value, values_fill = 0) %>%
	merge(annotated, by="PeakID")
	
mat$Category<-factor(mat$Category, levels=c("Distal_downstream","Terminator","Gene_body","Promoter","Distal_upstream"))

#
colorlist<-unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
i<-1
queries<-c()
listcolor<-c()
for (sampletype in types) {
	setcols<-colnames(mat)[grep(sampletype, colnames(mat))]
	combos<-map(seq_along(setcols), ~ combn(setcols, ., FUN = c, simplify = FALSE)) %>% 
			unlist(recursive = FALSE)
	
	tmpqueries<-map(combos, ~ upset_query(intersect = .x, color = colorlist[i], 
                               fill = colorlist[i], only_components = c('intersections_matrix')))
	queries<-append(queries, tmpqueries)
	listcolor<-append(listcolor, colorlist[i])
	i<-i+1
}
colmarks<-setNames(listcolor, types)
colmarks["Mix"] <- "black"

type_cols <- lapply(types, function(t) { grep(t, colnames(mat), value = TRUE) })
names(type_cols) <- types

inputable <- inputable %>% 
  mutate(
    exclusive_mark = case_when(
      !!!setNames(lapply(types, function(t) {
			rowSums(across(all_of(type_cols[[t]]))) == rowSums(across(all_of(sampleslist)))
		}), types),
		TRUE ~ as.character("Mix")
    )
  ) %>% 
  relocate(exclusive_mark, .after = Category)

plot<-upset(mat, sampleslist, name="Peaks", 
      mode='exclusive_intersection',
      n_intersections=30, 
      sort_sets=FALSE,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared peaks'=intersection_size(
          counts=FALSE, mapping=aes(fill=Category)) +
          scale_fill_manual(values=c("Distal_downstream"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","Distal_upstream"="#2e2e2e"),
                            name="Distance category")
      ),
      annotations = list(
        'Distance to closest gene' = (
          ggplot(mapping = aes(x=intersection, y=Distance, fill = exclusive_mark)) +
            geom_violin(scale="width", na.rm=TRUE, color = "black") +
            scale_y_continuous(trans = "log10",
                               labels=scales::label_number_si(accuracy = 1, unit = "bp")) +
            scale_fill_manual(values=colmarks, name="Exclusive marks"))
      ),
      queries = queries,
      set_sizes = (upset_set_size() + ylab("Total peaks") +
        theme(axis.text.x = element_text(angle = 45))),
      matrix = (intersection_matrix(geom = geom_point(shape = "circle",size = 3),
          segment = geom_segment(size = 1.5),
          outline_color = list(active = alpha("white", 0),inactive = alpha("white", 0))) +
          scale_color_manual(values = c("TRUE" = "black", "FALSE" = alpha("white", 0)),
            labels = c("TRUE" = "yes", "FALSE" = "no"),
            breaks = c("TRUE", "FALSE"),
            guide = FALSE) +
          theme(axis.ticks = element_blank(),
            panel.grid = element_blank())
      ),
      themes = upset_modify_themes(
        list(
          "default" = theme(
            panel.grid.major.x = element_blank(),
            axis.ticks.y = element_line(size = 0.25, color = "#2e2e2e")
          ),
          "intersections_matrix" = theme(
            panel.grid = element_blank(),
            panel.grid.major.y = element_line(color = c("#CFCCCF", "white"), size = 5)
          ),
          "Intersection size" = theme(
            panel.grid = element_blank(),
          ),
          "overall_sizes" = theme(
            panel.grid = element_blank(),
            axis.ticks.x = element_line(size = 0.25, color = "#2e2e2e")
          )
        )
      ),
      stripes = alpha("white", 0)
)

pdf(output,10,8)
print(plot)
dev.off()
