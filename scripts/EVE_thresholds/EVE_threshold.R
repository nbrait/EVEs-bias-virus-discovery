setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/THRESHOLD")
#install.packages("devtools")

library(seqinr)
library (devtools)
library (tidyverse)

##################################
#EVE_dataframe
EVE_matches <- read.table("EVE_matches_edited.csv", header = FALSE)
EVE_matches = subset(EVE_matches, select = -c(V4,V5,V6,V7,V8,V9))
colnames(EVE_matches) <- c("SAMPLE", "EVE","EVE_SIMILARITY", "EVE_BITSCORE", "FRAMESHIFT_EVE")
EVE_matches$EVE_SIMILARITY <- sub("^", "-", EVE_matches$EVE_SIMILARITY)
EVE_matches <- transform(EVE_matches, EVE_SIMILARITY = as.numeric(EVE_SIMILARITY))
EVE_matches$EVE_BITSCORE <- sub("^", "-", EVE_matches$EVE_BITSCORE)
EVE_matches <- transform(EVE_matches, EVE_BITSCORE = as.numeric(EVE_BITSCORE))
#EXO dataframe
EXO_matches <- read.table("EXO_matches_edited.csv", header = FALSE)
EXO_matches = subset(EXO_matches, select = -c(V4,V5,V6,V7,V8,V9))
colnames(EXO_matches) <- c("SAMPLE", "EXO","EXO_SIMILARITY", "EXO_BITSCORE", "FRAMESHIFT_EXO")
# merging dataframes
matches_table <- merge(EXO_matches, EVE_matches, by = 'SAMPLE' ) 
matches_table$DIFF_SIMILARITY <- matches_table$EXO_SIMILARITY+(matches_table$EVE_SIMILARITY)
matches_table$DIFF_BITSCORE <- matches_table$EXO_BITSCORE+(matches_table$EVE_BITSCORE)
write.table(matches_table, file = "EVE_THRESHOLDS.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
matches_table <- read.table("EVE_THRESHOLDS_NEW", header = TRUE, sep = "\t")
ggplot(matches_table) + 
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color=EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY), shape=FRAMESHIFT_EXO), size=4) +
  scale_color_viridis_c(name = "identity to known EVE") +
  scale_y_log10() +
  scale_shape_manual(values=c(20, 24), name = "frameshifts present")+
  theme(legend.position = c(0.7, 0.2),
        legend.direction = "vertical",
        legend.key=element_blank(),
        legend.title.align = 1,
        #legend.direction = legend.direction,
        #legend.box.just = "center",
        legend.background=element_blank(),
        panel.background = element_rect(fill = "white"), 
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +
          xlab("Δ BLAST bit score") +# for the x axis label
          ylab("maximum BLAST bit score per contig (log10)") 
#  guide_legend(title.position = "top")

