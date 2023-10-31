setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/THRESHOLD/genomic")
#install.packages("devtools")

library(seqinr)
library (devtools)
library (tidyverse)
library(ggplot2)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")


##################################
##### genomic ######
#EVE_dataframe
#EVE_matches <- read.csv("EVE_genomic_matches.csv", header = FALSE)
EVE_matches <- read.csv("EVE_gen_matches.tsv", header = FALSE, sep = "\t")
EVE_matches = subset(EVE_matches, select = -c(V4,V5,V6,V7,V8,V9,V11))
#colnames(EVE_matches) <- c("SAMPLE", "EVE","EVE_SIMILARITY", "EVE_BITSCORE", "FRAMESHIFT_EVE")
colnames(EVE_matches) <- c("SAMPLE", "EVE","EVE_SIMILARITY", "EVE_BITSCORE", "FRAMESHIFT_EVE")
EVE_matches$EVE_SIMILARITY <- sub("^", "-", EVE_matches$EVE_SIMILARITY)
EVE_matches <- transform(EVE_matches, EVE_SIMILARITY = as.numeric(EVE_SIMILARITY))
EVE_matches$EVE_BITSCORE <- sub("^", "-", EVE_matches$EVE_BITSCORE)
EVE_matches <- transform(EVE_matches, EVE_BITSCORE = as.numeric(EVE_BITSCORE))
#EXO dataframe
#EXO_matches <- read.csv("EXO_genomic_matches.csv", header = FALSE)
EXO_matches <- read.csv("EXO_gen_matches.tsv", header = FALSE, sep = "\t")
EXO_matches = subset(EXO_matches, select = -c(V4,V5,V6,V7,V8,V9,V11))
#colnames(EXO_matches) <- c("SAMPLE", "EXO","EXO_SIMILARITY", "EXO_BITSCORE", "FRAMESHIFT_EXO")
colnames(EXO_matches) <- c("SAMPLE", "EXO","EXO_SIMILARITY", "EXO_BITSCORE", "FRAMESHIFT_EXO")
# merging dataframes
matches_table <- merge(EXO_matches, EVE_matches, by = 'SAMPLE', all.x = TRUE, all.y = TRUE ) 
matches_table <- matches_table%>% replace(is.na(.), 0)
matches_table$DIFF_SIMILARITY <- matches_table$EXO_SIMILARITY+(matches_table$EVE_SIMILARITY)
matches_table$DIFF_BITSCORE <- matches_table$EXO_BITSCORE+(matches_table$EVE_BITSCORE)
#write.table(matches_table, file = "EVE_THRESHOLDS.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(matches_table, file = "gen.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

#matches_table <- read.table("EVE_THRESHOLDS_NEW", header = TRUE, sep = "\t")
genomic <- ggplot(matches_table) + 
  #geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color=EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY), shape=FRAMESHIFT_EXO), size=4) +
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color=EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY), shape=FRAMESHIFT_EXO), size=4) +
  scale_color_viridis_c(option = "plasma", name = "identity to known EVE", breaks=c(0.0,0.2,0.4,0.6,0.8),labels=c("0 (min)",0.2,0.4,0.6,"0.8 (max)"),limits=c(0,0.8)) +
  scale_y_log10(breaks=c(100, 250, 500,1000, 1500)) +
  #xlim(-500,1800) +
  scale_x_continuous(breaks = c(-500, -250,-100,0,100, 250, 500, 1000, 1500),limits=c(-500, 1800)) +
  ggtitle("B) EVEs in reference genome") +
  scale_shape_manual(values=c(20, 24), name = "frameshifts present")+
  theme(legend.position = c(0.7, 0.3),
        legend.direction = "vertical",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key=element_blank(),
        legend.title.align = 1,
        legend.key.size = unit(1, 'cm'),
        #legend.direction = legend.direction,
        #legend.box.just = "center",
        legend.background=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "white"), 
        plot.title = element_text(size = 22, hjust = 0.5),
        axis.line.x = element_line(color="black", linewidth = 1),
        axis.line.y = element_line(color="black", linewidth = 1)) +
  xlab("?? BLAST bit score") +# for the x axis label
  ylab("maximum BLAST bit score per contig (log10)") 
#  guide_legend(title.position = "top")
genomic

#### transcirptomic ####

#EVE_dataframe
#EVE_matches <- read.csv("EVE_trans_matches.csv", header = FALSE)
EVE_matches <- read.csv("EVE_trans_matches.tsv", header = FALSE, sep = "\t")
EVE_matches = subset(EVE_matches, select = -c(V4,V5,V6,V7,V8,V9,V11))
colnames(EVE_matches) <- c("SAMPLE", "EVE","EVE_SIMILARITY", "EVE_BITSCORE", "FRAMESHIFT_EVE")
#colnames(EVE_matches) <- c("SAMPLE", "EVE","EVE_SIMILARITY", "EVE_BITSCORE", "FRAMESHIFT_EVE")
EVE_matches$EVE_SIMILARITY <- sub("^", "-", EVE_matches$EVE_SIMILARITY)
EVE_matches <- transform(EVE_matches, EVE_SIMILARITY = as.numeric(EVE_SIMILARITY))
EVE_matches$EVE_BITSCORE <- sub("^", "-", EVE_matches$EVE_BITSCORE)
EVE_matches <- transform(EVE_matches, EVE_BITSCORE = as.numeric(EVE_BITSCORE))
#EXO dataframe
#EXO_matches <- read.csv("EXO_trans_matches.csv", header = FALSE)
EXO_matches <- read.csv("EXO_trans_matches.tsv", header = FALSE, sep = "\t")
EXO_matches = subset(EXO_matches, select = -c(V4,V5,V6,V7,V8,V9,V11))
#colnames(EXO_matches) <- c("SAMPLE", "EXO","EXO_SIMILARITY", "EXO_BITSCORE", "FRAMESHIFT_EXO")
colnames(EXO_matches) <- c("SAMPLE", "EXO","EXO_SIMILARITY", "EXO_BITSCORE", "FRAMESHIFT_EXO")
# merging dataframes
matches_table <- merge(EXO_matches, EVE_matches, by = 'SAMPLE', all.x = TRUE, all.y = TRUE ) 
matches_table <- matches_table%>% replace(is.na(.), 0)
matches_table$DIFF_SIMILARITY <- matches_table$EXO_SIMILARITY+(matches_table$EVE_SIMILARITY)
matches_table$DIFF_BITSCORE <- matches_table$EXO_BITSCORE+(matches_table$EVE_BITSCORE)
#write.table(matches_table, file = "EVE_THRESHOLDS_TRANS.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(matches_table, file = "trans.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#matches_table <- read.table("EVE_THRESHOLDS_NEW", header = TRUE, sep = "\t")



transcriptomic <- ggplot(matches_table) + 
  #geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color=EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY), shape=FRAMESHIFT_EXO), size=4) +
  geom_jitter(aes(DIFF_BITSCORE, pmax(EXO_BITSCORE, -EVE_BITSCORE), color=EVE_SIMILARITY/(EVE_SIMILARITY-EXO_SIMILARITY), shape=FRAMESHIFT_EXO), size=4) +
  #scale_color_viridis_c(option = "plasma", name = "identity to known EVE", breaks=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),labels=c("0 (min)",0.1,0.2,0.3,0.4,0.5,0.6,"0.7 (max)"),limits=c(0,0.7)) +
  scale_color_viridis_c(option = "plasma", name = "identity to known EVE", breaks=c(0.0,0.2,0.4,0.6,0.8),labels=c("0 (min)",0.2,0.4,0.6,"0.8 (max)"),limits=c(0,0.8)) +
  scale_y_log10(breaks=c(100, 250, 500,1000, 1500)) +
  #xlim(-500,1800) +
  scale_x_continuous(breaks = c(-500, -250,-100,0,100, 250, 500, 1000, 1500), limits=c(-500, 1800)) +
  scale_shape_manual(values=c(20, 24), name = "frameshifts present")+
  ggtitle("A) No EVEs in reference genome") +
  theme(legend.position = c(0.7, 0.3),
        legend.direction = "vertical",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.5, 'cm'),
        legend.key=element_blank(),
        legend.title.align = 1,
        legend.key.size = unit(1, 'cm'),
        #legend.direction = legend.direction,
        #legend.box.just = "center",
        legend.background=element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        panel.background = element_rect(fill = "white"), 
        axis.line.x = element_line(color="black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.line.y = element_line(color="black", linewidth = 1)) +
  xlab("?? BLAST bit score") +# for the x axis label +
  ylab("maximum BLAST bit score per contig (log10)") 

show(transcriptomic)

library(patchwork)
transcriptomic + genomic

#  guide_legend(title.position = "top")


#### statistical test 


EVE_matches <- read.csv("EVE_trans_matches.tsv", header = FALSE, sep = "\t")
EVE_matches = subset(EVE_matches, select = -c(V5,V6,V7,V8,V9,V11))
colnames(EVE_matches) <- c("SAMPLE", "EVE","EVE_SIMILARITY", "EVE_LENGTH", "EVE_BITSCORE", "FRAMESHIFT_EVE")
#EVE_matches$EVE_SIMILARITY <- sub("^", "-", EVE_matches$EVE_SIMILARITY)
#EVE_matches <- transform(EVE_matches, EVE_SIMILARITY = as.numeric(EVE_SIMILARITY))
EVE_matches$EVE_BITSCORE <- sub("^", "-", EVE_matches$EVE_BITSCORE)
EVE_matches <- transform(EVE_matches, EVE_BITSCORE = as.numeric(EVE_BITSCORE))
#EXO dataframe

EXO_matches <- read.csv("EXO_trans_matches.tsv", header = FALSE, sep = "\t")
EXO_matches = subset(EXO_matches, select = -c(V5,V6,V7,V8,V9,V11))
colnames(EXO_matches) <- c("SAMPLE", "EXO","EXO_SIMILARITY", "EXO_LENGTH", "EXO_BITSCORE", "FRAMESHIFT_EXO")
# merging dataframes
matches_table <- merge(EXO_matches, EVE_matches, by = 'SAMPLE', all.x = TRUE, all.y = TRUE ) 
matches_table <- matches_table%>% replace(is.na(.), 0)
matches_table$DIFF_SIMILARITY <- matches_table$EXO_SIMILARITY+(matches_table$EVE_SIMILARITY)
matches_table$DIFF_BITSCORE <- matches_table$EXO_BITSCORE+(matches_table$EVE_BITSCORE)

library(dplyr)

# Assign EVE or exogenous dependent on diff bitscore 
matches_table <- matches_table %>%
  mutate(
    GROUP = ifelse(DIFF_BITSCORE <= 0, "EVE", "exogenous")
  )

# Remove all rows where GROUP is "EVE"
matches_table <- matches_table %>%
  filter(GROUP != "EVE")

# Test 1: Amino Acid Similarities (EVE_SIMILARITY vs. EXO_SIMILARITY)
result_similarity <- wilcox.test(matches_table$EVE_SIMILARITY, matches_table$EXO_SIMILARITY, paired = TRUE)

# Test 2: Sequence Lengths (EVE_LENGTH vs. EXO_LENGTH)
result_length <- wilcox.test(matches_table$EVE_LENGTH, matches_table$EXO_LENGTH, paired = TRUE)

# Combine the p-values using Fisher's method
combined_p_value <- 1 - pchisq(-2 * sum(log(c(result_similarity$p.value, result_length$p.value))), df = 2 * 2)

# Print the test results
print(result_similarity)
print(result_length)

