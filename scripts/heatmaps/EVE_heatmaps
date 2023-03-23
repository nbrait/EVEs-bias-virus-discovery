## Dependencies 
#####################################################
#install.packages("plyr")
#install.packages("ggplot2")
library(plyr)
library(ggplot2)
library(tidyr)
library(patchwork)
#####################################################
setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/EVE_COVERAGE")
table <- read.csv("genomic_EVES_count.csv", header = TRUE, sep = "\t")
############### geom_tile heatmap ###################

table$EVE_Segment <- factor(table$EVE_Segment, levels = c("PB1","PB2","PA","NP","GP","HP1","HP2","HP3"))
table$WGS <- factor(table$WGS, levels = c("AAWU01","JADCQB01","JADDTP01","AAGE02","JABWQA01","JAFDOQ01","JAJPPM01","JXPU01", "JXUM01","LMAV01","MNAF02","NIGP01","SWKZ01"))

average.heatmap <- ggplot(data = table, mapping = aes(x = EVE_Segment,
                                                      y = WGS,
                                                      fill = Average)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label=format(round(Average,0))), colour="white", size = 5) +
  scale_fill_gradient(low = "white", high = "#006d2c", name = "") +
  coord_equal() +
  theme(panel.background = element_rect(fill = "white"),text = element_text(size = 15,face = "bold"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title.align=0.2, legend.position = "bottom") +
  ggtitle(label = "Average length coverage (%)") 

largest.heatmap <- ggplot(data = table, mapping = aes(x = EVE_Segment,
                                                      y = WGS,
                                                      fill = Largest)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label=format(round(Largest,0))), colour="white", size = 5) +
  scale_fill_gradient(low = "white", high = "#cc4c02", name = "") +
  coord_equal() +
  theme(panel.background = element_rect(fill = "white"),text = element_text(size = 15,face = "bold"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title.align=0.2, legend.position = "bottom") +
  ggtitle(label = "Maximum length coverage (%)") 

count.heatmap <- ggplot(data = table, mapping = aes(x = EVE_Segment,
                                                    y = WGS,
                                                    fill = Count)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label=Count), colour="black", size = 5) +
  ylab(label = "Segment") + xlab(label = "WGS sample") +
  scale_fill_gradient2(low = "white", high = "#6a51a3", name = "") +
  coord_equal() +
  theme(panel.background = element_rect(fill = "white"),text = element_text(size = 15,face = "bold"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90),axis.title.x = element_blank(), axis.title.y = element_blank(),  legend.title.align=0.2, legend.position = "bottom") +
  ggtitle(label = "EVE count")

all.heatmap <- count.heatmap+average.heatmap+largest.heatmap

all.heatmap
