#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("forcats")
library(plyr)
library(ggplot2)
library(forcats)
setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/HEATMAP")

## coverage calculations 
#table <- read.csv("NP_all_segments.csv")
#Xlength <- lengths(regmatches(table$Sequence, gregexpr("N", table$Sequence)))
#totallength <- nchar(table$Sequence)
#Percentage <- 100-(Xlength/totallength*100)
#table$Percentage <- Percentage
#table <- table[order(table$Percentage),]
#write.csv(table, file = "tableforcoverage.csv")

#colnames(table)[3] <- "HP3"
#table= subset(table, select = -c(Sequence) )
#table <- table[order(table$Name),]
#table$Name <- ave(table$Name, table$Name, FUN = function(i) paste0(i, '.', seq_along(i)))
#row.names(table) <- table$Name

#table3 <- table

#table3 <- rbind.fill(table3,table)
#table3 <- merge(table,table3, all=TRUE)
#table3 <- table3[order(table3$Name),]
#table2 <- table3
#table3 <- table2
# changing NAs to 0
#table3[is.na(table3)] <- 0

############### geom_tile heatmap ###################

table3 <- read.csv("Aedini_new.tsv", sep = "\t", header = TRUE)
table3 <- table3[order(as.character(table3$Sample.1)),]
table3$Segment <- factor(table3$Segment, levels = c("PB1","PB2","PA","NP","GP","HP1","HP2","HP3"))

Aedes.heatmap <- ggplot(data = table3, mapping = aes(x = Segment,
                                                       y = fct_rev(as_factor(Sample.1)),
                                                       fill = Percentage)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  #ylab(label = "Segment") + xlab(label = "SRA sample") +
  scale_fill_gradient(low = "#d9d9d9", high = "#cb181d", na.value = 'white', name = "average cov\nlength %" ) +
  # keep them squared
  coord_fixed(0.5) +
  scale_x_discrete(position = "top") +
  #scale_y_continuous(minor_breaks = c(1.7, 2.3, 4.1)) +
  # change the width and the height of the legend color bar
  #guides(fill = guide_colourbar(barwidth = 0.5,
                                #barheight = 5, title = "Coverage %")) +
  #theme_bw() 
  theme(axis.text.y = element_text(hjust = 0), panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y =element_blank()) 
  #ggtitle(label = "Orthomyxoviridae Segments") 
  #scale_fill_continuous(na.value = 'salmon')
  
# add the values over the tiles
# geom_text(aes(label = value), color = "white", size = 4) +

Aedes.heatmap

#Culex table
table3 <- read.csv("Culicini.tsv", sep = "\t", header = TRUE)
table3 <- table3[order(as.character(table3$Sample.1)),]
table3$Segment <- factor(table3$Segment, levels = c("PB1","PB2","PA","NP","GP","HP1","HP2","HP3"))
#write.table(table3, file = "test2.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
Culex.heatmap <- ggplot(data = table3, mapping = aes(x = Segment,
                                                     y = fct_rev(as_factor(Sample.1)),
                                                     fill = Percentage)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  #ylab(label = "Segment") + xlab(label = "SRA sample") +
  scale_fill_gradient(low = "#d9d9d9", high = "#2171b5", na.value = 'white', name = "average cov\nlength %" ) +
  # keep them squared
  coord_fixed(0.5) +
  scale_x_discrete(position = "top") +
  # change the width and the height of the legend color bar
  #guides(fill = guide_colourbar(barwidth = 0.5,
  #barheight = 5, title = "Coverage %")) +
  #theme_bw() +
  theme(axis.text.y = element_text(hjust = 0), panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y =element_blank()) +
  ggtitle(label = "Orthomyxoviridae Segments") 
#scale_fill_continuous(na.value = 'salmon')

# add the values over the tiles
# geom_text(aes(label = value), color = "white", size = 4) +

Culex.heatmap

#Other table
table3 <- read.csv("Culicinae.tsv", sep = "\t", header = TRUE)
table3 <- table3[order(as.character(table3$Sample.1)),]
table3$Segment <- factor(table3$Segment, levels = c("PB1","PB2","PA","NP","GP","HP1","HP2","HP3"))
#write.table(table3, file = "test1.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
other.heatmap <- ggplot(data = table3, mapping = aes(x = Segment,
                                                     y = fct_rev(as_factor(Sample.1)),
                                                     fill = Percentage)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradient(low = "#d9d9d9", high = "purple4", na.value = 'white', name = "average cov\nlength %" ) +
  coord_fixed(0.5) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 0),panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90), axis.title.x=element_blank(), axis.title.y =element_blank()) +
  ggtitle(label = "Orthomyxoviridae Segments") 
other.heatmap

#Aedes.heatmap + Culex.heatmap + other.heatmap
