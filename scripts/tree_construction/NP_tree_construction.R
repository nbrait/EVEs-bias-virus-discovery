setwd("C:/Users/nadja/Documents/R/trees_for_quaranja_new/TREES") # state working directory

##################################################
# Dependencies
##################################################
#install.packages("ggplot2")
library (ggplot2)
library(ggrepel)
library (phytools)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree", force = TRUE) 
library(BiocManager)
#install("ggtree", force = TRUE)
#creating tree
library(ggtree)
library(tidytree)
library(ape)
library(treeio)
##################################################

# following files are needed for this workflow:
# contree file from iqtree (or other nex/nwk tree files- with bootstrap values)
# csv file with samples and their country/continent/file origin - acquired from sra metadata
# csv file with reference accession numbers and scientific names 

## reading in tree
tree <- read.tree("NP_tree.contree") # created with iqtree
## replacing the label descriptions
replace_words <- read.csv("NP_references.txt", sep = "\t", header = FALSE, col.names = c("Accession", "Virus"))
for(j in seq_along(replace_words$Accession)){
  tree$tip.label <- gsub(replace_words$Accession[j], replace_words$Virus[j], tree$tip.label)
}
tree <- as.phylo(tree)
#outgroup rooting
nodes <- as_tibble(tree) # check parent node for Isavirus 
tree <- root(tree, node = 323, edgelabel = TRUE)
#### converting the bootstrap value outside ggtree
q <- ggtree(tree) + xlim(NA, 7.5) + expand_limits(y = 130) + coord_cartesian(clip="off") #+ geom_text2(aes(label=node), hjust=-.3, size=2)  
q <- flip(q, 330, 493)
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,] # for significant bootstrap values
d$label <- replace(d$label, d$label > 94,'*') 
q <- q + geom_text(data=d, aes(label=label), hjust = 0.0, vjust = 0.1, size = 3.0) 

## annotation location and samples
tipcategories = read.csv("location_type_NP_new.txt", 
                         sep = "\t",
                         col.names = c("seq", "Origin", "file_old", "file" ), 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

dd = as.data.frame(tipcategories)

q <- q %<+% dd +
  #geom_tiplab(aes(), hjust = -0.2, size=2.5, alpha=.75) +
  geom_cladelab(node=633, label="Influenzavirus", align=FALSE,  
                offset = .1, textcolor='black', barcolor='black', angle=0, vjust=0.5, hjust=0, barsize = 1,fontface = 2) +
  geom_cladelab(node=636, label="Thogotovirus", align=FALSE,  
                offset = .1, textcolor='black', barcolor='black', angle=0, vjust=0.5, hjust=0, barsize = 1,fontface = 2) +
  geom_cladelab(node=626, label="Quaranjavirus (7 segments)", align=FALSE,  
                offset = .1, textcolor='black', barcolor='black', angle=0, vjust=0.5, hjust=0, barsize = 1,fontface = 2) +
  geom_cladelab(node=329, label="Quaranjavirus (8 segments)", align=FALSE,  
                offset = -0.5, textcolor='black', barcolor='black', angle=90, vjust=1, hjust=0.5, barsize = 1, fontface = 2) +
  geom_cladelab(node=571, label="Aedes detritus \n orthomyxo-like virus", align=FALSE,  
                offset = .5, textcolor='black', barcolor='black', angle=90, vjust=1, hjust=0.5) +
  geom_cladelab(node=375, label="Usinis virus", align=FALSE,  
                offset = .25, textcolor='black', barcolor='black', angle=90, vjust=0.5, hjust=0.5) +
  geom_cladelab(node=74, label="HP RP20_CCG015879 - ", align=FALSE,  
                offset = -0.3, offset.text= -0.7, textcolor='dark red', barcolor='dark red', fontsize = 3) +
  geom_cladelab(node=441, label="Wuhan Mosquito virus 4", align=FALSE,  
                offset = .15, textcolor='black', barcolor='black', angle=90, vjust=0.5, hjust=0.5) +
  geom_cladelab(node=458, label="	Culex \n orthomyxo-like virus ", align=FALSE,  
                offset = .3, textcolor='black', barcolor='black', angle=90, vjust=1, hjust=0.5) +
  geom_cladelab(node=335, label="	Wuhan mosquito virus 6", align=FALSE,  
                offset = 0.05, textcolor='black', barcolor='black', angle=90, vjust=0.5, hjust=0.5) +
  geom_cladelab(node=526, label="Guadeloupe mosquito \n quaranja-like virus 1", align=FALSE,  
                offset = .1, textcolor='black', barcolor='black', angle=90, vjust=1, hjust=0.5) +
  geom_cladelab(node=1, label="Isavirus", align=FALSE,  
                offset = .1, textcolor='black', barcolor='black',fontface = 2) +
  geom_cladelab(node=2, label="Sardinovirus", align=FALSE,  
                offset = .4, textcolor='black', barcolor='black',fontface = 2) +
  geom_cladelab(node=497, label="Aedes alboannulatus \n orthomyxo-like virus ", align=FALSE,  
                offset = .3, textcolor='black', barcolor='black', angle=90, vjust=1, hjust=0.5) +
  geom_cladelab(node=295, label="HP RP20_CCG016659 - ", align=FALSE,  
                offset = -0.15, offset.text = -0.7, textcolor='dark red', barcolor='dark red', fontsize = 3) +
  geom_cladelab(node=618, label="Astopletus virus ", align=FALSE,  
                offset = .3, textcolor='black', barcolor='black', angle=90, vjust=0.5, hjust=0.5) +
  geom_cladelab(node=224, label="SRR6155880 -  ", align=FALSE,  
                offset = -0.5, offset.text = 0, textcolor='dark blue', barcolor='white', angle=0, fontsize = 3) +
  geom_cladelab(node=229, label="SRR6155879 -  ", align=FALSE,  
                offset = -0.5, offset.text = 0, textcolor='dark blue', barcolor='white', angle=0, fontsize = 3) +
  geom_tippoint(aes(color=file), size=1.5, alpha=.75) +
  scale_color_manual(values = c("EVE" = "darkred", "SRR" = "dark blue", "EVET"="orange"), na.translate=F, name = "", labels = c("endogenous viral elements (EVE)", "detected exogenous viral segments", "EVE + transcriptomic detection")) +
  theme(legend.position= c(0.15, 0.2), legend.text = element_text(size = 15), title = element_text(size = 22)) +
  theme(legend.key=element_blank(),legend.background=element_blank()) 
show(q)

### closer look at specific clade for checking
#viewClade(tree_view = NULL, 335, xmax_adjust = 0)
#show (q)
