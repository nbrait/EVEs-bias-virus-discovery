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
tree <- read.tree("HP2_tree.contree") # insert nwk tree with bootstrap
## replacing the label descriptions
replace_words <- read.csv("HP2_references.txt", sep = "\t", header = FALSE, col.names = c("Accession", "Virus"))
for(j in seq_along(replace_words$Accession)){
  tree$tip.label <- gsub(replace_words$Accession[j], replace_words$Virus[j], tree$tip.label)
}

tree <- as.phylo(tree)
tree <- midpoint.root(tree)
#outgroup rooting
#nodes <- as_tibble(tree) 
#tree <- root(tree, node = 68, edgelabel = TRUE)
#### converting the bootstrap value outside ggtree
q <- ggtree(tree) + xlim(NA, 5) + expand_limits(y = 70) + coord_cartesian(clip="off") #+ geom_text2(aes(label=node), hjust=-.3, size=2) 
show(q)
#q <- flip(q, 330, 493)
#q <- rotate(q,138) %>% rotate(150)
d <- q$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d <- d[d$label > 95,]
d$label <- replace(d$label, d$label > 95,'*') 
q <- q + geom_text(data=d, aes(label=label), hjust = 0.0, vjust = 0.1, size = 3.0) + geom_treescale(x=0, y=3) 

## annotation location and samples
tipcategories = read.csv("location_type_NP_new.txt", 
                         sep = "\t",
                         col.names = c("seq", "Origin", "file_old", "file" ), 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

dd = as.data.frame(tipcategories)

q <- q %<+% dd +
  #geom_tiplab(aes(),hjust = -0.2, size=2.5, alpha=.75) +
  geom_cladelab(node=171, label="Usinis virus", align=FALSE,  
                offset = .25, textcolor='black', barcolor='black', angle=0, vjust=0.5, hjust=0) +
  geom_cladelab(node=1, label="Byreldi virus", align=FALSE,  
                 offset = .25, textcolor='black', barcolor='black') +
   geom_cladelab(node=178, label="Byreska virus", align=FALSE,  
                 offset = .3, textcolor='black', barcolor='black', angle=90, vjust=0.5, hjust=0.5) +
  geom_cladelab(node=144, label="	Wuhan mosquito virus 6", align=FALSE,  
                offset = 0.05, textcolor='black', barcolor='black', angle=90, vjust=0.5, hjust=0.5) +
  geom_cladelab(node=100, label="Guadeloupe mosquito \n quaranja-like virus 1", align=FALSE,  
                offset = .1, textcolor='black', barcolor='black', angle=90, vjust=1, hjust=0.5) +
  geom_cladelab(node=2, label="Astopletus virus ", align=FALSE,  
                offset = .3, textcolor='black', barcolor='black', angle=0, vjust=0.5, hjust=0) +
  # #geom_strip(176, 165, color='black',align=FALSE,offset = .6, label="unclassified \n Orthomyxoviridae", angle=90, vjust=1, hjust=0.5) +
  # #geom_strip(257, 262, color='black',align=FALSE,offset = .2, label="Culex pipiens orthomyxo-like virus",fontsize = 6) +
  geom_tippoint(aes(color=file), size=1.5, alpha=.75) +
  scale_color_manual(values = c("EVE" = "darkred", "SRR" = "blue", "EVET"="orange"), na.translate=F, name = "", labels = c("endogenous viral elements (EVE)", "detected exogenous viral segments", "EVE + transcriptomic detection")) +
  theme(legend.position= c(0.25, 0.9), legend.text = element_text(size = 15), title = element_text(size = 22)) +
  theme(legend.key=element_blank(),legend.background=element_blank()) 
show(q)
#viewClade(tree_view = NULL, 196, xmax_adjust = 0)

