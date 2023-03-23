#!/usr/bin/Rscript

### SEQUENCING DEPTH ###

setwd("/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/Aedes/diamond_blastx/Aedes/diamond_blast/EVE_ANALYSIS/R_folder")

# Sequence depth
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('stringr')) install.packages('stringr'); library('stringr')
#library(ggplot2)
#library(stringr)
# create data table
data.cov <- read.table(list.files(pattern =".cov"), header=FALSE, na.strings = "NA")
colnames(data.cov) = c("ref_name","position","depth")

# create ggplot
ggplot(data=data.cov, aes(x=position, y=log10(depth))) + geom_line() +
  theme_minimal() + theme(text = element_text(size=15)) +
  ylab(label="Sequencing depth (log10 scale)") +
  xlab(label="Position on reference genome")

# save plot in RESULTS
file=word(list.files(pattern =".cov"))
ggsave(file,device = "png", path = "/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/Aedes/diamond_blastx/Aedes/diamond_blast/EVE_ANALYSIS/COVERAGE/PROCESSED/DEPTH_VARIANTS")
# delete .cov file to process next one
file.remove(list.files(pattern =".cov"))

### LOFREQ ###
if (!require('readr')) install.packages('readr'); library('readr')
if (!require('plyr')) install.packages('plyr'); library('plyr')
library(readr)
library(stringr)
library(plyr)
vcf <- read_delim(list.files(pattern =".vcf"), "\t", escape_double = FALSE, 
    col_names = FALSE, trim_ws = TRUE, skip = 18, col_types = cols())
colnames(vcf) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

vcf <- ddply(vcf, .(POS), mutate, DP = word(word(INFO,1,sep=";"),2, sep="="),
             AF = as.numeric(word(word(INFO,2,sep=";"),2, sep="=")),
             SB = as.numeric(word(word(INFO,3,sep=";"),2, sep="=")),
             DP4_refFOR = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),1,sep=",")),
             DP4_refREV = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),2,sep=",")),
             DP4_altFOR = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),3,sep=",")),
             DP4_altREV = as.numeric(word(word(word(INFO,4,sep=";"),2, sep="="),4,sep=",")))
vcf$INFO <- NULL

write.csv(vcf,file = word(list.files(pattern ="vcf")) , row.names = FALSE)

library(ggplot2)
ggplot(data=vcf, aes(x=POS, y=AF)) + geom_point() + 
  geom_hline(yintercept = 0.5, color="red", linetype="dotted") +
  theme_minimal() + theme(text = element_text(size=15)) +
  ylab(label="Variant frequency") +
  xlab(label="Position on reference genome")
 
file=word(list.files(pattern =".vcf"), sep = ".vcf")
ggsave(file,device = "png", path = "/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/Aedes/diamond_blastx/Aedes/diamond_blast/EVE_ANALYSIS/COVERAGE/PROCESSED/DEPTH_VARIANTS")
