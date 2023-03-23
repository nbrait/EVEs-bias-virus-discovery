##!/usr/bin/Rscript

library(readr)
library(ggplot2);library(gridExtra)
library(dplyr)
library(reshape2)
library(gtools)
library(stringr)
library(data.table)
library(seqinr)

setwd("/mnt/c/Users/nadja/Documents/LaptopAsus/PhD/Quaranjaproject/Aedes/diamond_blastx/Aedes/diamond_blast/EVE_ANALYSIS/COVERAGE/PROCESSED/CORRECTION")

listsample <- read_csv("listsample.txt", col_names = FALSE)

listsample = as.vector(listsample$X1)

for(i in listsample){
  #loading chimera
  chimera = read.fasta(paste(i,".fa",sep=""),seqtype = c("DNA"),forceDNAtolower = FALSE)
 
  #loading variants
  data.var <- read_delim(paste(i,".vcf",sep=""), "\t", escape_double = FALSE, trim_ws = TRUE, skip = 17)
  colnames(data.var) = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  if(nrow(data.var) != 0){
  data.var = data.var %>% mutate(AF=as.numeric(word(word(INFO,2,sep=";"),2,sep="=")))
  data.var.consensus = subset(data.var, AF > 0.5)
  }else(data.var.consensus = data.var)
  
  #loading coverage
  data.cov <- read.table(paste(i,".cov",sep=""), header=F, na.strings = "NA")
  colnames(data.cov) = c("sample","POS","X")
  data.cov <- dplyr::filter(data.cov, grepl(i,sample))
  data.cov.tolow = subset(data.cov, X < 3) # change according to how stringent you want your cut-off!

  #coverage correction
  if(nrow(data.cov.tolow)!= 0){
  for(k in data.cov.tolow$POS){
    chimera[[paste(i)]][k] = "N"
  }}
 
  #sequence correction
  if(nrow(data.var.consensus) != 0){
  for(j in 1:nrow(data.var.consensus)){
    chimera[[paste(i)]][data.var.consensus[j,]$POS] = data.var.consensus[j,]$ALT
  }}

  write.fasta(chimera, paste(i), paste(i,"-3_corrected.fasta",sep=""), open = "w", nbchar = 60, as.string = FALSE)
} 
