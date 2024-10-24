## ---------------------------
##
##
## Purpose of script: volcano plot ATAC-seq MDS
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-14
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
require(ggplot2)
## ---------------------------
#read data
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/all_samples/cohesin_vs_nocohesin/")
deseq2res <- read.csv("all_samples_cohesin_vs_nocohesin_DESeq2_results.csv")
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
padjco=0.05
deseq2res <- deseq2res%>%mutate(DIFF=ifelse(padj<=padjco&!is.na(padj),
                                            ifelse(log2FoldChange>0,"UP","DOWN"),
                                            "NO"))
table(deseq2res$DIFF)

pdf("voclano_plot.pdf")
ggplot(deseq2res,aes(x=log2FoldChange,y=-log10(padj))) +geom_point(aes(col=DIFF))+
  theme_classic() + scale_color_manual(values=cbPalette[c(6,1,7)])
dev.off()
