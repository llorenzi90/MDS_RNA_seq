#generate mid-point (TSS-like) file for ATAC-seq peaks
options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
## ---------------------------

#load DESeq2 results
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/all_samples/cohesin_vs_nocohesin/")
deseq2res=read.csv("all_samples_cohesin_vs_nocohesin_DESeq2_results.csv",header = T)
padjcutoff=0.05
upPeaks <- deseq2res %>% filter(padj<=padjcutoff , !is.na(padj) ,log2FoldChange>0)
downPeaks <- deseq2res %>% filter(padj<=padjcutoff , !is.na(padj) ,log2FoldChange<0)

pl <- list(UPpeaks=upPeaks,
           DOWNpeaks=downPeaks)
for (nam in names(pl)) {
  pp=pl[[nam]]
  chrs <- sapply(strsplit(pp$X,":"),function(x)x[1])
  coords <- strsplit(sapply(strsplit(pp$X,":"),function(x)x[2]),
                     split = "-")
  midpoints <- sapply(coords, function(x){
    x=as.numeric(x)
    wid=x[2] - x[1]+1
    return(x[1] + round(wid/2)-1)
  }
                      )
  dftw <- data.frame(chr=chrs,
                     midpoints=midpoints,
                     id=pp$X) 
  write.table(dftw,paste0(nam,"_midpoints_coords.txt"),quote = F,row.names = F, sep = "\t")
}

#test length of fragments before and after midpoint
table(sapply(1:nrow(dftw),function(i){
  length(as.numeric(coords[[i]][1]):dftw$midpoints[i]) - length((dftw$midpoints[i]+1):as.numeric(coords[[i]][2]))
}))
