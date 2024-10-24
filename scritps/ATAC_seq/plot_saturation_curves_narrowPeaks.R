library(ggplot2)
library(ATACseqQC)

samp=commandArgs(trailingOnly = T)[1]
datadir=paste0("/Users/llorenzi/MDS/ATAC-seq/data/macs2_peaks/",samp,"/")
setwd("/Users/llorenzi/MDS/ATAC-seq/analyses/QC/saturation_curves")
peakFiles <- dir(datadir, "narrowPeak$")

summary_file=read.csv("/Users/llorenzi/MDS/ATAC-seq/analyses/QC/sample_stats/summary_QC_after_alignment.csv")
total_frags <- summary_file$primarychr_MAPQ20_proper_pairs[summary_file$sample_id==samp]

subsamplingFractions <- gsub("(AT-.*.name.sorted.)(0.[0-9]{1,2})(.coord.sorted_peaks.narrowPeak)","\\2",peakFiles)
subsamplingFractions[grep("primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak",subsamplingFractions)] <- 1
subsamplingSizes <- as.numeric(subsamplingFractions)*total_frags

names(subsamplingSizes) <- paste0(datadir,peakFiles)

peakStat <- saturationPlot(subsamplingPeakFiles=paste0(datadir,peakFiles), subsamplingSizes=subsamplingSizes, sep="\t", 
                           header= FALSE, fdr=0.05, fdrCol=9, startCol=2, 
                           endCol=3, skipLines=1, peakCaller="MACS2", outPrefix=samp)

pdf(paste0("Saturation curve for ",samp,".pdf"), width=8, height=8)

## total peak number-based curve
p <-ggplot(data= peakStat, aes(x= subsamplingSizes/10^6, y = numPeaks/10^3)) + geom_point(shape=16, size=3) + geom_smooth()
p + expand_limits(y=0) + theme(text=element_text(size=10)) + xlab(expression(Effective~fragments~x~10^6)) + ylab(expression(Peaks~x~10^3)) 

## total peak width-based 
p <-ggplot(data= peakStat, aes(x= subsamplingSizes/10^6, y = breadth/10^6)) + geom_point(shape=16, size=3) + geom_smooth()
p + expand_limits(y=0) + theme(text=element_text(size=10)) + xlab(expression(Effective~fragments~x~10^6)) + ylab(expression(Total~peaks~width~(Mb))) 

dev.off()
