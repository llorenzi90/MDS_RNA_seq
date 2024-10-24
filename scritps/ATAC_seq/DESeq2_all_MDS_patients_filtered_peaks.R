## ---------------------------
##
##
## Purpose of script: MDS DESeq2 LPS vs no-LPS, filtering low confidence peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-18
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: I will try different filtering cutoffs/criteria
##   
##
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## ---------------------------
require(tidyverse)
require(data.table)
require(DESeq2)
library("dplyr")
library("ggplot2")
library(RColorBrewer)
library(pheatmap)
## ---------------------------
countdata=read.table("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/merged_macs2peaksNoDups_MDS_all_samples.counts")
countdata=read.table("/home/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/merged_macs2peaksNoDups_MDS_all_samples.counts")

samples=colnames(countdata)
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1])
coldata=data.frame(samples,patient,LPS)

rownames(coldata) <- coldata$samples
table(rownames(coldata)==colnames(countdata))
coldata$LPS <- as.factor(coldata$LPS)

#Add data on mutations

mut_data <- fread("~/Nextcloud/MDS/mut_data_all_patients.txt")
#peaks=read.table("~/Nextcloud/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames",header = T)
peaks <- read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames',header = T )
peaks$peakID=paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)


Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]
coldata$Cohesin=coldata$patient%in%mut_data$Sample[mut_data$Pathway=="Cohesin"]
table(coldata$Cohesin)

#overview n samples per peak
total_peaks=nrow(peaks)
summary(peaks$nsamples)
#On averge each peak is called in 4.5 samples
table(peaks$nsamples)
table(peaks$nsamples)/total_peaks*100
#63% of the peaks were called in only one sample
#only 18 peaks are called in all samples
table(peaks$nsamples>=42)
table(peaks$nsamples>=42)/total_peaks*100
#  12048 peaks (2.7%) are called in at least half of the samples


#filtering strategies:
#1) peaks that are in half of the samples
#2) peaks that are in half of the samples + peaks that have a very low q-val
#NOTE: the column "meanscore" in the peaks file is not correctly named
#this value is the average of the (-log10) q-val across all samples
# that contribute to that peak, this q-val is retrieved as column 9 in narrowPeak macs2 file
max(10^(-peaks$meanscore))
#[1] 0.04999539
#this is ok, as I used the default cutoff for q-val ,i.e: 0.05

summary(10^(-peaks$meanscore))
peaks$qval=10^(-peaks$meanscore)

#check distribution of q-vals for those peaks that are called in only 1 sample
summary(peaks$qval[peaks$nsamples==1])
plot(density(peaks$qval[peaks$nsamples==1]))
points(density(peaks$qval[peaks$nsamples>=10]),col="red",type = "s")
plot(density(peaks$qval[peaks$nsamples>=10]))
summary(peaks$qval[peaks$nsamples>=10])

table(peaks$qval<0.0001&peaks$nsamples==1)

peaks$nsamps_factor=1
peaks$nsamps_factor[peaks$nsamples<=14] <- "1-14"  
peaks$nsamps_factor[peaks$nsamples>14&peaks$nsamps_factor<=28] <- "15-28"  
peaks$nsamps_factor[peaks$nsamples>28&peaks$nsamps_factor<=42] <- "29-42"  
peaks$nsamps_factor[peaks$nsamples>42&peaks$nsamps_factor<=56] <- "43-56"  
peaks$nsamps_factor[peaks$nsamples>57&peaks$nsamps_factor<=70] <- "57-70"  
peaks$nsamps_factor[peaks$nsamples>71&peaks$nsamps_factor<=84] <- "71-84"  

peaks$nsamps_factor <- as.factor(peaks$nsamps_factor)

ggplot(peaks,aes(x=nsamps_factor,y=qval))+geom_boxplot() + theme_classic() + xlab("# samples")
ggplot(peaks,aes(x=nsamps_factor,y=meanscore))+geom_boxplot() +theme_classic()+ xlab("# samples") + ylab("mean -log10(q-val)")

#I will try using as cutoff for peaks found in only one sample the mean q-val among 15-28 samples
mean(peaks$meanscore[peaks$nsamps_factor=="15-28"])
mean(peaks$qval[peaks$nsamps_factor=="15-28"])
coval=mean(peaks$meanscore[peaks$nsamps_factor=="15-28"])

table(peaks$meanscore>=coval&peaks$nsamples<42)
peaks$length=peaks$end - peaks$start +1

plot(peaks$meanscore,peaks$length)
ggplot(peaks,aes(x=meanscore, y=length))+ stat_binhex()
peaks$rank <- order(peaks$meanscore)
ggplot(peaks,aes(x=rank, y=length)) + geom_point()

list_peaks_to_keep <- list(peaks_in_half_of_samples=peaks$peakID[peaks$nsamples>=42],
                           peaks_in_half_of_samples_and_lowqval=peaks$peakID[peaks$nsamples>=42|peaks$meanscore>=coval])

length(list_peaks_to_keep$peaks_in_half_of_samples)

length(list_peaks_to_keep$peaks_in_half_of_samples_and_lowqval)

##comparisons:

#all samples:
    #LPS vs no LPS
    #cohesin vs no cohesin
#Cohesin patients:
    #LPS 1 vs LPS 2
#LPS1 samples:
    #cohesin vs no cohesin
#LPS2 samples:
    #cohesin vs no cohesin

setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples/")
dir.create("all_samples")
dir.create("Cohesin_patients")
dir.create("LPS1")
dir.create("LPS2")
setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/")
dir.create("all_samples")
dir.create("Cohesin_patients")
dir.create("LPS1")
dir.create("LPS2")

# #all samples Cohesin without filtering
# setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/")
# setwd("all_samples/")
# dir.create("cohesin_vs_nocohesin")
# setwd("cohesin_vs_nocohesin/")
# ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
#                                   colData = coldata,
#                                   design = ~ Cohesin)
# 
# keep <- rowSums(counts(ddsMat)) > 1
# ddsMat <- ddsMat[keep,]
# nrow(ddsMat)
# 
# ddsMat <- estimateSizeFactors(ddsMat)
# dds <- DESeq(ddsMat)
# resultsNames(dds)
# 
# res <- results(dds,alpha = 0.05)
# write(capture.output(resultsNames(dds)),"summary_res.txt")
# write(capture.output(summary(res)),"summary_res.txt",append = T)
# res <- as.data.frame(res)
# write.csv(res[order(res$padj),],"all_samples_cohesin_vs_nocohesin_results.csv")



for (filt in names(list_peaks_to_keep)) {
  setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
  setwd(filt)
  peaks_to_keep=list_peaks_to_keep[[filt]]
  
  #filter countdata
  Tcountdata=countdata[rownames(countdata)%in%peaks_to_keep,]
  
  
  #all samples
  #LPS vs no LPS
  setwd("all_samples/")
  dir.create("LPS_vs_noLPS")
  setwd("LPS_vs_noLPS/")
  TddsMat <- DESeqDataSetFromMatrix(countData = Tcountdata,
                                   colData = coldata,
                                   design = ~ LPS)

  keep <- rowSums(counts(TddsMat)) > 1
  TddsMat <- TddsMat[keep,]
  nrow(TddsMat)
  
  TddsMat <- estimateSizeFactors(TddsMat)
  Tdds <- DESeq(TddsMat)
  resultsNames(Tdds)
  
  res <- results(Tdds,alpha = 0.05)
  write(capture.output(resultsNames(Tdds)),"summary_res.txt")
  write(capture.output(summary(res)),"summary_res.txt",append = T)
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],"all_samples_LPS_vs_noLPS_DESeq2_results.csv")
  
  #cohesin vs no cohesin
  setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
  setwd(filt)
  setwd("all_samples/")
  dir.create("cohesin_vs_nocohesin")
  setwd("cohesin_vs_nocohesin/")
  
  TddsMat <- DESeqDataSetFromMatrix(countData = Tcountdata,
                                    colData = coldata,
                                    design = ~ Cohesin)
  
  keep <- rowSums(counts(TddsMat)) > 1
  TddsMat <- TddsMat[keep,]
  nrow(TddsMat)
  
  TddsMat <- estimateSizeFactors(TddsMat)
  Tdds <- DESeq(TddsMat)
  resultsNames(Tdds)
  
  res <- results(Tdds,alpha = 0.05)
  write(capture.output(resultsNames(Tdds)),"summary_res.txt")
  write(capture.output(summary(res)),"summary_res.txt",append = T)
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],"all_samples_cohesin_vs_nocohesin_DESeq2_results.csv")
  
  
  #Cohesin patients:
  #LPS 1 vs LPS 2
  setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
  setwd(filt)
  setwd("Cohesin_patients/")
  dir.create("LPS_vs_noLPS")
  setwd("LPS_vs_noLPS/")
  
  #filter countdata and coldata:
  Tcoldata=coldata[coldata$Cohesin==TRUE,]
  TTcountdata <- Tcountdata[,colnames(countdata)%in%Tcoldata$samples]
  
  TddsMat <- DESeqDataSetFromMatrix(countData = TTcountdata,
                                    colData = Tcoldata,
                                    design = ~ LPS)
  
  keep <- rowSums(counts(TddsMat)) > 1
  TddsMat <- TddsMat[keep,]
  nrow(TddsMat)
  
  TddsMat <- estimateSizeFactors(TddsMat)
  Tdds <- DESeq(TddsMat)
  resultsNames(Tdds)
  
  res <- results(Tdds,alpha = 0.05)
  write(capture.output(resultsNames(Tdds)),"summary_res.txt")
  write(capture.output(summary(res)),"summary_res.txt",append = T)
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],"Cohesin_patients_LPS_vs_noLPS_DESeq2_results.csv")
  
  #LPS1 samples:
  #cohesin vs no cohesin
  setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
  setwd(filt)
  setwd("LPS1/")
  dir.create("cohesin_vs_nocohesin")
  setwd("cohesin_vs_nocohesin/")
  
  #retain only LPS1 samples:
  Tcoldata <- coldata[coldata$LPS==1,]
  TTcountdata <- Tcountdata[,colnames(countdata)%in%Tcoldata$samples]
  TddsMat <- DESeqDataSetFromMatrix(countData = TTcountdata,
                                    colData = Tcoldata,
                                    design = ~ Cohesin)
  
  keep <- rowSums(counts(TddsMat)) > 1
  TddsMat <- TddsMat[keep,]
  nrow(TddsMat)
  
  TddsMat <- estimateSizeFactors(TddsMat)
  Tdds <- DESeq(TddsMat)
  resultsNames(Tdds)
  
  res <- results(Tdds,alpha = 0.05)
  write(capture.output(resultsNames(Tdds)),"summary_res.txt")
  write(capture.output(summary(res)),"summary_res.txt",append = T)
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],"noLPS_samples_cohesin_vs_noCohesin_DESeq2_results.csv")
  
  
  #LPS2 samples:
  #cohesin vs no cohesin
  setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/")
  setwd(filt)
  setwd("LPS2/")
  dir.create("cohesin_vs_nocohesin")
  setwd("cohesin_vs_nocohesin/")
  
  Tcoldata <- coldata[coldata$LPS==2,]
  TTcountdata <- Tcountdata[,colnames(countdata)%in%Tcoldata$samples]
  TddsMat <- DESeqDataSetFromMatrix(countData = TTcountdata,
                                    colData = Tcoldata,
                                    design = ~ Cohesin)
  
  keep <- rowSums(counts(TddsMat)) > 1
  TddsMat <- TddsMat[keep,]
  nrow(TddsMat)
  
  TddsMat <- estimateSizeFactors(TddsMat)
  Tdds <- DESeq(TddsMat)
  resultsNames(Tdds)
  
  res <- results(Tdds,alpha = 0.05)
  write(capture.output(resultsNames(Tdds)),"summary_res.txt")
  write(capture.output(summary(res)),"summary_res.txt",append = T)
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],"LPS_samples_cohesin_vs_noCohesin_DESeq2_results.csv")
  
  
  }

