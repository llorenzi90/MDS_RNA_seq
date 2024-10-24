## ---------------------------
##
##
## Purpose of script: PCA and clustering MDS samples
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-12
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
require(tidyverse)
library(edgeR)
library(DESeq2)
## ---------------------------

####PCA and clustering
countdata=read.table("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/merged_macs2peaksNoDups_MDS_all_samples.counts")

rownames(countdata) <- countdata$peakID
samples=colnames(countdata)
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1])
coldata=data.frame(samples,patient,LPS)

rownames(coldata) <- coldata$samples
table(rownames(coldata)==colnames(countdata))
coldata$LPS <- as.factor(coldata$LPS)

#add batch
samples=read.table("~/Nextcloud/MDS/ATAC-seq/analyses/list_of_samples.txt")
samples$sampleID <- gsub("AT-","",samples$V1)
table(samples$sampleID%in%coldata$samples)

samples$batch="2"
samples$batch[grep("AT-",samples$V1)]="1"

coldata$batch=samples$batch[match(coldata$samples,samples$sampleID)]

design <- model.matrix(~ coldata$LPS+ coldata$patient )
colnames(design) <- gsub("coldata\\$","",colnames(design))
design

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ LPS+ patient)

ddsMat <- estimateSizeFactors(ddsMat)
#keep <- rowSums(counts(ddsMat)) > 100

vsd<- vst(ddsMat)


pdf("~/Nextcloud/MDS/ATAC-seq/analyses/PCA_and_clustering/All_samples/PCA_batch.pdf")
g=plotPCA(vsd,intgroup=c("batch"))
g$labels$colour="batch"
g
dev.off()
