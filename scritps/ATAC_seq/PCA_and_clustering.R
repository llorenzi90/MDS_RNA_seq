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

design <- model.matrix(~ coldata$LPS+ coldata$patient )
colnames(design) <- gsub("coldata\\$","",colnames(design))
design

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ LPS+ patient)

ddsMat <- estimateSizeFactors(ddsMat)
#keep <- rowSums(counts(ddsMat)) > 100

vsd<- vst(ddsMat)

library("dplyr")
library("ggplot2")

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))


lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
colnames(df)[1:2] <- c("x","y")
#png(paste0(odir,"/vst_normalization_effect_vs_log2.png"))
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

#dev.off()


sampleDists <- dist(t(assay(vsd)))
sampleDists


sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
setwd("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/PCA_and_clustering")
pdf("pheatmap_sample_clustering.pdf")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6)
print(p)
dev.off()

#PCA
pdf("PCA_patient_LPS.pdf")
g=plotPCA(vsd,intgroup=c("patient","LPS"))
g$labels$colour="patient_LPS"
g
dev.off()

pdf("PCA_LPS.pdf")
g=plotPCA(vsd,intgroup=c("LPS"),)
g$labels$colour="LPS"
g
dev.off()

pdf("PCA_patient.pdf")

g=plotPCA(vsd,intgroup=c("patient"))
g$labels$colour="patient"
g
dev.off()




#Add data on mutations

setwd("/Users/llorenzi/Nextcloud/MDS/OncoPrint/OncoPrint")

#note that all data is in sheet 3
raw_data=readxl::read_xlsx("Resumen_analisis paneles Pame.xlsx",sheet = 3)

length(unique(raw_data$Sample))
length(unique(raw_data$Gene))
#for the type of plot I want to make, similar to the one in doi:10.1038/leu.2013.336, I don't need to use the type of mutation for now
#I just need the genes mutated per patient
mat=as.matrix(unclass(table(raw_data$Gene,raw_data$Sample)))
S3table <- read.csv("Suppl_TableS3_Leukemia_2014.csv")
table(rownames(mat)%in%S3table$Gene)
rownames(mat)[!rownames(mat)%in%S3table$Gene]
#there are 5 genes that are not in the paper supp table
# but U2AF1/1L5 are just U2AF1 in the paper
Pathway_annot <- data.frame(Gene=rownames(mat),Pathway=S3table$Pathway[match(rownames(mat),S3table$Gene)])
table(S3table$Pathway)
# add the missing annotations manually
#CSF3R is a receptor 
manual_annot <- c(CSF3R="Receptors/Kinases")
#CUX1 is a transcription factor "Transcription" group
manual_annot <- c(manual_annot,CUX1="Transcription")
#KMT2A is a methyl-tranferase
manual_annot <- c(manual_annot,KMT2A="Chromatin modification")
# STBP1 according to Wikipedia:
#The SETBP1 gene provides instructions for making a protein known as the SET binding protein 1, which is widely distributed throughout somatic cells. The protein is known to bind to another protein called SET. SETBP1 is a DNA-binding protein that forms part of a group of proteins that act together on histone methylation to make chromatin more accessible and regulate gene expression.[6] There is still more to learn about the overall function of the SETBP1 protein and the effect of SET binding.
#So it is a chromatin modifier
manual_annot <- c(manual_annot,SETBP1="Chromatin modification")
manual_annot <- c(manual_annot,`U2AF1;U2AF1L5`="RNA splicing")

Pathway_annot$Pathway[Pathway_annot$Gene%in%names(manual_annot)] <- manual_annot[match(Pathway_annot$Gene[Pathway_annot$Gene%in%names(manual_annot)],
                                                                                       names(manual_annot))]

grp=as.factor(Pathway_annot$Pathway)

raw_data$Pathway=Pathway_annot$Pathway[match(raw_data$Gene,Pathway_annot$Gene)]
table(raw_data$Pathway)
raw_data$Sample[raw_data$Pathway=="Cohesin"]
raw_data$Gene[raw_data$Pathway=="Cohesin"]

coldata$Cohesin=coldata$patient%in%raw_data$Sample[raw_data$Pathway=="Cohesin"]
table(coldata$Cohesin)

setwd("~/Nextcloud/MDS/ATAC-seq/analyses/PCA_and_clustering/")
dir.create("With_Cohesin")
setwd("With_Cohesin/")

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ LPS+ patient)

ddsMat <- estimateSizeFactors(ddsMat)
#keep <- rowSums(counts(ddsMat)) > 100

vsd<- vst(ddsMat)

library("dplyr")
library("ggplot2")

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))


lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
colnames(df)[1:2] <- c("x","y")
#png(paste0(odir,"/vst_normalization_effect_vs_log2.png"))
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

#dev.off()


sampleDists <- dist(t(assay(vsd)))
sampleDists


sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df=data.frame(Cohesin=coldata$Cohesin)
rownames(df)=coldata$samples
table(rownames(df)==rownames(sampleDistMatrix))

df$Cohesin=as.factor(df$Cohesin)
pdf("pheatmap_sample_clustering.pdf")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6,
           annotation_row = df)
print(p)
dev.off()

#PCA
pdf("PCA_Cohesin.pdf")
g=plotPCA(vsd,intgroup=c("Cohesin"),)
g$labels$colour="Cohesin"
g
dev.off()


pdf("PCA_Cohesin_LPS.pdf")
g=plotPCA(vsd,intgroup=c("Cohesin","LPS"))
g$labels$colour="Cohesin_LPS"
g
dev.off()



pdf("PCA_patient_Cohesin.pdf")

g=plotPCA(vsd,intgroup=c("Cohesin","patient"))
g$labels$colour="Cohesin:patient"
g
dev.off()

