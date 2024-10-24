## ---------------------------
##
##
## Purpose of script: MDS DESeq2
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-13
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
require(DESeq2)
library("dplyr")
library("ggplot2")
library(RColorBrewer)
library(pheatmap)
## ---------------------------
countdata=read.table("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/merged_macs2peaksNoDups_MDS_all_samples.counts")

samples=colnames(countdata)
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1])
coldata=data.frame(samples,patient,LPS)

rownames(coldata) <- coldata$samples
table(rownames(coldata)==colnames(countdata))
coldata$LPS <- as.factor(coldata$LPS)

#Add data on mutations

setwd("/Users/llorenzi/Nextcloud/MDS/OncoPrint/Analisis_Pame/")

#note that all data is in sheet 3
mut_data=readxl::read_xlsx("Resumen_analisis paneles Pame.xlsx",sheet = 3)
extra_samples=readxl::read_xlsx("Resumen_analisis paneles Pame.xlsx",sheet =4)
extra_samples=extra_samples[c(1:3,8),]
MDS3=extra_samples[1:3,]

#for now I am only interested in genes mutated
MDS3[,!tolower(colnames(MDS3))%in%tolower(colnames(mut_data))] <- NA
matching_cols=match(tolower(colnames(mut_data)),tolower(colnames(MDS3)))
restof_cols=(1:ncol(MDS3))[!(1:ncol(MDS3))%in%matching_cols]
matching_cols[is.na(matching_cols)] <- restof_cols
MDS3 <- MDS3[,matching_cols]
colnames(MDS3) <- colnames(mut_data)



MDS35=extra_samples[7:8,]
colnames(MDS35)=MDS35[1,]
MDS35[,!tolower(colnames(MDS35))%in%tolower(colnames(mut_data))] <- NA
matching_cols=match(tolower(colnames(mut_data)),tolower(colnames(MDS35)))
restof_cols=(1:ncol(MDS35))[!(1:ncol(MDS35))%in%matching_cols]
matching_cols[is.na(matching_cols)] <- restof_cols
MDS35 <- MDS35[,matching_cols]
colnames(MDS35) <- colnames(mut_data)

mut_data=rbind(mut_data,MDS3,MDS35[-1,])




length(unique(mut_data$Sample))
length(unique(mut_data$Gene))
#I just need the genes mutated per patient
mat=as.matrix(unclass(table(mut_data$Gene,mut_data$Sample)))
#Classify genes by pathwaty using the Leukemia 2014 paper as reference
S3table <- read.csv("../OncoPrint/Suppl_TableS3_Leukemia_2014.csv")
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

mut_data$Pathway=Pathway_annot$Pathway[match(mut_data$Gene,Pathway_annot$Gene)]
table(mut_data$Pathway)
mut_data$Sample[mut_data$Pathway=="Cohesin"]
mut_data$Gene[mut_data$Pathway=="Cohesin"]

coldata$Cohesin=coldata$patient%in%mut_data$Sample[mut_data$Pathway=="Cohesin"]
table(coldata$Cohesin)

#Separate LPS vs no LPS



for (lps in c(1,2)) {
  
  setwd("~/Nextcloud/MDS/ATAC-seq/analyses/")
  odir=paste0("PCA_and_clustering/","LPS_",lps)
  dir.create(odir)
  setwd(odir)
  sink("PCA_cluster.log")
  
  tcoldata=coldata[coldata$LPS==lps,]
  tcdata=countdata[,colnames(countdata)%in%tcoldata$samples]

  
  ddsMat <- DESeqDataSetFromMatrix(countData = tcdata,
                                   colData = tcoldata,
                                   design = ~ Cohesin)

  keep <- rowSums(counts(ddsMat)) > 1
  ddsMat <- ddsMat[keep,]
  nrow(ddsMat)
  
  ddsMat <- estimateSizeFactors(ddsMat)
  vsd<- vst(ddsMat)
  
  df <- bind_rows(
    as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  
  lvls <- c("log2(x + 1)", "vst")
  df$transformation <- factor(df$transformation, levels=lvls)
  colnames(df)[1:2] <- c("x","y")
  
  png("vst_normalization_effect_vs_log2.png")
  g=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  print(g)
  dev.off()
  
  
  sampleDists <- dist(t(assay(vsd)))
  sampleDists
  
  
  sampleDistMatrix <- as.matrix( sampleDists )
  colnames(sampleDistMatrix) <- NULL
  
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
  graphics.off()
  

  #PCA
  pdf("PCA_Cohesin.pdf")
  g=plotPCA(vsd,intgroup=c("Cohesin"),)
  g$labels$colour="Cohesin"
  print(g)
  dev.off()
  
  sink()
  
  ###DESeq2
  setwd("~/Nextcloud/MDS/ATAC-seq/analyses/")
  odir=paste0("DESeq2/","LPS_",lps)
  dir.create(odir)
  setwd(odir)
  
  sink("DESeq2.log")
  dds <- DESeq(ddsMat)
  resultsNames(dds)
  
  res <- results(dds,alpha = 0.05)
  write(capture.output(resultsNames(dds)),"summary_res.txt")
  write(capture.output(summary(res)),"summary_res.txt",append = T)

  
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],paste0("LPS_",lps,"_DESeq2_results.csv"))
  sink()
}



#### select only Cohesin mutants
odir="Cohesin_mutants"
#do it with
# 1)all peaks 
# 2) only peaks that are present in Cohesin mutants
# 3) only peaks that have enough counts in Cohesin mutants
setwd("~/Nextcloud/MDS/ATAC-seq/analyses/PCA_and_clustering")
dir.create(odir)
setwd(odir)


##################### 1) all peaks ####
tcoldata=coldata[coldata$Cohesin==TRUE,]
tcdata=countdata[,colnames(countdata)%in%tcoldata$samples]

ddsMat <- DESeqDataSetFromMatrix(countData = tcdata,
                                 colData = tcoldata,
                                 design = ~ LPS)

keep <- rowSums(counts(ddsMat)) > 1
table(keep)
ddsMat <- ddsMat[keep,]
nrow(ddsMat)

ddsMat <- estimateSizeFactors(ddsMat)

##PCA and clustering
setwd("~/Nextcloud/MDS/ATAC-seq/analyses/PCA_and_clustering")
dir.create(odir)
setwd(odir)
dir.create("all_peaks")
setwd("all_peaks/")

sink("PCA_cluster.log")
vsd<- vst(ddsMat)

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))


lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
colnames(df)[1:2] <- c("x","y")


png("vst_normalization_effect_vs_log2.png")
g=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
print(g)
dev.off()


sampleDists <- dist(t(assay(vsd)))
sampleDists


sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df=data.frame(LPS=tcoldata$LPS)
rownames(df)=tcoldata$samples
table(rownames(df)==rownames(sampleDistMatrix))

df$LPS=as.factor(df$LPS)

pdf("pheatmap_sample_clustering.pdf")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6,
           annotation_row = df)
print(p)
graphics.off()


#PCA
pdf("PCA_LPS.pdf")
g=plotPCA(vsd,intgroup=c("LPS"),)
g$labels$colour="LPS"
print(g)
dev.off()

sink()

##### DESeq2
setwd("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/")

dir.create(odir)
setwd(odir)
dir.create("all_peaks")
setwd("all_peaks/")

sink("DESeq2.log")
dds <- DESeq(ddsMat)
resultsNames(dds)

res <- results(dds,alpha = 0.05)

write(capture.output(resultsNames(dds)),"summary_res.txt")
write(capture.output(summary(res)),"summary_res.txt",append = T)


res <- as.data.frame(res)
write.csv(res[order(res$padj),],paste0("LPS_",lps,"_DESeq2_results.csv"))
sink()


############## 2) only peaks that are present in Cohesin mutants ####
peaks=read.table("~/Nextcloud/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames",header = T)
peaks$peakID=paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)
Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]

samps_as_vect=sapply(peaks$samples,function(x)strsplit(x,split = ","))
samps_as_vect <- sapply(samps_as_vect, function(x){t=gsub("\\.[1-2]","",x)
t=gsub("AT-","",t)})
peak_in_Cohesin_patient=sapply(samps_as_vect, function(x) any(x%in%Cohesin_patients))
class(peak_in_Cohesin_patient)
table(peak_in_Cohesin_patient)

Cohesin_peaks=peaks$peakID[peak_in_Cohesin_patient]

tcoldata=coldata[coldata$Cohesin==TRUE,]
tcdata=countdata[rownames(countdata)%in%Cohesin_peaks,colnames(countdata)%in%tcoldata$samples]

ddsMat <- DESeqDataSetFromMatrix(countData = tcdata,
                                 colData = tcoldata,
                                 design = ~ LPS)

keep <- rowSums(counts(ddsMat)) > 1
table(keep)
ddsMat <- ddsMat[keep,]
nrow(ddsMat)

ddsMat <- estimateSizeFactors(ddsMat)

##PCA and clustering
setwd("~/Nextcloud/MDS/ATAC-seq/analyses/PCA_and_clustering")
dir.create(odir)
setwd(odir)
dir.create("only_Cohesin_peaks")
setwd("only_Cohesin_peaks/")

sink("PCA_cluster.log")
vsd<- vst(ddsMat)

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))


lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
colnames(df)[1:2] <- c("x","y")


png("vst_normalization_effect_vs_log2.png")
g=ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
print(g)
dev.off()


sampleDists <- dist(t(assay(vsd)))
sampleDists


sampleDistMatrix <- as.matrix( sampleDists )
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df=data.frame(LPS=tcoldata$LPS)
rownames(df)=tcoldata$samples
table(rownames(df)==rownames(sampleDistMatrix))

df$LPS=as.factor(df$LPS)

pdf("pheatmap_sample_clustering.pdf")
p=pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           fontsize_row = 6,
           annotation_row = df)
print(p)
graphics.off()


#PCA
pdf("PCA_LPS.pdf")
g=plotPCA(vsd,intgroup=c("LPS"),)
g$labels$colour="LPS"
print(g)
dev.off()

sink()

##### DESeq2
setwd("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/")

dir.create(odir)
setwd(odir)
dir.create("only_Cohesin_peaks")
setwd("only_Cohesin_peaks/")

sink("DESeq2.log")
dds <- DESeq(ddsMat)
resultsNames(dds)

res <- results(dds,alpha = 0.05)

write(capture.output(resultsNames(dds)),"summary_res.txt")
write(capture.output(summary(res)),"summary_res.txt",append = T)


res <- as.data.frame(res)
write.csv(res[order(res$padj),],paste0("LPS_",lps,"_DESeq2_results.csv"))
sink()
