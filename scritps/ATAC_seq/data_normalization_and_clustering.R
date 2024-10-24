## ---------------------------
##
##
## Purpose of script: normalization and sample clustering MDS samples
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-07-27
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
library(csaw)
library(edgeR)
library(DESeq2)
## ---------------------------
wdir=("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/results/csaw/")
setwd(wdir)
#read counts on sliding windows 
window_data <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/win.data.MDS.RDS")
#read background filtering stats based on 10Kb bins
filter.stat <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/windows.filter.stat.globalBackground.RDS")

#add info on samples
bam.paths=window_data$bam.files
samples=gsub("(/scratch/llorenzi/MDS/ATAC-seq/AT-MDS[0-9]*.[12]/AT-)(.*)(.primary_chr.quality.concordant.sorted.markedDups.bam)",
             "\\2",
             bam.paths)

LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1])
coldata=data.frame(samples,patient,LPS)

coldata$LPS <- as.factor(coldata$LPS)

design <- model.matrix(~ coldata$LPS+ coldata$patient )
colnames(design) <- gsub("coldata\\$","",colnames(design))
design

#####what cutoff to choose to filter windows??
col <- list("2"="green","3"="red","5"="blue")

png("logFC_from_global_background_histogram.png")
hist(filter.stat$filter, xlab="Log-fold change from global background", 
       breaks=100, main="", col="grey80", xlim=c(0, 5))
for (fc in c(2,3,5)) {
    abline(v=log2(fc), col=col[[as.character(fc)]], lwd=2)
 
}
dev.off()
#####Analyze sample relationships

library("pheatmap")
library("RColorBrewer")
for(fc in c(2,3)){
  setwd(wdir)
  odir <- paste0("plots/filtering_cutoff_",fc)
  dir.create(odir)
  fil <- filter.stat$filter > log2(fc)
  
  filtered.data <- window_data[fil,]
  countData <- assay(filtered.data)
  colnames(countData) <- samples
  ddsMat <- DESeqDataSetFromMatrix(countData = assay(filtered.data),
                                   colData = coldata,
                                   design = ~ LPS+ patient)
  
  nrow(ddsMat)
  ## [1] 58294
  
  ddsMat <- estimateSizeFactors(ddsMat)
  
  #Data visualization
  #transform with vst
  vsd<- vst(ddsMat, blind = FALSE)
  #blind = FALSE, means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
  
  #compare log2 transformation with vst
  library("dplyr")
  library("ggplot2")
  
  df <- bind_rows(
    as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  
  lvls <- c("log2(x + 1)", "vst")
  df$transformation <- factor(df$transformation, levels=lvls)
  colnames(df)[1:2] <- c("x","y")
  png(paste0(odir,"/vst_normalization_effect_vs_log2.png"))
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  
  dev.off()
  
  
  #calculate sample-sample distances
  sampleDists <- dist(t(assay(vsd)))
  sampleDists
  
  
  sampleDistMatrix <- as.matrix( sampleDists )
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf(paste0(odir,"/pheatmap_sample_clustering.pdf"))
  p=pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
  print(p)
  dev.off()
  
  png(paste0(odir,"/pheatmap_sample_clustering.png"))
  p=pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
  print(p)
  dev.off()
  # #an alternative way to calculate distance:
  # library("PoiClaClu")
  # poisd <- PoissonDistance(t(counts(ddsMat)))
  # 
  # samplePoisDistMatrix <- as.matrix( poisd$dd )
  # rownames(samplePoisDistMatrix) <- paste( ddsMat$name)
  # colnames(samplePoisDistMatrix) <- NULL
  # pheatmap(samplePoisDistMatrix,
  #          clustering_distance_rows = poisd$dd,
  #          clustering_distance_cols = poisd$dd,
  #          col = colors)
  
  #PCA
  png(paste0(odir,"/PCA_patient_LPS.png"))
  g=plotPCA(vsd,intgroup=c("patient","LPS"))
  g$labels$colour="patient_LPS"
  g
  dev.off()
  
  png(paste0(odir,"/PCA_LPS.png"))
  g=plotPCA(vsd,intgroup=c("LPS"),)
  g$labels$colour="LPS"
  g
  dev.off()
  
  png(paste0(odir,"/PCA_patient.png"))
  
  g=plotPCA(vsd,intgroup=c("patient"))
  g$labels$colour="patient"
  g
  dev.off()
  }
  



###with CPM
adjc <- calculateCPM(filtered.data, use.offsets=FALSE)
colnames(adjc) <- samples
df <- bind_rows(
  as_data_frame(adjc[,1:2]) %>%
    mutate(transformation = "log2.CPM"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))


lvls <- c("log2.CPM", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
colnames(df)[1:2] <- c("x","y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)


#From csaw manual
# Check if there are trended biases between samples. 
#This refers to a systematic fold-difference in per-window coverage between samples 
#that changes according to the average abundance of the window. 

#first calculate the scaled average per filtered window
win.ab <- scaledAverage(filtered.data)
#calculate CPM for each window in each sample
adjc <- calculateCPM(filtered.data, use.offsets=FALSE)

# logfc <- adjc[,1] - adjc[,2]
# smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
#               xlab="Average abundance", ylab="Log-fold change")
# 
# lfit <- smooth.spline(logfc~win.ab, df=5)
# o <- order(win.ab)
# lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

#for each pair-wise comparison make a plot
ncol(adjc)
libsizes <- colSums(assay(filtered.data))/1000000
names(libsizes) <- samples
dev.new()

par(mfrow=c(4,6))
for (i in 2:ncol(adjc)) {
  logfc <- adjc[,1] - adjc[,i]
  smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
                xlab="Average abundance", ylab="Log-fold change",main=paste0(samples[1]," vs ",samples[i]))
  
  lfit <- smooth.spline(logfc~win.ab, df=5)
  o <- order(win.ab)
  lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)
}

filtered.data <- normOffsets(filtered.data)
head(assay(filtered.data, "offset"))

rm.adjc <- calculateCPM(filtered.data, use.offsets=TRUE)
norm.fc <- norm.adjc[,4]-norm.adjc[,1]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
              xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(norm.fc~win.ab, df=5)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

####Sample clustering and normalization methods for visualization
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~0+rep + time)

df <- bind_rows(
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

