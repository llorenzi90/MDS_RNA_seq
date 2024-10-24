## ---------------------------
##
##
## Purpose of script:Plot mean coverage as heatmap
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-09
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:Similar to what I did in plotmeanCoverage_aroundTSSs.R
##  But now I want to take for each gene the average across samples
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
library(rtracklayer)
## ---------------------------
#load TSS coordinates
tss.inf_genes=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_inflammatory_genes_hg38.txt",header = T)

#generate granges object with -3000 - +3000 bps around the 200 TSSs
extended_granges <- GRanges(tss.inf_genes$seqid,
                            IRanges(start = tss.inf_genes$start - 2999,
                                    width = 6000))


#for each sample, import the regions of the bigwig files matching the regions around TSSs 
datadir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"

bwfiles=list.files(datadir,full.names = T,pattern = "bam.CPM.bw")

mut_data=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/mut_data_all_patients.txt",sep = "\t",header = T)

samples=gsub(".sorted.markedDups.primary_chr.proper_pairs.minq2.bam.CPM.bw","",basename(bwfiles))
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=gsub("AT-","",sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1]))
Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]
Cohesin=patient%in%Cohesin_patients

list_meanCov=list()
splitted_ranges=list()
list_all_TSSposCPM <- list()
i=0
for (bw in bwfiles) {
  i=i+1
  tmp_tssCPM=import(bw,which=extended_granges)
  #extend the GRanges object so each position around -3000 - +3000 bp of each TSS has a CPM value
  all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
  list_all_TSSposCPM[[i]] <- all_TSSposCPM
}
#this list has one element for each sample
#each element has the CPM value for each position across the 200 genes

#take average from this
all_TSSpos_meanCPM <- apply(as.data.frame(list_all_TSSposCPM),1,mean)
Cohesin_TSSpos_meanCPM <- apply(as.data.frame(list_all_TSSposCPM[Cohesin]),1,mean)
NONCohesin_TSSpos_meanCPM <- apply(as.data.frame(list_all_TSSposCPM[!Cohesin]),1,mean)

Cohesin_meanCPM_all_genes <- as.data.frame(split(Cohesin_TSSpos_meanCPM, ceiling(seq_along(all_TSSposCPM)/6000)))
NONCohesin_meanCPM_all_genes <- as.data.frame(split(NONCohesin_TSSpos_meanCPM, ceiling(seq_along(all_TSSposCPM)/6000)))

data_to_plot=list(Non_Cohesin_patients=NONCohesin_meanCPM_all_genes,
                  Cohesin_patients=Cohesin_meanCPM_all_genes)

#sort data to plot by total CPM
total_CPMs <- lapply(data_to_plot, function(x)apply(x, 2,sum))

mean_totalCPM <- apply(as.data.frame(total_CPMs),1,mean)

data_to_plot <- lapply(data_to_plot, function(x){
  return(t(x[,order(mean_totalCPM,decreasing = T)]))
  
})
ncol(data_to_plot$Non_Cohesin_patients)
data_to_plot <- lapply(data_to_plot,function(x) {colnames(x) = -2999:3000
return(x)})
colnames(data_to_plot$Cohesin_patients)

library(pheatmap)
getwd()
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/inflammatory_genes/")
dir.create("heatmaps")
setwd("heatmaps/")
pdf("test.pdf")
pheatmap(data_to_plot$Non_Cohesin_patients,cluster_rows = F,cluster_cols = F,show_colnames = F,show_rownames = F)
dev.off()

pdf("test.Cohesin_patients.pdf")
pheatmap(data_to_plot$Cohesin_patients,cluster_rows = F,cluster_cols = F,show_colnames = F,show_rownames = F)
dev.off()


pdf("test.colnames.pdf")
pheatmap(data_to_plot$Non_Cohesin_patients,cluster_rows = F,cluster_cols = F,show_colnames = T,show_rownames = F)
dev.off()

#require(grid)
labcol <- rep("",6000)
labcol[1] <- "-3kb"
labcol[1001] <- "-2kb"
labcol[2001] <- "-1kb"
labcol[3000] <- "TSS"
labcol[4000] <- "1kb"
labcol[5000] <- "2kb"
labcol[6000] <- "3kb"

heat <- pheatmap(data_to_plot$Non_Cohesin_patients,cluster_rows = F,cluster_cols = F,
                 show_colnames = T,labels_col = labcol, show_rownames = F)

dev.off()
pdf("test.colnames.pdf")
heat
dev.off()

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")
hmcols<-colorRampPalette(c("#0072B2","white","#D55E00"))(256)

heat <- ComplexHeatmap::pheatmap(data_to_plot$Non_Cohesin_patients,cluster_rows = F,cluster_cols = F,
                         show_colnames = T,labels_col = labcol, show_rownames = F,
                         color = hmcols, 
                         use_raster=F,
                         fontsize_number = 9, column_names_side = c("top"),
                         angle_col = c("0"))

heat2 <- ComplexHeatmap::pheatmap(data_to_plot$Cohesin_patients,cluster_rows = F,cluster_cols = F,
                                          show_colnames = T,labels_col = labcol, show_rownames = F,
                                          color = hmcols, 
                                          use_raster=F,
                                          fontsize_number = 9, column_names_side = c("top"),
                                          angle_col = c("0"))

pdf("test.colnames.pdf")
heat
dev.off()

graphics.off()
pdf("test.mfrow.pdf")
par(mfrow=c(1,2))
heat
heat2
dev.off()

pdf("test.concatenate.pdf")
heat + heat2
dev.off()

#chequear si es log scale 

heat <- ComplexHeatmap::pheatmap(log2(data_to_plot$Non_Cohesin_patients+0.01),cluster_rows = F,cluster_cols = F,
                                 show_colnames = T,labels_col = labcol, show_rownames = F,
                                 color = hmcols, 
                                 use_raster=F,
                                 fontsize_number = 9, column_names_side = c("top"),
                                 angle_col = c("0"))
heat
#plot entre -1.5 to 1.5
#check color scale between plots
#usar escala con cero como minimo valor 

heat <- ComplexHeatmap::pheatmap(log2(data_to_plot$Non_Cohesin_patients+0.01),cluster_rows = F,cluster_cols = F,
                                 show_colnames = T,labels_col = labcol, show_rownames = F,
                                  
                                 use_raster=F,
                                 fontsize_number = 9, column_names_side = c("top"),
                                 angle_col = c("0"))

heat

#hacer con atacseq peaks diferentes entre muestras de cohesina y el resto


cbPalette <- c("#999999", 
               "#E69F00", 
               "#56B4E9", 
               "#009E73", 
               "#F0E442", 
               "#0072B2", 
               "#D55E00", 
               "#CC79A7") #source http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

hmcols<-colorRampPalette(c("white",cbPalette[6],cbPalette[7]))(256)
ComplexHeatmap::pheatmap(log2(data_to_plot$Non_Cohesin_patients+0.01),cluster_rows = F,cluster_cols = F,
                         show_colnames = T,labels_col = labcol, show_rownames = F,
                         color = hmcols,
                         use_raster=F,
                         fontsize_number = 9, column_names_side = c("top"),
                         angle_col = c("0"))

ComplexHeatmap::pheatmap(data_to_plot$Non_Cohesin_patients,cluster_rows = F,cluster_cols = F,
                         show_colnames = T,labels_col = labcol, show_rownames = F,
                         color = hmcols,
                         use_raster=F,
                         fontsize_number = 9, column_names_side = c("top"),
                         angle_col = c("0"))
boxplot(log2(as.vector(as.matrix(data_to_plot$Non_Cohesin_patients))+0.01))
plot(density(log2(as.vector(as.matrix(data_to_plot$Non_Cohesin_patients))+0.01)))

heat2 <- ComplexHeatmap::pheatmap(data_to_plot$Cohesin_patients,cluster_rows = F,cluster_cols = F,
                                  show_colnames = T,labels_col = labcol, show_rownames = F,
                                  color = hmcols, 
                                  use_raster=F,
                                  fontsize_number = 9, column_names_side = c("top"),
                                  angle_col = c("0"),
                                    legend=F)

heat + heat2
