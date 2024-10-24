## ---------------------------
##
##
## Purpose of script: MDS DESeq2 RNA-seq Cohesin vs other, LPS vs no-LPS
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-18
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



## ---------------------------
require(tidyverse)
require(data.table)
require(DESeq2)
library("dplyr")
library("ggplot2")
library(RColorBrewer)
library(pheatmap)
## ---------------------------
conversion_ENSMBL_NCBI <- read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/biomart_tables/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv",
                                     sep = "\t",
                                     header = T)
conversion_ENSMBL_NCBI <- conversion_ENSMBL_NCBI[!is.na(conversion_ENSMBL_NCBI$Gene.name),]
conversion_ENSMBL_NCBI <- conversion_ENSMBL_NCBI[conversion_ENSMBL_NCBI$Gene.name!="",]


setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/")                                 

# if (!requireNamespace("BiocManager", quietly = TRUE))                                                                                   
#   install.packages("BiocManager")                                                                                                       
#                                                                                                                                         
# BiocManager::install("DESeq2")                                                                                                          

library(DESeq2)                                                                                                                           
#extrafont::loadfonts()                                                                                                                   
#countdata=read.csv("htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv") 
countdata=read.csv("htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv") 

coldata <- read.csv("../../processed_samples.metadata.csv")
colnames(countdata) <- gsub("X","",colnames(countdata)  )            
rownames(countdata) <- countdata$gene_id                                                                                                  
countdata <- countdata[,-1]                                                                                                               
coldata$patient.LPS <- paste(coldata$patient,coldata$LPS,sep = ".")
table(colnames(countdata)%in%coldata$patient.LPS)
#filter coldata: 
#coldata rows have to have the same order than columns in countdata
coldata <- coldata[match(colnames(countdata),coldata$patient.LPS),]

#filter all-zero genes
countdata <- countdata[rowSums(countdata)!=0,]


# patient and LPS are numeric but should be factors
coldata$patient <- as.factor(coldata$patient)
coldata$LPS <- as.factor(coldata$LPS)

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,                                                                                   
                                 colData = coldata,                                                                                       
                                 design = ~0+ patient + LPS )   

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

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/")
dir.create("DESeq2")
deseqdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/DESeq2/"
setwd(deseqdir)
dir.create("all_samples")
dir.create("Cohesin_patients")
dir.create("LPS1")
dir.create("LPS2")



run_DESeq_analysis <- function(ddsmat,prefix=""){
  keep <- rowSums(counts(ddsmat)) > 1
  ddsMat <- ddsmat[keep,]
  
  ddsMat <- estimateSizeFactors(ddsMat)
  dds <- DESeq(ddsMat)
  resultsNames(dds)
  res <- results(dds,alpha = 0.05)
  write(capture.output(resultsNames(dds)),paste0(prefix,"_summary_res.txt"))
  write(capture.output(summary(res)),paste0(prefix,"_summary_res.txt"),append = T)
  res <- as.data.frame(res)
  write.csv(res[order(res$padj),],paste0(prefix,"DESeq2_results.csv"))
  
  rank <- as.data.frame(res[order(res$stat,decreasing = T),])
  rank$gene <- rownames(rank)
  rank$gene[rank$gene%in%conversion_ENSMBL_NCBI$Gene.stable.ID.version] <- 
  conversion_ENSMBL_NCBI$Gene.name[match(rank$gene[rank$gene%in%conversion_ENSMBL_NCBI$Gene.stable.ID.version],
                                              conversion_ENSMBL_NCBI$Gene.stable.ID.version)]
  dir.create("rank_files")
  write.table(rank[,c("gene","stat")],paste0("rank_files/",prefix,"_DESeq2.rnk"),quote = F,col.names = F,row.names = F,sep = "\t")
  
}


 
  #all samples
  #LPS vs no LPS
setwd(deseqdir)
setwd("all_samples/")
dir.create("LPS_vs_noLPS")
setwd("LPS_vs_noLPS/")
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                    colData = coldata,
                                    design = ~ LPS)
run_DESeq_analysis(ddsmat = ddsMat,prefix = "all_samples_LPSvsNoLPS")

#cohesin vs no cohesin
setwd(deseqdir)
setwd("all_samples/")
dir.create("cohesin_vs_nocohesin")
setwd("cohesin_vs_nocohesin/")
  
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                    colData = coldata,
                                    design = ~ Cohesin)
run_DESeq_analysis(ddsmat = ddsMat,prefix = "all_samples_CohesinvsNoCohesin")
  
  
#Cohesin patients:
#LPS 1 vs LPS 2
setwd(deseqdir)
setwd("Cohesin_patients/")
dir.create("LPS_vs_noLPS")
setwd("LPS_vs_noLPS/")
  
#filter countdata and coldata:
Tcoldata=coldata[coldata$Cohesin==TRUE,]
Tcountdata <- countdata[,colnames(countdata)%in%Tcoldata$patient.LPS]
  
TddsMat <- DESeqDataSetFromMatrix(countData = Tcountdata,
                                    colData = Tcoldata,
                                    design = ~ LPS)
run_DESeq_analysis(TddsMat,"Cohesinpatients_LPSvsNoLPS") 

#LPS1 samples:
#cohesin vs no cohesin
setwd(deseqdir)
setwd("LPS1/")
dir.create("cohesin_vs_nocohesin")
setwd("cohesin_vs_nocohesin/")
  
#retain only LPS1 samples:
Tcoldata <- coldata[coldata$LPS==1,]
Tcountdata <- countdata[,colnames(countdata)%in%Tcoldata$patient.LPS]
TddsMat <- DESeqDataSetFromMatrix(countData = Tcountdata,
                                    colData = Tcoldata,
                                    design = ~ Cohesin)
  
run_DESeq_analysis(TddsMat,"LPS1_CohesinvsNoCohesin")  
  
#LPS2 samples:
#cohesin vs no cohesin
setwd(deseqdir)
setwd("LPS2/")
dir.create("cohesin_vs_nocohesin")
setwd("cohesin_vs_nocohesin/")
  
Tcoldata <- coldata[coldata$LPS==2,]
Tcountdata <- countdata[,colnames(countdata)%in%Tcoldata$patient.LPS]
TddsMat <- DESeqDataSetFromMatrix(countData = Tcountdata,
                                    colData = Tcoldata,
                                    design = ~ Cohesin)
  
  
run_DESeq_analysis(TddsMat,"LPS2_CohesinvsNoCohesin")  


