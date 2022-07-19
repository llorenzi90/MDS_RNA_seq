## ---------------------------
##
##
## Purpose of script: PCA of MDS RNA-seq
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-10-19
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

## ---------------------------
## ---------------------------                                                                                                            
#setwd("~/MDS/RNA-seq/analyses/")  
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


ddsMat <- DESeqDataSetFromMatrix(countData = countdata,                                                                                   
                                 colData = coldata,                                                                                       
                                 design = ~0+ patient + LPS )                                                                  

#check this warning 
# the design formula contains one or more numeric variables with integer values,
# specifying a model with increasing fold change for higher values.
# did you mean for this to be a factor? if so, first convert
# this variable to a factor using the factor() function
# the design formula contains one or more numeric variables that have mean or
# standard deviation larger than 5 (an arbitrary threshold to trigger this message).
# Including numeric variables with large mean can induce collinearity with the intercept.
# Users should center and scale numeric variables in the design to improve GLM convergence.

#in fact patient is numeric but should be a factor, same thing for LPS
coldata$patient <- as.factor(coldata$patient)
coldata$LPS <- as.factor(coldata$LPS)

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,                                                                                   
                                 colData = coldata,                                                                                       
                                 design = ~0+ patient + LPS )   

#filter non-informative rows                                                                                                              
#removing rows of the DESeqDataSet that have no counts, or only a single count across all samples. Additional weighting/filtering to impro
#ve power is applied at a later step in the workflow.                                                                                      
nrow(ddsMat)                                                                                                                              
## [1] 58294                                                                                                                              
keep <- rowSums(counts(ddsMat)) > 1                                                                                                       
ddsMat <- ddsMat[keep,]                                                                                                                   
nrow(ddsMat)                                                                                                                              
## [1] 31604                                                                                                                              
ddsMat <- estimateSizeFactors(ddsMat)                                                                                                     

#Data visualization                                                                                                                       
#transform with vst                                                                                                                       
vsd<- vst(ddsMat, blind = TRUE)                                                                                                          
head(assay(vsd), 3)                                                                                                                       

#blind = FALSE, means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected 
#variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global a
#mount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).             

#compare log2 transformation with vst                                                                                                     
library("dplyr")                                                                                                                          
library("ggplot2")                                                                                                                        

df <- bind_rows(                                                                                                                          
  as_data_frame(log2(counts(ddsMat, normalized=TRUE)[, 1:2]+1)) %>%                                                                       
    mutate(transformation = "log2(x + 1)"),                                                                                               
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))                                                                    


lvls <- c("log2(x + 1)", "vst")                                                                                                           
df$transformation <- factor(df$transformation, levels=lvls)                                                                               
colnames(df)[1:2] <- c("sample1","sample2")
ggplot(df, aes(x = sample1, y = sample2)) + geom_hex(bins = 80) +                                                                     
  coord_fixed() + facet_grid( . ~ transformation)                                                                                         


#calculate sample-sample distances                                                                                                        
sampleDists <- dist(t(assay(vsd)))                                                                                                        
sampleDists                                                                                                                               

library("pheatmap")                                                                                                                       
library("RColorBrewer")                                                                                                                   
sampleDistMatrix <- as.matrix( sampleDists )                                                                                              
colnames(sampleDistMatrix) <- NULL                                                                                                        
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)                                                                            
pm_euc <- pheatmap(sampleDistMatrix,                                                                                                                
         clustering_distance_rows = sampleDists,                                                                                          
         clustering_distance_cols = sampleDists,                                                                                          
         col = colors)                                                                                                                    
pdf("plots/Feb2022_MDS_RNA_seq_pheatmap.Eucledian.pdf")
print(pm_euc)
dev.off()

#an alternative way to calculate distance:                                                                                                
library("PoiClaClu")                                                                                                                      
poisd <- PoissonDistance(t(counts(ddsMat)))                                                                                               

samplePoisDistMatrix <- as.matrix( poisd$dd )                                                                                             
rownames(samplePoisDistMatrix) <- paste( ddsMat$name)                                                                                     
colnames(samplePoisDistMatrix) <- NULL                                                                                                    
pm_poisson <- pheatmap(samplePoisDistMatrix,                                                                                                            
         clustering_distance_rows = poisd$dd,                                                                                             
         clustering_distance_cols = poisd$dd,                                                                                             
         col = colors)                                                                                                                    

pdf("plots/Feb2022_MDS_RNA_seq_pheatmap.Poisson.pdf")
print(pm_poisson)
dev.off()

#PCA  

pca_both <- plotPCA(vsd,intgroup=c("patient","LPS"))+scale_color_discrete(name = "patient.LPS")                                                                                              
print(pca_both)
pca_patient <- plotPCA(vsd,intgroup=c("patient"))+scale_color_discrete(name = "patient")                                                                                                   
print(pca_patient)
pca_LPS <- plotPCA(vsd,intgroup=c("LPS"))+scale_color_discrete(name = "LPS")                                                                                              
print(pca_LPS)
pca_Cohesin <- plotPCA(vsd,intgroup=c("Cohesin"))+scale_color_discrete(name = "Cohesin")                                                                                              
print(pca_Cohesin)

list_to_plot <- list(patient.LPS=pca_both,Cohesin=pca_Cohesin,
                     LPS=pca_LPS,patient=pca_patient)


for (nam in names(list_to_plot)) {
  tpca=list_to_plot[[nam]]
  pdf(paste0("plots/Feb2022_MDS_RNA_seq_PCA_colour.",nam,".pdf"))
  print(tpca + theme_classic())
  dev.off()
}
