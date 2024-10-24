#I'll try PCA with 3 types of normalization: 
#A)CPM
#B)DESeq2’s median of ratios method (see https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html)
#C)variance stabilizing transformation (also DESeq2)

#A)
libsizes <- colSums(countdata)
scaled_counts_CPM <-t( apply(counts(ddsMat,normalized=F), 1, function(x)x/libsizes*1000000))

#B)median of ratios method from DESeq2:
#equivalent to do:
scaled_counts_DESeq2 <- counts(ddsMat,normalized=TRUE)

#C)
vsd<- vst(ddsMat, blind = FALSE) 
scaled_counts_vst <- assay(vsd)


library("pheatmap")                                                                                                                       
library("RColorBrewer")  

pca_both <- plotPCA(vsd,intgroup=c("patient","LPS"))+scale_color_discrete(name = "patient.LPS")                                                                                              
print(pca_both) +theme_classic()
pca_both$data

#to check what plotPCA does:
DESeq2:::plotPCA.DESeqTransform

pca_manual=prcomp(t(scaled_counts_vst))
pca_manual$x
View(pca_both$data)
rv <- rowVars(scaled_counts_vst)
ntop=500
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                   length(rv)))]
pca_manual=prcomp(t(scaled_counts_vst[select,]))


#
pcadf <- as.data.frame(pca_manual$x[,1:2])
pcadf <- cbind(pcadf,coldata[match(rownames(pcadf),coldata$patient.LPS),])
ggplot(pcadf,aes(x=PC1,y=PC2,color=patient.LPS))+geom_point(size=2)+
  theme_classic()
ggplot(pcadf,aes(x=PC1,y=PC2,color=Cohesin))+geom_point(size=2)+
  theme_classic()
#
pca_manual_CPM=prcomp(t(scaled_counts_CPM[select,]),scale. = T)
pcadf <- as.data.frame(pca_manual_CPM$x[,1:2])
pcadf <- cbind(pcadf,coldata[match(rownames(pcadf),coldata$patient.LPS),])
ggplot(pcadf,aes(x=PC1,y=PC2,color=patient.LPS))+geom_point(size=2)+
  theme_classic()
ggplot(pcadf,aes(x=PC1,y=PC2,color=Cohesin))+geom_point(size=2)+
  theme_classic()

pca_manual_log2CPM <- prcomp(t(log2(scaled_counts_CPM[select,]+0.01)))
pcadf <- as.data.frame(pca_manual_log2CPM$x[,1:2])
pcadf <- cbind(pcadf,coldata[match(rownames(pcadf),coldata$patient.LPS),])
ggplot(pcadf,aes(x=PC1,y=PC2,color=patient.LPS))+geom_point(size=2)+
  theme_classic()
ggplot(pcadf,aes(x=PC1,y=PC2,color=Cohesin))+geom_point(size=2)+
  theme_classic()
#variables: number of genes to take into account 
#DEseq2 takes the top 500 more variable
#normalization method

pca_cpm=prcomp(t(scaled_counts_CPM))
ggplot(as.data.frame(pca_cpm$x),aes(x=PC1,y=PC2))+geom_point(size=2)+
  theme_classic()

pca_cpm=prcomp(t(scaled_counts_CPM),scale. = T)
ggplot(as.data.frame(pca_cpm$x),aes(x=PC1,y=PC2))+geom_point(size=2)+
  theme_classic()

pca_cpm=prcomp(t(log2(scaled_counts_CPM+0.01)))
ggplot(as.data.frame(pca_cpm$x),aes(x=PC1,y=PC2))+geom_point(size=2)+
  theme_classic()


pca_deseq2_log=prcomp(t(log2(scaled_counts_DESeq2+0.01)))
ggplot(as.data.frame(pca_deseq2_log$x),aes(x=PC1,y=PC2))+geom_point(size=2)+
  theme_classic()
