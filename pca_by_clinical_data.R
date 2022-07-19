#######################  2. PCA plot   ###########################
#I will make PCAs based on gene expression
# 1) Read in expression data (raw counts)

counts <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv")
#firstbatch <- read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/first_batch_gencodev38.counts.csv")
#secondbatch  <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/second_batch_gencodev38.counts.csv")
colnames(counts)
colnames(counts) <- gsub("X","",colnames(counts))
#convert counts into a numeric data frame
#for this assign gene ids as rownames and remove this column
rownames(counts) <- counts$gene_id
counts <- counts[,-1]
# 2) keep only  genes that have more than 1 count 
#(this is a very basic filtering though)
rowSums(counts)
libsize=colSums(counts)
sort(libsize/1000000)

table(rowSums(counts)>1)

filtered_counts <- counts[rowSums(counts)>=1,]
#we can overwrite counts with filtered_counts as we won't use filtered_counts
counts <- filtered_counts
#we could do both steps in one:
counts <- counts[rowSums(counts)!=0,]

# 3) Normalize data. 
#Useful reading: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

#in order to account for differences in library sizes (one of the main
#non-biological differences between samples) and other
#non-biological factors of variance we need to normalize the raw counts
#for PCA I usually use the normalizing method recommended by the
#DESeq2 package, variance stabilizing transformation
#This already accounts for lib sizes and is a kind of log transformation
#that takes into account the dispersion trend with higher expression
#and corrects for this. 
#Although now I have seen some people argue 
#this is not a good method: https://support.bioconductor.org/p/108321/

library(DESeq2)
class(counts)
vst_counts <- vst(as.matrix(counts))

# 4) compute principal component analysis using base R function prcomp
#this is the function that DESeq2 actually uses
rv <- rowVars(vst_counts)
ntop=500 #select only the top 500 genes with highest variance
#This is optional, is the default of the DESeq2 package
select <- order(rv, decreasing = TRUE)[1:ntop]
vst_top500variable <- vst_counts[select,]

pca_results <- prcomp(t(vst_top500variable)) #IMPORTANT: we have to take the transpose of our data
#becasue this function expects variables as columns!!!
class(pca_results)
pca_results$x
View(pca_results$x)
summary(pca_results)
data_to_plot <- as.data.frame(pca_results$x[,1:2])

# 5) We are ready to do our basic ugly plot
# plot(data_to_plot$PC1,
#      data_to_plot$PC2,
#      pch=20)

# 6) read sample metadata to add meaningful colors to sample features
processed_samples.metadata <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv")
View(processed_samples.metadata)
#let's create a variable to match our sample names in counts with metadata
processed_samples.metadata$patient.LPS <- paste(processed_samples.metadata$patient,
                                                processed_samples.metadata$LPS,sep = ".")
table(colnames(vst_counts)%in%processed_samples.metadata$patient.LPS)
#7) bind metadata to our data_to_plot data frame
data_to_plot <- cbind(data_to_plot,
                      processed_samples.metadata[match(rownames(data_to_plot),
                                                       processed_samples.metadata$patient.LPS),])

#add clinical metadata
clinical_data <- readxl::read_xlsx("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/Analisis_Pame/Datos clinicos_pacientes SC.xlsx",sheet = 2)
table(processed_samples.metadata$patient%in%clinical_data$`Case ID`)
data_to_plot <- cbind(data_to_plot,clinical_data[match(data_to_plot$patient,
                                                       clinical_data$`Case ID`),])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/plots/")
dir.create("PCA_clinical_data")
setwd("PCA_clinical_data/")
data_to_plot$Cohesin[data_to_plot$Cohesin] <- 1 



library(ggplot2)
data_to_plot$Cohesin <- as.factor(data_to_plot$Cohesin)
data_to_plot$LPS <- as.factor(data_to_plot$LPS)
ggplot(data_to_plot,aes(x = PC1, y = PC2)) + geom_point(aes(col=Cohesin,shape=LPS)) +
  theme_classic() + scale_color_manual(values = c("#F0E442", "darkgreen")) +
  theme(line = element_line(size = 1), plot.title = element_text(hjust = 0.5) ) + ylim(c(-75,30)) + xlim(c(-75,30)) + 
  ggtitle("Hola PCA")

feat="`Cytogenetc Risk`"
feat="`IPSS-R`"
feat="`WHO 2017`"

colnames(data_to_plot) <- gsub(" \\[.*\\]","",colnames(data_to_plot))
colnames(data_to_plot) <- gsub(" ","_",colnames(data_to_plot))
colnames(data_to_plot) <- gsub("-","_",colnames(data_to_plot))
data_to_plot$Gender <- tolower(data_to_plot$Gender)
data_to_plot$Cytogenetc_Risk <- tolower(data_to_plot$Cytogenetc_Risk)
data_to_plot$IPSS_R <- tolower(data_to_plot$IPSS_R)
data_to_plot$CPSS <- tolower(data_to_plot$CPSS)
data_to_plot$combn_riskscores <- data_to_plot$IPSS_R
cond=is.na(data_to_plot$IPSS_R)|data_to_plot$IPSS_R=="-"
data_to_plot$combn_riskscores[cond] <- data_to_plot$CPSS[cond]

feats_to_plot <- colnames(data_to_plot)[c(9,10,13:19,21:25)]
data_to_plot$IPSS_R_Score <- as.numeric(data_to_plot$IPSS_R_Score)

pdf("plots_vst_top500moreVariableGenes.pdf",onefile = T)
for (feat in feats_to_plot) {
  g=ggplot(data_to_plot,aes(x = PC1, y = PC2)) + 
    geom_point(aes_string(col=feat,shape="LPS")) +
    theme_classic() 
  print(g)
}
dev.off()
# + scale_color_manual(values = c("#F0E442", "darkgreen")) +
#   theme(line = element_line(size = 1), plot.title = element_text(hjust = 0.5) ) + 
#   ylim(c(-75,30)) + xlim(c(-75,30)) + 
#   ggtitle(feat)


####################Try log2 normalized data###############
library(tidyverse)
libsizes <- colSums(counts)
cpm_counts <- t(apply(counts, 1, function(x)x/libsizes*1000000))

df <- bind_rows(as_data_frame(log2(cpm_counts +1)) %>%mutate(transformation = "log2(x + 1)"),
                as_data_frame(vst_counts) %>% mutate(transformation = "vst"))


lvls <- c("log2(x + 1)", "vst")
df$transformation <- factor(df$transformation, levels=lvls)
colnames(df)[1:2] <- c("x","y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

log2_xplus1 <- log2(cpm_counts+1)

rv <- rowVars(log2_xplus1)
ntop=500 #select only the top 500 genes with highest variance
#This is optional, is the default of the DESeq2 package
select <- order(rv, decreasing = TRUE)[1:ntop]
log2_xplus1_top500variable <- log2_xplus1[select,]

pca_results_log2 <- prcomp(t(log2_xplus1_top500variable)) #IMPORTANT: we have to take the transpose of our data

data_to_plot_log2cpm <- data_to_plot
data_to_plot_log2cpm[,1:2] <- as.data.frame(pca_results_log2$x[,1:2])

pdf("plots_log2cpm_top500moreVariableGenes.pdf",onefile = T)
for (feat in feats_to_plot) {
  g=ggplot(data_to_plot_log2cpm,aes(x = PC1, y = PC2)) + 
    geom_point(aes_string(col=feat,shape="LPS")) +
    theme_classic() 
  print(g)
}
dev.off()
