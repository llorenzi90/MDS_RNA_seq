options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------
countdata=read.csv("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/htseq_counts/all_samples_gencodev19.counts.gene_name.csv")
rownames(countdata) <- countdata$gene_name
countdata <- countdata[,-1]

spcols=strsplit(colnames(countdata),split = "_")
patient <- sapply(spcols, function(x)x[1])
treatment <- sapply(spcols, function(x)x[3])


samples_info <- fread("~/shares/BDCuartero/RAW_DATA/EGA/info.txt")
table(samples_info$Gene.refGene)
samples_info_muts <- samples_info[samples_info$Gene.refGene!="WHO_subtype",]
samples_info_whost <- samples_info[samples_info$Gene.refGene=="WHO_subtype",]
WHO_subtype <- samples_info_whost$ExonicFunc.refGene[match(patient,samples_info_whost$sample)]



#Classify genes by pathwaty using the Leukemia 2014 paper as reference
S3table <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/OncoPrint/Suppl_TableS3_Leukemia_2014.csv")
S3table_cohesin_genes <- (S3table %>% filter(Pathway=="Cohesin") %>%select(Gene))[,1]
Cohesin_samples <- unique(samples_info$sample[samples_info$Gene.refGene%in%S3table_cohesin_genes])


#But not all genes in EGA info table are in the S3 table
unique(samples_info_muts$Gene.refGene)[!unique(samples_info_muts$Gene.refGene)%in%S3table$Gene]
genes_not_in_S3_table <- unique(samples_info_muts$Gene.refGene)[!unique(samples_info_muts$Gene.refGene)%in%S3table$Gene]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


bm <- getBM(attributes = c("phenotype_description",
                           "external_gene_name",
                           "hgnc_symbol","external_synonym"
),
mart = ensembl)




table(genes_not_in_S3_table%in%c(bm$hgnc_symbol,bm$external_synonym))
genes_symbol_not_in_S3_table <- sapply(genes_not_in_S3_table, function(z){
  ifelse(z%in%bm$hgnc_symbol,
         return(bm$hgnc_symbol[match(z,bm$hgnc_symbol)]),
         ifelse(z%in%bm$external_synonym,return(bm$external_synonym[match(z,bm$external_synonym)]),return(NA)))
})
genes_symbol_not_in_S3_table <- genes_symbol_not_in_S3_table[!is.na(genes_symbol_not_in_S3_table)]
possible_gene_fusions=genes_not_in_S3_table[!genes_not_in_S3_table%in%c(bm$hgnc_symbol,bm$external_synonym)]
possible_gene_fusions <- strsplit(possible_gene_fusions,split = "-")
possible_gene_fusions <- possible_gene_fusions[sapply(possible_gene_fusions, function(x)length(x)==2)]
tt=lapply(possible_gene_fusions,function(x){
  sapply(x, function(z){
    ifelse(z%in%bm$hgnc_symbol,
           return(bm$hgnc_symbol[match(z,bm$hgnc_symbol)]),
           ifelse(z%in%bm$external_synonym,return(bm$external_synonym[match(z,bm$external_synonym)]),return(NA)))
  })
})

# go_list <- getBM(attributes=c("go_id", "name_1006", "namespace_1003","hgnc_symbol"),
#                  
#                  mart=ensembl
# )
#BiocManager::install("GOfuncR")
library("GOfuncR")
#annocat=get_anno_categories(genes_symbol_not_in_S3_table,  database = 'Homo.sapiens')
#BiocManager::install('Homo.sapiens')
#library("Homo.sapiens")
#annocat=get_anno_categories(genes_symbol_not_in_S3_table,database = "org.Hs.eg.db")
annocat=get_anno_categories(c(bm$hgnc_symbol,bm$external_synonym),database = "org.Hs.eg.db")

length(genes_symbol_not_in_S3_table)

tt_gene_symbol <- lapply(possible_gene_fusions,function(x){
  sapply(x, function(z){
    ifelse(z%in%bm$hgnc_symbol,
           return(bm$hgnc_symbol[match(z,bm$hgnc_symbol)]),
           ifelse(z%in%bm$external_synonym,return(bm$hgnc_symbol[match(z,bm$external_synonym)]),return(NA)))
  })
})
tt_gene_symbol
genes_symbol_not_in_S3_table_symbol <- sapply(genes_not_in_S3_table, function(z){
  ifelse(z%in%bm$hgnc_symbol,
         return(bm$hgnc_symbol[match(z,bm$hgnc_symbol)]),
         ifelse(z%in%bm$external_synonym,return(bm$hgnc_symbol[match(z,bm$external_synonym)]),return(NA)))
})

genes_symbol_not_in_S3_table_symbol <- genes_symbol_not_in_S3_table_symbol[!is.na(genes_symbol_not_in_S3_table_symbol)]

#assign gene ontology to each gene
#any cohesin?
table(annocat$name[annocat$gene%in%c(genes_symbol_not_in_S3_table_symbol,
                                     unlist(tt_gene_symbol))])
gene_ontologies_unclassified_genes <- annocat[annocat$gene%in%c(genes_symbol_not_in_S3_table_symbol,
                                                                unlist(tt_gene_symbol)),]
gene_ontologies_ug_matching_cohesin <- gene_ontologies_unclassified_genes[grep("cohesin",
                                                                               gene_ontologies_unclassified_genes$name,
                                                                               ignore.case = T),]

#there is one potential gene that might be cohesin related
#"RB1" although is not a component of cohesin complex
#we DO NOT consider this gene as Cohesin mutant

mutations <- sapply(patient,function(x)paste(samples_info_muts$Gene.refGene[samples_info_muts$sample==x],collapse = ","))
cohesin <- patient%in%Cohesin_samples

coldata=data.frame(cnames=colnames(countdata),
                   patient=patient,
                   treatment=treatment,
                   WHO_subtype=WHO_subtype,
                   mutations=mutations,
                   cohesin_mut=cohesin)
coldata$N_muts <- sapply(coldata$mutations,function(x)length(strsplit(x,split = ",")[[1]]))
#conversion WHO_subtypes old new
WHO_old_new <- c("RCUD"="MDS-SLD",
                 "RARS"="MDS-RS-SLD",
                 "RCMD"="MDS-MLD",
                 "RCMD-RS"="MDS-RS-MLD",
                 "RAEB-1"="MDS-EB-1",
                 "RAEB-2"="MDS-EB-2",
                 "RARS-T"="MDS/MPN-RS-T")

table(coldata$WHO_subtype%in%names(WHO_old_new))
coldata$WHO_subtype_old <- coldata$WHO_subtype
coldata$WHO_subtype[coldata$WHO_subtype%in%names(WHO_old_new)] <- WHO_old_new[match(coldata$WHO_subtype[coldata$WHO_subtype%in%names(WHO_old_new)],
                                                                                    names(WHO_old_new))]

#clean up data:

#A) General clean-up
#     remove patients wihtout WHO_subtype and those with "treatment" (not sure what they are)

coldata$WHO_subtype[grep("normal",coldata$cnames)] <- "healthy"

write.csv(coldata,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/metadata_RNA_seq_counts.csv",row.names = F)
#remove NA for WHO_subtype
filtered_coldata <- coldata[!is.na(coldata$WHO_subtype),]
filtered_coldata <- filtered_coldata[!(!is.na(filtered_coldata$treatment)&filtered_coldata$treatment=="CHX"),]
length(unique(filtered_coldata$patient))


filtered_countdata <- countdata[,colnames(countdata)%in%filtered_coldata$cnames]
table(colnames(filtered_countdata)%in%coldata$cnames)


counts <- filtered_countdata[rowSums(filtered_countdata)!=0,]

library(DESeq2)
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
table(filtered_coldata$cnames==rownames(data_to_plot))

data_to_plot <- cbind(data_to_plot,filtered_coldata[,c("patient","WHO_subtype","cohesin_mut","N_muts")])

library(ggplot2)
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")
data_to_plot$cohesin_mut <- as.factor(data_to_plot$cohesin_mut)

dir.create("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/plots")
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/plots")

pdf("PCA_CD34_samples_CohesinMutandWHO_subtype.pdf")
g=ggplot(data_to_plot,aes(x = PC1, y = PC2)) + geom_point(aes(col=WHO_subtype,shape=cohesin_mut)) +
  theme_classic()  + scale_color_manual(values=c(cbPalette,"red","black","green","blue"))+
  theme(line = element_line(size = 1), plot.title = element_text(hjust = 0.5) ) 
print(g)
dev.off()

pdf("PCA_CD34_samples_CohesinMutandNmutations.pdf")
g=ggplot(data_to_plot,aes(x = PC1, y = PC2)) + geom_point(aes(col=as.factor(N_muts),shape=cohesin_mut)) +
  theme_classic()  + scale_color_manual(values=c(cbPalette,"black")) +
  theme(line = element_line(size = 1), plot.title = element_text(hjust = 0.5) )
print(g)
dev.off()

#B) comparisons to do: 1) remove healthy and remove all patients that have no annotated genes mutated
#                   2) remove only healthy

#1) remove healthy and remove all patients that have no annotated genes mutated

colda <- filtered_coldata[filtered_coldata$N_muts!=0,]
countda <- counts[,colnames(counts)%in%colda$cnames]
table(colnames(countda)==colda$cnames)
ddsmat <- DESeqDataSetFromMatrix(countData = countda,                                                                                   
                                 colData = colda,                                                                                       
                                 design = ~ cohesin_mut )   

dir.create("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/DESEq2")
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/DESEq2")
dir.create("CohesinMut_vs_nonCohesinMut")
setwd("CohesinMut_vs_nonCohesinMut/")
keep <- rowSums(counts(ddsmat)) > 1
ddsMat <- ddsmat[keep,]

ddsMat <- estimateSizeFactors(ddsMat)
dds <- DESeq(ddsMat)
resultsNames(dds)
res <- results(dds,alpha = 0.05)
write(capture.output(resultsNames(dds)),"summary_res.txt")
write(capture.output(summary(res)),"summary_res.txt",append = T)
res <- as.data.frame(res)
write.csv(res[order(res$padj),],"DESeq2_results.csv")

rank <- as.data.frame(res[order(res$stat,decreasing = T),])
rank$gene <- rownames(rank)
# rank$gene[rank$gene%in%conversion_ENSMBL_NCBI$Gene.stable.ID.version] <- 
#   conversion_ENSMBL_NCBI$Gene.name[match(rank$gene[rank$gene%in%conversion_ENSMBL_NCBI$Gene.stable.ID.version],
#                                          conversion_ENSMBL_NCBI$Gene.stable.ID.version)]
dir.create("rank_files")
write.table(rank[,c("gene","stat")],"rank_files/DESeq2.rnk",quote = F,col.names = F,row.names = F,sep = "\t")
