## ---------------------------
##
##
## Purpose of script: genertae TPMs and FPKMs from counts 
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-11
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
count_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv"
gene_lengths_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/annotation/gencode/gencode.v38.gene_lengths.tsv"
conversion_table_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/biomart_tables/Human_EnsemblIDs.104.GRCh38.p13_to_gene_names.tsv"

out_file_prefix <- gsub(".counts.csv","",count_file)
countdata <- read.csv(count_file)

### Normalize data ###
#remove rownames to work for now
gene_ids <- countdata$gene_id
countdata <- countdata[,-1]
colnames(countdata) <- gsub("X","",colnames(countdata))

#read gene lengths and gene ids/gene names table
gene_lengths <- read.table(gene_lengths_file)
table(gene_ids==gene_lengths$V1)
geneid_gene_names <- read.table(conversion_table_file)
table(gene_ids%in%geneid_gene_names$Gene.stable.ID.version)

gene_names <- geneid_gene_names$Gene.name[match(gene_ids,
                                                geneid_gene_names$Gene.stable.ID.version)]
length(unique(gene_names[!is.na(gene_names)]))
length(gene_names[!is.na(gene_names)])

gene_names[is.na(gene_names)] <- gene_ids[is.na(gene_names)]
anyNA(gene_names)
duplicated_gene_names <- gene_names[duplicated(gene_names)]


ord=order(duplicated_gene_names)
ndups=table(duplicated_gene_names)
table(unique(duplicated_gene_names[ord])==names(ndups))

dup_vect <- sapply(ndups, function(x)return(paste0("_",seq(1:x))))
modified_names <- duplicated_gene_names
modified_names[ord] <- paste0(duplicated_gene_names[ord],unlist(dup_vect))

gene_names[duplicated(gene_names)] <- modified_names
anyNA(gene_names)
any(duplicated(gene_names))


# 
##TPM
# 
##FPKM
#These three metrics attempt to normalize for sequencing depth and gene length. Here’s how you do it for RPKM (or FPKM):
#   
#1)Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
pmsf <- apply(countdata, 2,sum)/1000000
#2)Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
rpm <- t(apply(countdata,1,function(x)x/pmsf))
head(rpm)
#3)Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.
kilobases <- gene_lengths$V2/1000
FPKM <- apply(rpm,2,function(x)x/kilobases )
# 
# ########compare to result using fpkm function from DESeq2#######
# library(DESeq2)
# #extrafont::loadfonts()
# countdata=read.csv("data/all_samples_gencodevM27.counts.csv")
# rownames(countdata) <- countdata$gene_id
# countdata <- countdata[,-1]
# name=colnames(countdata)
# tr=strsplit(name,split = "\\.")
# vector=sapply(tr, function(x)return(x[1]))
# rep=sapply(tr, function(x)return(x[length(x)]))
# LPS=sapply(tr, function(x)return(any(grepl("LPS",x))))
# coldata=data.frame(name,vector,LPS,rep)
# rownames(coldata) <- coldata$name
# table(rownames(coldata)==colnames(countdata))
# coldata$vector <- as.factor(coldata$vector)
# coldata$vector
# class(coldata$LPS)
# coldata$LPS <- as.factor(coldata$LPS)
# coldata$LPS
# 
# ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
#                                  colData = coldata,
#                                  design = ~0+ vector + LPS + vector:LPS)
# 
# mcols(ddsMat)$basepairs <- gene_lengths$V2
# test <- fpkm(object=ddsMat,robust = F)
# 

FPKM_gene_name <- cbind(gene_name=gene_names,FPKM)
write.csv(FPKM_gene_name,paste0(out_file_prefix,".gene_name.FPKM.csv"),row.names = F,quote = F)

FPKM_gene_id <- cbind(gene_id=gene_ids,FPKM)
write.csv(FPKM_gene_id,paste0(out_file_prefix,".geneID.FPKM.csv"),row.names = F,quote = F)
          
# 
# 
#####TPM####
# countdata <- read.csv("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/CEBPa/RNA-seq/data/expression_files/all_samples_gencodevM27.htseq_counts.gene_name.csv")
# rownames(countdata) <- countdata$gene_name
# countdata <- countdata[,-1]
kilobases <- gene_lengths$V2/1000

#1)Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
rpk <- apply(countdata, 2, function(x)x/kilobases)
#2)Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
pmsf.1 <- apply(rpk, 2, sum)/1000000
#3)Divide the RPK values by the “per million” scaling factor.
TPM <- t(apply(rpk, 1, function(x)x/pmsf.1))
# 
TPM_gene_name <- cbind(gene_name=gene_names, TPM)
TPM_gene_id <- cbind(gene_id=gene_ids,TPM)

write.csv(TPM_gene_name,paste0(out_file_prefix,".gene_name.TPM.csv"),row.names = F,quote = F)

write.csv(TPM_gene_id,paste0(out_file_prefix,".geneID.TPM.csv"),row.names = F,quote = F)
