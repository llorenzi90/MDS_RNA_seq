## ---------------------------
##
##
## Purpose of script: Classify EGA patients and select data to copy to csuc
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-09
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
samples_info <- fread("~/shares/BDCuartero/RAW_DATA/EGA/info.txt")

list_of_bam_files <- read.table("~/shares/BDCuartero/RAW_DATA/EGA/list_bam_files.txt",sep = "/")
samp_elements <-  unlist(strsplit(list_of_bam_files$V3,split = "_"))
table(samp_elements[!samp_elements%in%samples_info$sample])
list_of_bam_files <- list_of_bam_files %>% separate(col = V3,sep = "_",into = c("sample","cell_treat"),remove = F,extra = "merge")
table(unique(list_of_bam_files$sample)%in%samples_info$sample)

View(as.data.frame(table(samples_info$Gene.refGene,samples_info$sample)))
length(unique(samples_info$sample))
#Classify genes by pathwaty using the Leukemia 2014 paper as reference
S3table <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/OncoPrint/Suppl_TableS3_Leukemia_2014.csv")
table(unique(samples_info$Gene.refGene)%in%S3table$Gene)
#there are 44 genes that are in the paper table
#how many of these are Cohesin genes?
S3table_cohesin_genes <- (S3table %>% filter(Pathway=="Cohesin") %>%select(Gene))[,1]
unique(samples_info$sample[samples_info$Gene.refGene%in%S3table_cohesin_genes])
Cohesin_samples <- unique(samples_info$sample[samples_info$Gene.refGene%in%S3table_cohesin_genes])
table(Cohesin_samples%in%list_of_bam_files$sample)
View(list_of_bam_files[list_of_bam_files$sample%in%Cohesin_samples,])
unique(samples_info$Gene.refGene)[!unique(samples_info$Gene.refGene)%in%S3table$Gene]
View(list_of_bam_files[!list_of_bam_files$sample%in%samples_info$sample,])

#There are some genes that are not in the paper table.
#I will try to classify them using biomaRt
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listEnsembl()
searchFilters(ensembl,pattern = "phenotype")
searchAttributes(ensembl,pattern = "phenotype")
searchFilters(ensembl,pattern= "symbol")
searchAttributes(ensembl,pattern= "gene_name")
searchAttributes(ensembl,pattern= "synon")
searchAttributes(ensembl,pattern= "GO")


bm <- getBM(attributes = c("phenotype_description",
                           "external_gene_name",
                           "hgnc_symbol","external_synonym"
                           ),
            mart = ensembl)




genes_not_in_S3_table <- unique(samples_info$Gene.refGene)[!unique(samples_info$Gene.refGene)%in%S3table$Gene]
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

go_list <- getBM(attributes=c("go_id", "name_1006", "namespace_1003","hgnc_symbol"),
                 
                 mart=ensembl
)
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
gene_ontologies_ug_matching_cohesin

#there is one potential gene that might be cohesin related
#"RB1" although is not a component of cohesin complex
#would it be considered a cohesin mutant? Ask Sergi

#besides genes, there are other types of information in this table,
#including the WHO_subtype for each patient
#check what other fields, different to genes there are:
table(samples_info$Gene.refGene[!samples_info$Gene.refGene%in%c(bm$hgnc_symbol,bm$external_synonym)])
table(samples_info$ExonicFunc.refGene[samples_info$Gene.refGene=="WHO_subtype"])
#what type are Cohesin samples?
table(samples_info$ExonicFunc.refGene[samples_info$Gene.refGene=="WHO_subtype"&samples_info$sample%in%Cohesin_samples])
#how many of the Cohesin patients have cd34?
table(grepl("CD34",list_of_bam_files$cell_treat[list_of_bam_files$sample%in%Cohesin_samples]))
#7 samples
table(grepl("BMMNC",list_of_bam_files$cell_treat[list_of_bam_files$sample%in%Cohesin_samples]))
#all of them have BMMNC
#conversion WHO_subtypes old new
WHO_old_new <- c("RCUD"="MDS-SLD",
                 "RARS"="MDS-RS-SLD",
                 "RCMD"="MDS-MLD",
                 "RCMD-RS"="MDS-RS-MLD",
                 "RAEB-1"="MDS-EB-1",
                 "RAEB-2"="MDS-EB-2",
                 "RARS-T"="MDS/MPN-RS-T")


#select bams to copy
bams_to_copy <- list_of_bam_files %>%filter(sample%in%samples_info$sample )
sort(table(bams_to_copy$cell_treat))
View(bams_to_copy)
bams_to_copy <- bams_to_copy %>%filter(cell_treat%in%c("CD34.bam", "BMMNC.bam"))
length(unique(bams_to_copy$sample))
#199 samples there are 11 samples that are not in the selected bams
#why is that check what type of data had those samples 
bams_no_selected <- list_of_bam_files %>%filter(sample%in%samples_info$sample, !cell_treat%in%c("CD34.bam", "BMMNC.bam"))
#check those samples that are only in bams_no_selected
table(bams_no_selected$cell_treat[!bams_no_selected$sample%in%bams_to_copy$sample])
#ok so the non-selected samples are all samples that have only DNA bams available

write.table(apply(bams_to_copy[,1:3],1,function(x)paste(x,collapse = "/")),"~/shares/BDCuartero/RAW_DATA/EGA/files_to_copy.txt",quote = F,col.names = F,row.names = F)

#read sizes of files to copy:
sizes <- read.table("~/shares/BDCuartero/RAW_DATA/EGA/sizes_files_to_copy.txt")
table(grepl("G",sizes$V1))
sizes$size <- as.numeric(gsub(",",".",gsub("G","",sizes$V1)))
sum(sizes$size)
#I have only ~400 GB left in scratch to copy data and perform analyses
#So I should split the data in pieces of ~ 200 GB
sum(sizes$size)/9

#What if I only copy the samples that have CD34?
list_of_bam_files_CD34 <- list_of_bam_files[grepl("CD34",list_of_bam_files$cell_treat),]
sizes_CD34 <- read.table("shares/BDCuartero/RAW_DATA/EGA/sizes_CD34_bams.txt")
sizes_CD34 <- separate(sizes_CD34,col = 2,into = c(NA,"dir","subdir","bam"),sep = "/")
list_of_bam_files_CD34$sizes <- sizes_CD34$V1[match(list_of_bam_files_CD34$V3,
                                                 sizes_CD34$bam)]
list_of_bam_files_CD34$sizes <- gsub("G","",gsub(",",".",list_of_bam_files_CD34$sizes))
list_of_bam_files_CD34$sizes <- as.numeric(list_of_bam_files_CD34$sizes )
sum(list_of_bam_files_CD34$sizes)
set.seed(0303456)
ibams <- sample(seq_len(nrow(list_of_bam_files_CD34)),
               round(nrow(list_of_bam_files_CD34)/2))

sum(list_of_bam_files_CD34$sizes[ibams])
sum(list_of_bam_files_CD34$sizes[-ibams])
bams2sync1 <- unite(list_of_bam_files_CD34[ibams,1:3],col = "path",sep = "/")
bams2sync2 <- unite(list_of_bam_files_CD34[-ibams,1:3],col = "path",sep = "/")

write.table(bams2sync1,"~/shares/BDCuartero/RAW_DATA/EGA/bams2sync1.txt",quote = F,col.names = F,row.names = F)
write.table(bams2sync2,"~/shares/BDCuartero/RAW_DATA/EGA/bams2sync2.txt",quote = F,col.names = F,row.names = F)


#What to do with the data??
#These files are mapped to human genome hg19
#Should we just quantify using hg19 annotation?

