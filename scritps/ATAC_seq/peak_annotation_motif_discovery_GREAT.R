## ---------------------------
##
##
## Purpose of script: perform peak annotation of ATAC-seq MDS
##        and prepare data for motif enrichment and GREAT analysis  
## 
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-18
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
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## ---------------------------

peakfile <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks=read.table(peakfile,header = T)
peaks$peakID=paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)

BiocManager::install("org.Hs.eg.db")
peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb,annoDb = "org.Hs.eg.db" )
View(as.data.frame(peakAnno))

#add some more info: translate geneID, add p-values
#gene_translation_file <- "~/Nextcloud/references/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
#gene_translation_file <- "~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/biomart_tables/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
#gene_translation <- fread(gene_translation_file)
peakAnnodf <- as.data.frame(peakAnno)
# table(peakAnnodf$geneId%in%gene_translation$`NCBI gene (formerly Entrezgene) ID`)
# peakAnnodf$gene_name <- gene_translation$`Gene name`[match(peakAnnodf$transcriptId,
#                                                            gene_translation$`Transcript stable ID version`)]
#I omit this step now as I added the option annoDb = "org.Hs.eg.db":
# org.Hs.egACCNUM Map Entrez Gene identifiers to GenBank Accession Numbers
# Description
# org.Hs.egACCNUM is an R object that contains mappings between Entrez Gene identifiers and
# GenBank accession numbers

peakAnnodf$peakID=paste0(peakAnnodf$seqnames,":",peakAnnodf$start - 1,"-",peakAnnodf$end)

## Write bed files for up and down genes
#Load deseq2 results
#3) differential peaks between Cohesin vs non-Cohesin patients

#deseqres=read.csv("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/all_samples/cohesin_vs_nocohesin/all_samples_cohesin_vs_nocohesin_DESeq2_results.csv")
deseqres=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/all_samples/cohesin_vs_nocohesin/all_samples_cohesin_vs_nocohesin_DESeq2_results.csv")

up_peaks <- deseqres$X[deseqres$padj<=0.05&deseqres$log2FoldChange>0]
down_peaks <- deseqres$X[deseqres$padj<=0.05&deseqres$log2FoldChange<0]

bgpeaks <- deseqres$X[!deseqres$X%in%c(up_peaks,down_peaks)]

#write input files
target_peaks_UP <- peaks[peaks$peakID%in%up_peaks,]
target_peaks_DOWN <- peaks[peaks$peakID%in%down_peaks,]
background_peaks <- peaks[peaks$peakID%in%bgpeaks,]

target_out_UP=target_peaks_UP[,c(1,2,3,7,5,4)]
target_out_UP$nsamples <- "."

target_out_DOWN=target_peaks_DOWN[,c(1,2,3,7,5,4)]
target_out_DOWN$nsamples <- "."

bg_out=background_peaks[,c(1,2,3,7,5,4)]
bg_out$nsamples <- "."

#write.table(target_out_UP,"in_data/879UP_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")
#write.table(target_out_DOWN,"in_data/430DOWN_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")
#write.table(bg_out,"in_data/bakground_30054nondiff_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")

#with peaks as background
#/Users/llorenzi/software/homer/bin/findMotifsGenome.pl in_data/879UP_peaks_Cohesin.bed hg38 homer_out_UPCohesin_peakBG -size given -bg in_data/bakground_30054nondiff_peaks_Cohesin.bed
#/Users/llorenzi/software/homer/bin/findMotifsGenome.pl in_data/430DOWN_peaks_Cohesin.bed hg38 homer_out_DOWNCohesin_peakBG -size given -bg in_data/bakground_30054nondiff_peaks_Cohesin.bed

#random background
#/Users/llorenzi/software/homer/bin/findMotifsGenome.pl in_data/879UP_peaks_Cohesin.bed hg38 homer_out_UPCohesin_randomBG -size given
#/Users/llorenzi/software/homer/bin/findMotifsGenome.pl in_data/430DOWN_peaks_Cohesin.bed hg38 homer_out_DOWNCohesin_randomBG -size given

#I will try the same analysis with GREAT. For great the background must contain the target peaks:

target_out_UP_GREAT=target_peaks_UP[,c(1,2,3,7,6,4)]
target_out_UP_GREAT$nsamples <- "."
target_out_UP_GREAT$meanscore <- round(target_out_UP_GREAT$meanscore)

target_out_DOWN_GREAT=target_peaks_DOWN[,c(1,2,3,7,6,4)]
target_out_DOWN_GREAT$nsamples <- "."
target_out_DOWN_GREAT$meanscore <- round(target_out_DOWN_GREAT$meanscore)

great_background=peaks[peaks$peakID%in%deseqres$X,]
bg_out=great_background[,c(1,2,3,7,6,4)]
bg_out$nsamples <- "."
bg_out$meanscore <- round(bg_out$meanscore)

# getwd()
# setwd("~/Nextcloud/MDS/ATAC-seq/analyses/")
# dir.create("GREAT")
# setwd("GREAT/")
# dir.create("in_data")
# head(bg_out)
# write.table(bg_out,"in_data/bakground_alltested_peaks_Cohesin.GREAT.bed",row.names = F,col.names = F,quote = F,sep = "\t")
# write.table(target_out_UP_GREAT,"in_data/879UP_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(target_out_DOWN_GREAT,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/GREAT/in_data/430DOWN_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")

##Generate GREAT http requests
#1st step is to upload bed files to dropbox. This is necessary to get url links to the files

#I will generate two types of link:
#1 batch output to a tsv file
#2 web output

setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/GREAT/")

dblinks <- read.table("in_data/DB_links.txt",header = T)
#pat=dblinks$V1[1]
#link=dblinks$V2[1]
cmm1 <- "http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=" #batch or web
cmm2 <- "&requestSpecies=hg38&requestName="
cmm3 <- "&requestURL="
cmm4 <- "&bgURL="
comms <- c()
bglink=dblinks[dblinks$type_set=="background","Dropbox_link"]
rownames(dblinks) <- dblinks$Peak_set_name
for(optype in c("batch","web")){
  for (nam in dblinks$Peak_set_name[dblinks$type_set=="test"]) {
    #samp_counts=samp_counts+1
    dblink <- dblinks[nam,"Dropbox_link"]
    comm=paste0(cmm1,
                optype,
                cmm2,
                nam,
                cmm3,
                gsub("/","%2F",gsub(":","%3A",dblink)),
                cmm4,
                gsub("/","%2F",gsub(":","%3A",bglink))
                )
    
    if(optype=="batch"){
      comm <- paste0("wget -O ",nam,".GREAT.tsv ",'"',comm,'"')
      
      
    }
    comms <- c(comms,
               comm)
  }
  write(comms,paste0("run_GREAT_Cohesin_UPandDOWN_peaks_HTTP_request_",optype,".txt"))
  comms=c()
}


#Analize peak annotation
filteredPeaksAnno <- peakAnnodf[peakAnnodf$peakID%in%deseqres$X,]
filteredPeaksAnno$DA <- "UNCHANGED"
filteredPeaksAnno$DA[filteredPeaksAnno$peakID%in%up_peaks] <- "UP"
filteredPeaksAnno$DA[filteredPeaksAnno$peakID%in%down_peaks] <- "DOWN"

table(filteredPeaksAnno$DA)
filteredPeaksAnno$simplified_annot <- filteredPeaksAnno$annotation
filteredPeaksAnno$simplified_annot[grep("Exon",filteredPeaksAnno$annotation)] <- "exon"
filteredPeaksAnno$simplified_annot[grepl("Exon",filteredPeaksAnno$annotation)&
                                           grepl("exon 1 of",filteredPeaksAnno$annotation)] <- "exon (1st)"

filteredPeaksAnno$simplified_annot[grep("Intron",filteredPeaksAnno$annotation)] <- "intron"

filteredPeaksAnno$simplified_annot[grepl("Intron",filteredPeaksAnno$annotation)&
                                     grepl("intron 1 of",filteredPeaksAnno$annotation)] <- "intron (1st)"

cbPalette <- c("#999999", 
               "#E69F00", 
               "#56B4E9", 
               "#009E73", 
               "#F0E442", 
               "#0072B2", 
               "#D55E00", 
               "#CC79A7") #source http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

library(RColorBrewer)
levels(as.factor(filteredPeaksAnno$simplified_annot))
topie <- table(filteredPeaksAnno$simplified_annot)
pie(topie,col = colorRampPalette( cbPalette)(length(topie)))

all_filtered_peaks <- filteredPeaksAnno$peakID

peak_set <- list(ALL_filtered_peaks=all_filtered_peaks,
                 UP_peaks=up_peaks,
                 DOWN_peaks=down_peaks)
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/")
dir.create("peak_annotation")
setwd("peak_annotation/")

dir.create("filtered_peaks")
setwd("filtered_peaks/")
dir.create("peaks_in_half_of_samples_and_lowqval")
setwd("peaks_in_half_of_samples_and_lowqval/")

for (set in names(peak_set)) {
  topie <- table(filteredPeaksAnno$simplified_annot[filteredPeaksAnno$peakID%in%
                                                      peak_set[[set]]])
  
  pdf(paste0(set,"peak_annot_pie.pdf"))
  pie(topie,col = colorRampPalette( cbPalette)(length(topie)),main = set)
  dev.off()
  
}
#the same but plot all together
topie=table(filteredPeaksAnno$simplified_annot)
cols=colorRampPalette( cbPalette)(length(topie))
names(cols)=names(topie)
pdf("filtered_peak_annot_pie.pdf")
par(mfrow=c(1,3))
for (set in names(peak_set)) {
  topie <- table(filteredPeaksAnno$simplified_annot[filteredPeaksAnno$peakID%in%
                                                      peak_set[[set]]])
  pie(topie,col = cols[names(topie)],main = set)
  
  
  
}
dev.off()
peakAnno