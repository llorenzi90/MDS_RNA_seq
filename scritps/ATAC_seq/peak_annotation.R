## ---------------------------
##
##
## Purpose of script:peak annotation MDS 
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-13
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

#peakfile <- "~/Nextcloud/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL"
peakfile <- "~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks=read.table(peakfile)
colnames(peaks) <- c( "seqnames",
                      "start",
                      
                      "end",
                      "nsamples",
                      "samples",
                      "meanscore")
#write.table(peaks,"~/Nextcloud/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames",sep = "\t",row.names = F,quote = F)
#write.table(peaks,paste0(peakfile,".withcolnames"),sep = "\t",row.names = F,quote = F)

#peakfile <- "~/Nextcloud/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames"
peakfile <- paste0(peakfile,".withcolnames")

peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
View(as.data.frame(peakAnno))

#add some more info: translate geneID, add p-values
#gene_translation_file <- "~/Nextcloud/references/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
gene_translation_file <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/biomart_tables/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
gene_translation <- fread(gene_translation_file)
peakAnnodf <- as.data.frame(peakAnno)
table(peakAnnodf$geneId%in%gene_translation$`NCBI gene (formerly Entrezgene) ID`)
peakAnnodf$gene_name <- gene_translation$`Gene name`[match(peakAnnodf$geneId,
                                                           gene_translation$`NCBI gene (formerly Entrezgene) ID`)]
peakAnnodf$peakID=paste0(peakAnnodf$seqnames,":",peakAnnodf$start - 1,"-",peakAnnodf$end)
write.csv(peakAnnodf,paste0(peakfile,".annotated.csv"),row.names = F)

#read in DESeq2 results
for (lps in c(1,2)) {
  deseq2resLPS=read.csv(paste0("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/LPS_",lps,"/LPS_",lps,"_DESeq2_results.csv"))
  #add annotation
  deseq2resLPS <- cbind(deseq2resLPS,peakAnnodf[match(deseq2resLPS$X,peakAnnodf$peakID),c(6:18)])
  #write annotated results
  write.csv(deseq2resLPS,paste0("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/LPS_",lps,"/LPS_",lps,"_DESeq2_results.annotated.csv"),row.names = F)
}

write.csv(deseq2resLPS[1:100,],paste0("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/LPS_",lps,"/LPS_",lps,"_DESeq2_results.annotated.top100peaks.csv"),row.names = F)



#filter annotated peaks
# peaks <- read.table('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames',header = T )
# peaks$peakID=paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)
# 
# #check distribution of q-vals for those peaks that are called in only 1 sample
# summary(peaks$meanscore[peaks$nsamples==1])
# plot(density(peaks$meanscore[peaks$nsamples==1]))
# points(density(peaks$meanscore[peaks$nsamples>=10]),col="red",type = "s")
# plot(density(peaks$meanscore[peaks$nsamples>=10]))
# summary(peaks$meanscore[peaks$nsamples>=10])
# 
# table(peaks$meanscore<0.0001&peaks$nsamples==1)
# 
# peaks$nsamps_factor=1
# peaks$nsamps_factor[peaks$nsamples<=14] <- "1-14"  
# peaks$nsamps_factor[peaks$nsamples>14&peaks$nsamps_factor<=28] <- "15-28"  
# peaks$nsamps_factor[peaks$nsamples>28&peaks$nsamps_factor<=42] <- "29-42"  
# peaks$nsamps_factor[peaks$nsamples>42&peaks$nsamps_factor<=56] <- "43-56"  
# peaks$nsamps_factor[peaks$nsamples>57&peaks$nsamps_factor<=70] <- "57-70"  
# peaks$nsamps_factor[peaks$nsamples>71&peaks$nsamps_factor<=84] <- "71-84"  
# 
# peaks$nsamps_factor <- as.factor(peaks$nsamps_factor)
# 
# ggplot(peaks,aes(x=nsamps_factor,y=meanscore))+geom_boxplot() + theme_classic() + xlab("# samples")
# ggplot(peaks,aes(x=nsamps_factor,y=meanscore))+geom_boxplot() +theme_classic()+ xlab("# samples") + ylab("mean -log10(q-val)")
# 
# coval=mean(peaks$meanscore[peaks$nsamps_factor=="15-28"])
#coval
#9.122563
coval=9.122563

peaks_to_keep <- peakAnnodf$nsamples>=42|peakAnnodf$meanscore>=coval
table(peaks_to_keep)
peakAnnodf_filtered <- peakAnnodf[peaks_to_keep,]
write.csv(peakAnnodf_filtered,"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/peaks_in_half_of_samples_and_lowqval/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames.filtqValHalfsamps.annotated.csv",row.names = F)
