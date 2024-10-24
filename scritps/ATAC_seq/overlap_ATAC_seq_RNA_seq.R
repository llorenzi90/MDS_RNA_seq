annot_peaks <- read.csv('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/peaks_in_half_of_samples_and_lowqval/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames.filtqValHalfsamps.annotated.csv')
library(tidyverse)
#read DESeq2 results
deseq2res <- read.csv("shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/all_samples/cohesin_vs_nocohesin/all_samples_cohesin_vs_nocohesin_DESeq2_results.csv")

annot_peaks <- cbind(annot_peaks,deseq2res[match(annot_peaks$peakID,deseq2res$X),])

padjco=0.05
annot_peaks <- annot_peaks %>% mutate(DIFF=ifelse(padj<=padjco,ifelse(log2FoldChange>0,"UP","DOWN"),"NO"))
annot_peaks <- annot_peaks %>%mutate(simplified_annot=gsub('(.*)( \\()(.*)(\\))','\\1',annotation))

summary(annot_peaks %>% filter(DIFF=="DOWN") %>%select(distanceToTSS))
summary(annot_peaks %>% filter(DIFF=="UP") %>%select(distanceToTSS))
table(annot_peaks %>% filter(DIFF=="DOWN") %>%select(simplified_annot))
table(annot_peaks %>% filter(DIFF=="UP") %>%select(simplified_annot))

annot_peaks %>% filter(DIFF=="DOWN") %>% group_by(simplified_annot) %>% group_map(~ summary(.x$distanceToTSS))

annot_peaks %>% filter(DIFF=="DOWN") %>% split(.$simplified_annot) %>% map(function(x)summary(x$distanceToTSS))
annot_peaks %>% filter(DIFF=="DOWN") %>% split(.$simplified_annot) %>% map(function(x)length(x$distanceToTSS))

annot_peaks %>% filter(DIFF=="UP") %>% split(.$simplified_annot) %>% map(function(x) return(c(summary(x$distanceToTSS),length(x$distanceToTSS))))

#integrate with RNA-seq data
rnaseqres <- read.csv("shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/DESeq2/all_samples/cohesin_vs_nocohesin/all_samples_CohesinvsNoCohesinDESeq2_results.gene_names.csv")
rnaseqres <- rnaseqres %>% mutate(DIFF=ifelse(padj<=padjco&!is.na(padj),ifelse(log2FoldChange>0,"UP","DOWN"),"NO"))
table(rnaseqres$DIFF)
#venn
library(ggVennDiagram)
#install.packages("ggVennDiagram")
library(ggvenn)

table(rnaseqres$gene_name[rnaseqres$DIFF!="NO"]%in%annot_peaks$gene_name)
ggvenn(list(up_rna=rnaseqres$gene_name[rnaseqres$DIFF=="UP"],
            up_ATAC_seq=annot_peaks$gene_name[annot_peaks$DIFF=="UP"]))


rnaseqres <- rnaseqres %>% mutate(proximal_ATAC_seq_peak= ifelse(gene_name%in% annot_peaks$gene_name,"YES","NO"))
rnaseqres <- rnaseqres %>% mutate(ATAC_seq_diff=annot_peaks$DIFF[match(rnaseqres$gene_name,annot_peaks$gene_name)])
table(rnaseqres$proximal_ATAC_seq_peak,rnaseqres$DIFF)
table(rnaseqres$proximal_ATAC_seq_peak,rnaseqres$DIFF,rnaseqres$ATAC_seq_diff)

annot_peaks <- annot_peaks %>% mutate(rnaseq_diff=rnaseqres$DIFF[match(annot_peaks$gene_name,rnaseqres$gene_name)])

table(annot_peaks$DIFF,annot_peaks$rnaseq_diff)

annot_peaks %>% split(.$simplified_annot) %>% map(function(x)table(x$DIFF,x$rnaseq_diff))

#there are 13 up peaks that are distal intergenic and the closest gene is up

up_peaks_up_gene <- annot_peaks$peakID[annot_peaks$DIFF=="UP"&annot_peaks$rnaseq_diff=="UP"&annot_peaks$simplified_annot=="Distal Intergenic"]
up_peaks_up_gene <- up_peaks_up_gene[!is.na(up_peaks_up_gene)]
View(annot_peaks[annot_peaks$peakID%in%up_peaks_up_gene,])
