## ---------------------------
##
##
## Purpose of script: analyse consensus and per sample peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-08-12
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
require(ChIPseeker)
## ---------------------------

#' Consensus peaks were generated as follows:  
#' -First MACS2 narrow peaks were generated for each sample removing duplicates _(macs2 callpeak -t $1.primary_chr.quality.concordant.sorted.markedDups.bam -n $1.primary_chr.quality.concordant.sorted.markedDups.nodups -f BAMPE -g hs --nomodel --nolambda --keep-dup 1 --call-summits --verbose 3)_  
#' -Then columns 1-4, 9 and 6 from all files were concatenated in a single file _(cut -f1-4,6,9 $s >> /scratch/llorenzi/MDS/ATAC-seq/merged_peaks/all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak done < narrowPeak_files_paths.txt)_  
#' and sorted by coordinates *(sort -k1,1 -k2,2n all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak > all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.narrowPeak)*    
#' -Finally, bedtools merge was used to merge overlapping peaks, count the number of samples that contributed to each merged peak,  
#' get the list of samples that contributed to each merged peak and report the mean of merged peaks q-values: *(sed "s/.primary_chr.quality.concordant.sorted.markedDups.nodups_peak_[0-9a-z]\*//"  all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.narrowPeak | bedtools merge -c 4,4,6 -o count_distinct,distinct,mean -i - > merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed)* 
conspeaks=fread("~/MDS/ATAC-seq/data/macs2_peaks/merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed")  
sample_info <- read.csv("/Users/llorenzi/MDS/ATAC-seq/analyses/QC/sample_stats/summary_QC_after_alignment.csv")
DT::datatable(sample_info)

libsizes <- sort(sample_info$primarychr_MAPQ20_proper_pairs)
names(libsizes) <- sample_info$sample_id[order(sample_info$primarychr_MAPQ20_proper_pairs)]
n_total_samps=length(libsizes)

#for each peak generate list of samples in which it is found:
samples_lists <- sapply(conspeaks$V5, function(x) strsplit(x,","))

#some statistics:
#Consensus peak length:
peak_len <- conspeaks$V3 - conspeaks$V2

summary(peak_len)

par(las=1)
plot(density(peak_len), main="Consensus peaks length distribution")
plot(density(peak_len), main="Consensus peaks length distribution (Zoom-in)",xlim=c(0,1000))
plot(density(log10(peak_len)),main="Consensus peaks length distribution (log10)")

#Total number of peaks
n_total_peaks=nrow(conspeaks)
n_total_peaks
#Number of peaks found in all samples:
sum(conspeaks$V4==n_total_samps)
#What percentage of consensus peaks are seen in all samples?
#% of peaks common to all samples:
table(conspeaks$V4)/n_total_peaks*100
#less than 1% of the peaks are found in all samples
#and 60% of the peaks are found in only one sample!

#how many peaks are in X% of the samples?
n_peaks_in_Xpercent_samps=c()

frac_samples <- seq(0.1,1,by=0.1)
for(perc in frac_samples){
  n_peaks_in_Xpercent_samps=c(n_peaks_in_Xpercent_samps,sum(conspeaks$V4>=perc*n_total_samps))

}
plot(n_peaks_in_Xpercent_samps, frac_samples*100,main="Number of peaks in X% of the samples",
     ylab="%samples",xlab="N peaks")
plot(n_peaks_in_Xpercent_samps/n_total_peaks*100,frac_samples*100,
     main="% of total peaks in X% of the samples",
     ylab="%samples",xlab="%total_peaks")

#add peak annotation
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

peakfile <- "~/MDS/ATAC-seq/data/macs2_peaks/merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 3000), TxDb=txdb)
#View(as.data.frame(peakAnno))
#knitr::kable(as.data.frame(peakAnno))
#add some more info: translate geneID, add p-values
gene_translation_file <- "/Users/llorenzi/references/human/biomart_tables/Human_Ensembl_genes_104_GRCh38.p13_to_NCBI.tsv"
library(data.table)
gene_translation <- fread(gene_translation_file)
peakAnnodf <- as.data.frame(peakAnno)
table(peakAnnodf$geneId%in%gene_translation$`NCBI gene (formerly Entrezgene) ID`)
peakAnnodf$gene_name <- gene_translation$`Gene name`[match(peakAnnodf$geneId,
                                                           gene_translation$`NCBI gene (formerly Entrezgene) ID`)]

peakAnnodf$simplified_annot <- gsub("(\\.*)( \\(.*)","\\1",peakAnnodf$annotation)


pie(table(peakAnnodf$simplified_annot), main="Annotation all peaks")


#
par(mfrow=c(1,3))
for (fr in frac_samples[c(2,5,8)]) {
 pie(table(peakAnnodf$simplified_annot[peakAnnodf$V4>=fr*n_total_samps]),
     main=paste0("peaks found in \n>=",fr*100," % of samples"))
  
}

#Annotation unique peaks
par(mfrow=c(1,1))
pie(table(peakAnnodf$simplified_annot[peakAnnodf$V4==1]), main="Annotation peaks found in only 1 sample")

#per sample annotations:

# pdf("~/MDS/ATAC-seq/analyses/peak_annotation/peak_annot_per_sample1.pdf")
# par(mfrow=c(3,4))
# for (samp in sort(names(libsizes))[1:12]) {
#   cond=sapply(samples_lists, function(x)samp%in%x)
#   pie(table(peakAnnodf$simplified_annot[cond]), main=samp)
# }
# dev.off()

# pdf("~/MDS/ATAC-seq/analyses/peak_annotation/peak_annot_per_sample2.pdf")
# par(mfrow=c(3,4))
# for (samp in sort(names(libsizes))[13:24]) {
#   cond=sapply(samples_lists, function(x)samp%in%x)
#   pie(table(peakAnnodf$simplified_annot[cond]), main=samp)
# }
# dev.off()

#number of peaks vs libsize
names(libsizes)
total_peaks_per_sample <- c()

for(samp in names(libsizes)){
  total_peaks_per_sample=c(total_peaks_per_sample,sum(sapply(samples_lists,function(x)samp%in%x)))
         
}

#plot(libsizes/1000000,total_peaks_per_sample)
library(RColorBrewer)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cols_12 <- colorRampPalette(cbPalette)(12)
names(cols_12) <- unique(gsub("\\.[0-9]*","",names(libsizes)))
cols_24 <- cols_12[match(gsub("\\.[0-9]*","",names(libsizes)),
                         names(cols_12))]

LPS <- gsub("(.*)(\\.)([0-9]*)","\\3",names(libsizes))
LPS_point_type <- c("1"=15,"2"=17) 

par(mfrow=c(1,1))
plot(libsizes/1000000,total_peaks_per_sample,col=cols_24,pch=LPS_point_type[match(LPS,names(LPS_point_type))],main="N peaks vs lib size\n(color=patient)",
     xlab="M reads",ylab="N peaks")
legend("topleft",legend = c("no LPS","LPS"),pch = c(15,17))

cor(libsizes/1000000,total_peaks_per_sample)
cor(libsizes/1000000,total_peaks_per_sample,method = "spearman")

plot(libsizes[1:22]/1000000,total_peaks_per_sample[1:22],col=cols_24,
     pch=LPS_point_type[match(LPS,names(LPS_point_type))][1:22],main="N peaks vs lib size\n(color=patient) - 2 outliers excluded",
     xlab="M reads",ylab="N peaks")
legend("topright",legend = c("no LPS","LPS"),pch = c(15,17))

cor(libsizes[1:22]/1000000,total_peaks_per_sample[1:22])
cor(libsizes[1:22]/1000000,total_peaks_per_sample[1:22],method = "spearman")


#number unique peaks vs libsize
total_unique_peaks_per_sample <- c()

for(samp in names(libsizes)){
  total_unique_peaks_per_sample=c(total_unique_peaks_per_sample,
                                  sum(conspeaks$V4==1&conspeaks$V5==samp))
}
#
plot(libsizes/1000000,total_unique_peaks_per_sample,col=cols_24,
     pch=LPS_point_type[match(LPS,names(LPS_point_type))],
     main="N unique peaks vs lib size\n(color=patient)",
     xlab="M reads",ylab="N unique peaks")
legend("topleft",legend = c("no LPS","LPS"),pch = c(15,17))

cor(libsizes/1000000,total_unique_peaks_per_sample)
cor(libsizes/1000000,total_unique_peaks_per_sample,method = "spearman")

plot(libsizes[1:22]/1000000,total_unique_peaks_per_sample[1:22],col=cols_24[1:22],
     pch=LPS_point_type[match(LPS,names(LPS_point_type))][1:22],
     main="N unique peaks vs lib size\n(color=patient) - 2 outliers excluded",
     xlab="M reads",ylab="N unique peaks")
legend("topleft",legend = c("no LPS","LPS"),pch = c(15,17))

cor(libsizes[1:22]/1000000,total_unique_peaks_per_sample[1:22])
cor(libsizes[1:22]/1000000,total_unique_peaks_per_sample[1:22],method = "spearman")

#ideas:
# plot de numero de picos en todas las muestras vs quitar secuencialmente la muestra
# que tiene menos total counts (para ver si el hecho de que haya pocos picos comunes 
# a todas las muestras esta relacionado con la gran diferencia en libsize)

#The idea is: if most regions of open chromatin are expected to be seen in all samples,
#then probably the samples with low total reads miss some peaks because of low counts

samps_to_remove=c()
total_common_peaks <- c()
n_total_samps=length(libsizes)
for(samp in names(libsizes)){
  samps_to_remove <- c(samps_to_remove,samp)
  samps_to_keep <- names(libsizes)[!names(libsizes)%in%samps_to_remove]
  cond <- sapply(samples_lists,function(x)all(samps_to_keep%in%x))
  total_common_peaks <- c(total_common_peaks,sum(cond))

}



total_common_peaks <- c(sum(conspeaks$V4==n_total_samps),total_common_peaks[-length(total_common_peaks)])

#plot Number of common peaks vs lib size of the sample with lowest number of peaks
par(las=1)
plot(libsizes,total_common_peaks, xlab="min lib size",ylab="Number of common peaks",
     main="Number of common peaks\nafter removing the sample with smallest lib size")


plot(seq_along(libsizes)[-1],as.numeric(table(conspeaks$V4))[-1],
     main="Number of peaks in 2 or more samples",xlab="N samples",ylab="N peaks")

plot(seq_along(libsizes),log10(as.numeric(table(conspeaks$V4))),
     main="Number of peaks in X samples",xlab="N samples",ylab="log10(N peaks)",
     xlim=c(1,24),axes=FALSE)
axis(side=2,at=c(2,3,4,5),)
axis(side=1,at=seq_along(libsizes))
box()

