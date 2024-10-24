## ---------------------------
##
##
## Purpose of script: generate scripts to run homer for Cohesin only peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-19
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: 1) homer for ~450 peaks that are only found in two or more Cohesin patients
##        background: A)peaks that are found in more than 1 sample but not cohesin
##                    B)auto, let homer take random genomic regions
##        2) idea: as a comparison do motif finding for peaks that are 
##        found in 6 random samples
##        3) Also take the differential peaks between Cohesin vs non-Cohesin patients
##           A) up peaks in Cohesin
##           B) down peaks in Cohesin
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## ---------------------------

mut_data <- fread("~/Nextcloud/MDS/mut_data_all_patients.txt")
peaks=read.table("~/Nextcloud/MDS/ATAC-seq/data/peak_files/macs2_MDS_JAN2022/merged_peaks/filtered_BlackList/merged_peaks.MDS_JAN2022.1.q_0.05_peaks.sorted.narrowPeak.bed.filtBL.withcolnames",header = T)
peaks$peakID=paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)


Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]

#overview n samples per peak
total_peaks=nrow(peaks)
summary(peaks$nsamples)
#On averge each peak is called in 4.5 samples
table(peaks$nsamples)
table(peaks$nsamples)/total_peaks*100
#63% of the peaks were called in only one sample
#only 18 peaks are called in all samples
table(peaks$nsamples>=42)
table(peaks$nsamples>=42)/total_peaks*100
#  12048 peaks (2.7%) are called in at least half of the samples

#Convert list of comma separated samples in list of vectors and retain
samps_as_vect=sapply(peaks$samples,function(x)strsplit(x,split = ","))
#remove "AT-" from old samples and remove ".1" ".2"
samps_as_vect <- sapply(samps_as_vect, function(x){t=gsub("\\.[1-2]","",x)
t=gsub("AT-","",t)})

#how many peaks are present in at least 1 Cohesin patient?
is_in_Cohesin_patient=sapply(samps_as_vect, function(x) any(x%in%Cohesin_patients))
table(is_in_Cohesin_patient)
table(is_in_Cohesin_patient)/total_peaks*100
#65789 peaks (15%)
peaks_in_Cohesin=peaks[is_in_Cohesin_patient,]

table(peaks_in_Cohesin$nsamples)
summary(peaks_in_Cohesin$nsamples)

#1) from these, peaks that are found ONLY in Cohesin patients
is_ONLY_in_Cohesin_patient=sapply(samps_as_vect[is_in_Cohesin_patient],function(x)all(x%in%Cohesin_patients))
table(is_ONLY_in_Cohesin_patient)
table(is_ONLY_in_Cohesin_patient)/total_peaks*100
#8578 peaks (2%) are exclusively found in cohesin patients
peaks_ONLY_in_Cohesin=peaks_in_Cohesin[is_ONLY_in_Cohesin_patient,]
table(peaks_ONLY_in_Cohesin$nsamples)
summary(peaks_ONLY_in_Cohesin$nsamples)
# 19 peaks are found in 3 Cohesin samples only
#95% of the cohesin only peaks (8121) are found in only one sample (60% of peaks are found in only one sample, so this is not very informative)
#438 are found in 2 cohesin samples
#this is probably not statistically significant... 
#What could be done is to randomly select 6 other patients and see what numbers do we get 


#Peaks that are not found in Cohesin
peaks_NOT_in_Cohesin=peaks[!is_in_Cohesin_patient,]
table(peaks_NOT_in_Cohesin$nsamples)
table(peaks_NOT_in_Cohesin$nsamples)/nrow(peaks_NOT_in_Cohesin)*100
summary(peaks_NOT_in_Cohesin$nsamples)

target_peaks <- peaks_ONLY_in_Cohesin[peaks_ONLY_in_Cohesin$nsamples>1,]
background_peaks <- peaks_NOT_in_Cohesin[peaks_NOT_in_Cohesin$nsamples>1,]

#write data and generate script

setwd("~/Nextcloud/MDS/ATAC-seq/analyses/")
dir.create("HOMER")
setwd("HOMER/")
dir.create("in_data")

# BED files should have at minimum 6 columns (separated by TABs, additional columns will be ignored)
# Column1: chromosome
# Column2: starting position
# Column3: ending position
# Column4: Unique Peak ID
# Column5: not used
# Column6: Strand (+/- or 0/1, where 0="+", 1="-")
colnames(target_peaks)
target_out=target_peaks[,c(1,2,3,7,5,4)]
target_out$nsamples <- "."

colnames(background_peaks)
bg_out=background_peaks[,c(1,2,3,7,5,4)]
bg_out$nsamples <- "."

write.table(target_out,"in_data/Cohesin_only_peaks.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(bg_out,"in_data/background_nonCohesin_peaks.bed",row.names = F,col.names = F,quote = F,sep = "\t")

#script:
#/Users/llorenzi/software/homer/bin/findMotifsGenome.pl in_data/Cohesin_only_peaks.bed hg38 homer_out_CohesinA -size given -bg in_data/background_nonCohesin_peaks.bed


####2) idea: as a comparison do motif finding for peaks that are 
##        found in 6 random samples

all_patients=mut_data$Sample
all_patients <- all_patients[!is.na(all_patients)]
random6sampls=sample(all_patients[!all_patients%in%Cohesin_patients],6)

#peaks found in 6 random samples
is_in_6random_patient=sapply(samps_as_vect, function(x) any(x%in%random6sampls))
table(is_in_6random_patient)
table(is_in_Cohesin_patient)

# ONLY in 6 random patients
is_ONLY_in_6random_patient=sapply(samps_as_vect[is_in_6random_patient],function(x)all(x%in%random6sampls))
table(is_ONLY_in_6random_patient)
table(is_ONLY_in_Cohesin_patient)
peaks_ONLY_in_6random_patients=peaks[is_in_6random_patient,][is_ONLY_in_6random_patient,]
summary(peaks_ONLY_in_6random_patients$nsamples)
table(peaks_ONLY_in_6random_patients$nsamples)
table(peaks_ONLY_in_Cohesin$nsamples)
peaks_NOT_in_6random_patients=peaks[!is_in_6random_patient,]

target_peaks <- peaks_ONLY_in_6random_patients[peaks_ONLY_in_6random_patients$nsamples>1,]
background_peaks <- peaks_NOT_in_6random_patients[peaks_NOT_in_6random_patients$nsamples>1,]

colnames(target_peaks)
target_out=target_peaks[,c(1,2,3,7,5,4)]
target_out$nsamples <- "."

colnames(background_peaks)
bg_out=background_peaks[,c(1,2,3,7,5,4)]
bg_out$nsamples <- "."

write.table(target_out,"in_data/6random_only_peaks.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(bg_out,"in_data/background_non6random_peaks.bed",row.names = F,col.names = F,quote = F,sep = "\t")

#/Users/llorenzi/software/homer/bin/findMotifsGenome.pl in_data/6random_only_peaks.bed hg38 homer_out_6randomtest -size given -bg in_data/background_non6random_peaks.bed


#3) differential peaks between Cohesin vs non-Cohesin patients

deseqres=read.csv("~/Nextcloud/MDS/ATAC-seq/analyses/DESeq2/filtered_peaks/peaks_in_half_of_samples_and_lowqval/all_samples/cohesin_vs_nocohesin/all_samples_cohesin_vs_nocohesin_DESeq2_results.csv")

up_peaks <- deseqres$X[deseqres$padj<0.05&deseqres$log2FoldChange>0]
down_peaks <- deseqres$X[deseqres$padj<0.05&deseqres$log2FoldChange<0]

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

write.table(target_out_UP,"in_data/879UP_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(target_out_DOWN,"in_data/430DOWN_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(bg_out,"in_data/bakground_30054nondiff_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")

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

great_background=peaks[peaks$peakID%in%deseqres$X,]
bg_out=great_background[,c(1,2,3,7,6,4)]
bg_out$nsamples <- "."
bg_out$meanscore <- round(bg_out$meanscore)

getwd()
setwd("~/Nextcloud/MDS/ATAC-seq/analyses/")
dir.create("GREAT")
setwd("GREAT/")
dir.create("in_data")
head(bg_out)
write.table(bg_out,"in_data/bakground_alltested_peaks_Cohesin.GREAT.bed",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(target_out_UP_GREAT,"in_data/879UP_peaks_Cohesin.bed",row.names = F,col.names = F,quote = F,sep = "\t")
