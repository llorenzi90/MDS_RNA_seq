## ---------------------------
##
##
## Purpose of script: check per sample peak calling and peak coverage  
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-08-30
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

#require(tidyverse)
#require(data.table)
#sample_info <- read.csv("/Users/llorenzi/MDS/ATAC-seq/analyses/QC/sample_stats/summary_QC_after_alignment.csv")
sample_info <- read.csv("/scratch/llorenzi/MDS/ATAC-seq/summary_QC_after_alignment.csv")
#conspeaks=read.table("~/MDS/ATAC-seq/data/macs2_peaks/merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed")  
conspeaks=read.table("/scratch/llorenzi/MDS/ATAC-seq/merged_peaks/macs2_peaks/merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed")  

## ---------------------------
#header for info to retrieve:
colns <- c(#some numbers on total sequnced reads and reads mapped to peaks
           "Total_raw_reads_no_dups",
           "FRiP_own",
           "FRiP_cons",
           
           #some numbers on coverage of own peaks
           "N_own_peaks",
           "fraction_peaks_0.8", #fraction top peaks that account for 80% of TPM
           "fraction_peaks_0.6", #fraction top peaks that account for 60% of TPM
           "fraction_peaks_0.4", #fraction top peaks that account for 40% of TPM (I will calculate this also to plot cumulative distributions)
           "mean_FPM",
           "fraction_peaks_over_10_FPM", #number of peaks over X FPM (decide cutoff later)
           "fraction_peaks_over_20_FPM",
           "fraction_peaks_over_30_FPM",
           "fraction_peaks_over_40_FPM",
           "fraction_peaks_over_50_FPM",
           "fraction_peaks_over_60_FPM",
           "fraction_peaks_over_70_FPM",
           "fraction_peaks_over_80_FPM",
           "fraction_peaks_over_90_FPM",
           "fraction_peaks_over_100_FPM",
           #some stats on number of peaks and coverage over consensus peaks:
           "N_cons_peaks_min_coverage", # >= X TPM (calculate TPM considering only covered length for each peak)
           "mean_length_covered_cons_peaks_min_coverage",
           
           #some numbers on coverage of unique peaks (are these mainly noise?)
           "N_unique_peaks",
           "Total_counts_unique_peaks",
           "Max_FPM_unique_peaks",
           "Total_FPM_unique_peaks",
           
           #some numbers comparing the two treatments for same patient (for this read coverage of alt treatment)
           "N_own_peaks_intersect_alt_treat", #(check unique peaks that overlap )
           "fraction_own_peaks_intersect_alt_treat",
           "corr_cons_peaks_alt_treat",
           "corr_cons_peaks_alt_treat_third_quantile",
           #"unique_peaks_treatment", #peaks not in the alternative treatment
           #"unique_peaks_treatment_over_X_FPM",
           "number_cons_peaks_common__third_quantile", #those cons peaks that have more than X TPM in both samples 
           "fraction_cons_peaks_common__third_quantile"
           
           
           )

#sample="AT-MDS20.2"
#sample="AT-MDS35.2"
sample=commandArgs(trailingOnly = T)[1]
#datadir="~/MDS/ATAC-seq/data/bedtools_output/"
datadir=paste0("/scratch/llorenzi/MDS/ATAC-seq/",sample)
#intersect_conspeaks <- read.table(paste0(datadir,"intersect/",sample,".intersect_conspeaks.bed"))
intersect_conspeaks <- read.table(paste0(datadir,"/",sample,".intersect_conspeaks.bed"))
#intersect_alt_treat <- read.table(list.files(paste0(datadir,"intersect/"),pattern = paste0("*vs_",sample,"*"),full.names = T))
intersect_alt_treat <- read.table(list.files(datadir,pattern = paste0("*vs_",sample,"*"),full.names = T))
#coverage_conspeaks <- read.table(paste0(datadir,"coverage/",sample,".consensusPeaks.fragment.counts"))
coverage_conspeaks <- read.table(paste0(datadir,"/",sample,".consensusPeaks.fragment.counts"))
#coverage_ownpeaks <- read.table(paste0(datadir,"coverage/",sample,".narrowPeak.fragment.counts"))
coverage_ownpeaks <- read.table(paste0(datadir,"/",sample,".narrowPeak.fragment.counts"))

intersect_conspeaks$id <- paste0(intersect_conspeaks$V1,":",intersect_conspeaks$V2,"-",intersect_conspeaks$V3)
intersect_alt_treat$id <- paste0(intersect_alt_treat$V1,":",intersect_alt_treat$V2,"-",intersect_alt_treat$V3)
coverage_conspeaks$id <- paste0(coverage_conspeaks$V1,":",coverage_conspeaks$V2,"-",coverage_conspeaks$V3)
coverage_ownpeaks$id <- paste0(coverage_ownpeaks$V1,":",coverage_ownpeaks$V2,"-",coverage_ownpeaks$V3)

coverage_ownpeaks_unique <- coverage_ownpeaks[!duplicated(coverage_ownpeaks$id),]


out_info <- rep(NA,length(colns))
names(out_info) <- colns



#"Total_raw_reads_no_dups"
sample_raw_reads <- sample_info$primarychr_MAPQ20_proper_pairs[sample_info$sample_id==sample]
frac_dups <- sample_info$percent_duplicates[sample_info$sample_id==sample]/100
sample_raw_reads_nodups <- sample_raw_reads - sample_raw_reads*frac_dups

out_info["Total_raw_reads_no_dups"] <- sample_raw_reads_nodups

#FRiP_own
FRiP_own <- sum(coverage_ownpeaks_unique$V11)/sample_raw_reads_nodups
out_info["FRiP_own"] <- FRiP_own

#FRiP_cons
#first filter out those consensus peaks that have too low overlap:

summary(coverage_conspeaks$V8)

coverage_conspeaks <- coverage_conspeaks[coverage_conspeaks$V8>=min(coverage_ownpeaks$V13),]

FRiP_cons <- sum(coverage_conspeaks$V7)/sample_raw_reads_nodups
out_info["FRiP_cons"] <- FRiP_cons


#N_own_peaks
N_own_peaks <- nrow(coverage_ownpeaks_unique)
out_info["N_own_peaks"] <- N_own_peaks


#some peaks are longer than others. To compare the peak coverage between
#peaks I need to normalize for differences in peak length
#I will take a similar approach to the calculation of TPM for genes

#1) calculate reads per kilobase
rpkb <- coverage_ownpeaks_unique$V11/(coverage_ownpeaks_unique$V13/1000)
#2) Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
pmsf <- sum(rpkb)/1000000
#3) Divide the RPK values by the “per million” scaling factor.
#I will call it fragments per million
FPM <- rpkb/pmsf


coverage_ownpeaks_unique$FPM <- FPM

pdf(paste0(sample,".densityFPM.pdf"))
par(las=1)
plot(density(log2(coverage_ownpeaks_unique$FPM)),main=paste0(sample, " log2(FPM) density"))
dev.off()



#cumulative FPM :
#sort peaks by abundance and calculate sum of FPM covered by those peaks
cumulative_FPM <- sapply(seq_along(coverage_ownpeaks_unique$FPM)[-1], function(x){
  sum(sort(coverage_ownpeaks_unique$FPM,decreasing = T)[1:x])
})

pdf(paste0(sample,".cumulativeFPM.pdf"))
par(las=1)
plot(seq_along(cumulative_FPM)/length(cumulative_FPM),cumulative_FPM/1000000,
     type="s",
     xlab="fraction top peaks",
     ylab= "fraction total FPM")
dev.off()


mean_FPM=mean(coverage_ownpeaks_unique$FPM)

out_info["mean_FPM"] <- round(mean_FPM,2)

#"fraction_peaks_0.8", #N top peaks that account for 80% of FPM
fraction_peaks_0.8 <- round(which(cumulative_FPM>=800000)[1]/N_own_peaks,2)
out_info["fraction_peaks_0.8"] <- fraction_peaks_0.8
#"fraction_peaks_0.6", #N top peaks that account for 60% of FPM
fraction_peaks_0.6 <- round(which(cumulative_FPM>=600000)[1]/N_own_peaks,2)
out_info["fraction_peaks_0.6"] <- fraction_peaks_0.6
#"fraction_peaks_0.4", #N top peaks that account for 40% of FPM  
fraction_peaks_0.4 <- round(which(cumulative_FPM>=400000)[1]/N_own_peaks,2)
out_info["fraction_peaks_0.4"] <- fraction_peaks_0.4
#points(c(fraction_peaks_0.8,fraction_peaks_0.6,fraction_peaks_0.4),c(0.8,0.6,0.4),col="red")

#
FPMs<- seq(10,100,by=10)
for (x in FPMs) {
  out_info[paste0("fraction_peaks_over_",x,"_FPM")] <- round(sum(coverage_ownpeaks_unique$FPM>=x)/N_own_peaks,2)
    
}

pdf(paste0(sample,".ecdf_FPM.pdf"))
plot(ecdf(log10(FPM)))
dev.off()



#N_cons_peaks_min_coverage 
#I will use (a bit arbitrarily) 10 FPM as cutoff


#Calculate TPM for consensus peaks (over covered length):
#1) calculate reads per kilobase
rpkb <- coverage_conspeaks$V7/(coverage_conspeaks$V8/1000)
#2) Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
pmsf <- sum(rpkb,na.rm = T)/1000000
#3) Divide the RPK values by the “per million” scaling factor.
#I will call it fragments per million
FPM <- rpkb/pmsf

coverage_conspeaks$FPM <- FPM

cutoff=10
coverage_conspeaks_min <- coverage_conspeaks[!is.na(coverage_conspeaks$FPM)&(coverage_conspeaks$FPM>cutoff),]

N_cons_peaks_min_coverage <- nrow(coverage_conspeaks_min)
out_info["N_cons_peaks_min_coverage"] <- N_cons_peaks_min_coverage

out_info["mean_length_covered_cons_peaks_min_coverage"] <- round(mean(coverage_conspeaks_min$V10),2)

coverage_conspeaks_min_unique <- coverage_conspeaks_min[coverage_conspeaks_min$V5==sample,]
out_info["N_unique_peaks"] <- nrow(coverage_conspeaks_min_unique)
out_info["Total_counts_unique_peaks"] <- sum(coverage_conspeaks_min_unique$V7)
out_info["Max_FPM_unique_peaks"] <- max(coverage_conspeaks_min_unique$FPM)
out_info["Total_FPM_unique_peaks"] <- sum(coverage_conspeaks_min_unique$FPM)


#N_own_peaks_intersect_alt_treat 
out_info["N_own_peaks_intersect_alt_treat"] <- length(unique(intersect_alt_treat$id))
out_info["fraction_own_peaks_intersect_alt_treat"] <- length(unique(intersect_alt_treat$id))/N_own_peaks

treatment <- gsub("(AT-MDS[0-9]*\\.)(.)","\\2",sample)
patient <- gsub("(AT-MDS[0-9]*)(.*)","\\1",sample)

if(treatment==1){alt_treat <- paste0(patient,".2")}else(alt_treat <- paste0(patient,".1"))
  

coverage_conspeaks_alt <- read.table(paste0("/scratch/llorenzi/MDS/ATAC-seq/",
                                       alt_treat,"/",alt_treat,".consensusPeaks.fragment.counts"))
#coverage_conspeaks_alt <- read.table(paste0("~/MDS/ATAC-seq/data/bedtools_output/coverage/",alt_treat,".consensusPeaks.fragment.counts"))

coverage_conspeaks_alt$id <- paste0(coverage_conspeaks_alt$V1,":",coverage_conspeaks_alt$V2,"-",coverage_conspeaks_alt$V3)

#check peaks in common for 80% top peaks in both samples

coverage_conspeaks_alt <- coverage_conspeaks_alt[coverage_conspeaks_alt$V8>=min(coverage_ownpeaks$V13),]

#Calculate TPM for consensus peaks (over covered length):
#1) calculate reads per kilobase
rpkb <- coverage_conspeaks_alt$V7/(coverage_conspeaks_alt$V8/1000)
#2) Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
pmsf <- sum(rpkb,na.rm = T)/1000000
#3) Divide the RPK values by the “per million” scaling factor.
#I will call it fragments per million
FPM <- rpkb/pmsf

#top 80% peaks

coverage_conspeaks_alt$FPM <- FPM
#plot(density(log2(coverage_conspeaks_alt$FPM[!is.na(coverage_conspeaks_alt$FPM)])))

cutquantile <- quantile(coverage_conspeaks$FPM[!is.na(coverage_conspeaks$FPM)],0.75)
cutquantile_alt <- quantile(coverage_conspeaks_alt$FPM[!is.na(coverage_conspeaks_alt$FPM)],0.75)

coverage_conspeaks_up_75_quantile <- coverage_conspeaks[coverage_conspeaks$FPM>cutquantile&!is.na(coverage_conspeaks$FPM),]
coverage_conspeaks_alt_up_75_quantile <- coverage_conspeaks_alt[coverage_conspeaks_alt$FPM>cutquantile_alt&!is.na(coverage_conspeaks_alt$FPM),]
table(coverage_conspeaks_up_75_quantile$id%in%coverage_conspeaks_alt_up_75_quantile$id)

out_info["corr_cons_peaks_alt_treat"] <- cor(coverage_conspeaks$FPM[coverage_conspeaks$id%in%coverage_conspeaks_alt$id],
                                             coverage_conspeaks_alt$FPM[coverage_conspeaks_alt$id%in%coverage_conspeaks$id])

out_info["corr_cons_peaks_alt_treat_third_quantile"] <-  cor(coverage_conspeaks_up_75_quantile$FPM[coverage_conspeaks_up_75_quantile$id%in%coverage_conspeaks_alt_up_75_quantile$id],
                                                               coverage_conspeaks_alt_up_75_quantile$FPM[coverage_conspeaks_alt_up_75_quantile$id%in%coverage_conspeaks_up_75_quantile$id])

out_info["number_cons_peaks_common__third_quantile"] <- sum(coverage_conspeaks_up_75_quantile$id%in%coverage_conspeaks_alt_up_75_quantile$id)

out_info["fraction_cons_peaks_common__third_quantile"] <- sum(coverage_conspeaks_up_75_quantile$id%in%coverage_conspeaks_alt_up_75_quantile$id)/nrow(coverage_conspeaks_up_75_quantile)

outdf <- data.frame(info=names(out_info),sample=out_info)
colnames(outdf)[2] <- sample
write.csv(outdf,paste0(sample,".sample_peaks_info.csv"),row.names = F)

