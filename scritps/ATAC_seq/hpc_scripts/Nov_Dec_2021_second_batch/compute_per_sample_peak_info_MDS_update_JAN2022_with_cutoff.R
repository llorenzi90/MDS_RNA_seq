## ---------------------------
##
##
## Purpose of script: check per sample peak calling and peak coverage  
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-01-17
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: Modified on 21 Jan 2022 afeter fileterin peaks
## based on cutoff
##   
##
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)
cutoff=9
datadir="/scratch/llorenzi/macs2_MDS_JAN2022/"
outdir=paste0(datadir,"per_sample_peak_info_-qval_",cutoff)
intersect_dir=paste0(datadir,"per_sample_coverage_intersect_-log10qval.",cutoff)

#common files
sample_info <- read.csv("/scratch/llorenzi/macs2_MDS_JAN2022/MDS_summary_QC_after_alignment.JAN2022.csv")
conspeaks_file=paste0("/scratch/llorenzi/macs2_MDS_JAN2022/merged_peaks/merged_peaks.MDS_JAN2022.1.-log10qval.",cutoff,".sorted.narrowPeak")
coverage_conspeaks_file=paste0("/scratch/llorenzi/macs2_MDS_JAN2022/count_files/cutoff_-logqval_",cutoff,"/merged_macs2peaksNoDups_MDS_all_samples.-log10_qval.",cutoff,".counts")
conspeaks=read.table(conspeaks_file)
colnames(conspeaks) <- c("seqnames",	"start",	"end",	"nsamples","samples",	"meanscore")
coverage_conspeaks <- read.table(coverage_conspeaks_file)
#sample files
sample=commandArgs(trailingOnly=TRUE)[1]
sampleID=gsub("AT-","",sample)
sample_peaks_file <- paste0(datadir,"filtered_peaks_cutoff_-log10qval.",cutoff,"/",sample,".nodups.filtered_-log10qval.",cutoff)
intersect_conspeaks_file=paste0(intersect_dir,"/",sample,".intersect_conspeaks.-log10qval.",cutoff,".bed")
coverage_ownpeaks_file=paste0(datadir,"counts_own_peaks/",sample,".ownPeaks.counts")

sample_peaks <- read.table(sample_peaks_file)
sample_peaks$id <- paste0(sample_peaks$V1,":",sample_peaks$V2,"-",sample_peaks$V3) 
intersect_conspeaks <- read.table(intersect_conspeaks_file)
coverage_ownpeaks <- read.table(coverage_ownpeaks_file)
intersect_conspeaks$id <- paste0(intersect_conspeaks$V1,":",intersect_conspeaks$V2,"-",intersect_conspeaks$V3)
coverage_conspeaks$id <- rownames(coverage_conspeaks)
coverage_ownpeaks$id <- rownames(coverage_ownpeaks)

#filter coverage file
coverage_ownpeaks <- coverage_ownpeaks[coverage_ownpeaks$id%in%sample_peaks$id,]

altfile=list.files(intersect_dir,pattern = paste0("vs_",sample),full.names = T)
if (length(altfile)!=0) {
  if(file.exists(altfile)){
    intersect_alt_treat <- read.table(altfile)
    
  }else intersect_alt_treat=NULL
  
}else intersect_alt_treat=NULL

#
if(!is.null(intersect_alt_treat)){
  intersect_alt_treat$id <- paste0(intersect_alt_treat$V1,":",intersect_alt_treat$V2,"-",intersect_alt_treat$V3)
  
}
setwd(outdir)

## ---------------------------
#header for info to retrieve:
colns <- c(#some numbers on total sequnced reads and reads mapped to peaks
  "Total_read_pairs_no_dups",
  "FRiP_own",
  "FRiP_cons",
  
  #some numbers on coverage of own peaks
  "N_own_peaks",
  "mean_peak_length",
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
  "N_cons_peaks_with_overlap", #exclude unique peaks from here
  "mean_fraction_overlap", #from those calculate on average what fraction of the consensus peak is covered by this sample's peaks
  #"N_cons_peaks_min_coverage", # >= X TPM (calculate TPM considering only covered length for each peak)
  #"mean_length_covered_cons_peaks_min_coverage",
  
  #some numbers on coverage of unique peaks (are these mainly noise?)
  "N_unique_peaks",
  "Total_counts_unique_peaks",
  
  
  #some numbers comparing the two treatments for same patient (for this read coverage of alt treatment)
  "N_own_peaks_intersect_alt_treat", #(check unique peaks that overlap )
  "fraction_own_peaks_intersect_alt_treat"
  
  
)


out_info <- rep(NA,length(colns))
names(out_info) <- colns



#"Total_reads_no_dups"
sample_read_pairs_nodups <- sum(sample_info[sample_info$samples==sampleID,c("good_uniq_nonpeaks","good_uniq_in_peaks")])/2
sample_read_pairs_nodups
out_info["Total_read_pairs_no_dups"] <- sample_read_pairs_nodups

#FRiP_own
FRiP_own <- sum(coverage_ownpeaks[,1])/sample_read_pairs_nodups
out_info["FRiP_own"] <- FRiP_own

#FRiP_cons

summary(coverage_conspeaks[,sampleID])

FRiP_cons <- sum(coverage_conspeaks[,sampleID])/sample_read_pairs_nodups
out_info["FRiP_cons"] <- FRiP_cons


#N_own_peaks
N_own_peaks <- nrow(coverage_ownpeaks)
out_info["N_own_peaks"] <- N_own_peaks

#mean peak length
peak_lengths=sapply(coverage_ownpeaks$id,function(x){
  sp1=strsplit(x,split = "-")
  sp2=as.numeric(strsplit(sp1[[1]][1],split = ":")[[1]][2])
  return(as.numeric(sp1[[1]][2]) - sp2 +1)})

out_info["mean_peak_length"] <- round(mean(peak_lengths),2)

#some peaks are longer than others. To compare the peak coverage between
#peaks I need to normalize for differences in peak length
#I will take a similar approach to the calculation of TPM for genes

#1) calculate reads per kilobase
rpkb <- coverage_ownpeaks[,1]/(peak_lengths/1000)
#2) Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
pmsf <- sum(rpkb)/1000000
#3) Divide the RPK values by the “per million” scaling factor.
#I will call it fragments per million
FPM <- rpkb/pmsf


coverage_ownpeaks$FPM <- FPM

pdf(paste0(sample,".densityFPM.pdf"))
par(las=1)
plot(density(log2(coverage_ownpeaks$FPM)),main=paste0(sample, " log2(FPM) density"))
dev.off()



#cumulative FPM :
#sort peaks by abundance and calculate sum of FPM covered by those peaks
cumulative_FPM <- sapply(seq_along(coverage_ownpeaks$FPM)[-1], function(x){
  sum(sort(coverage_ownpeaks$FPM,decreasing = T)[1:x])
})

pdf(paste0(sample,".cumulativeFPM.pdf"))
par(las=1)
plot(seq_along(cumulative_FPM)/length(cumulative_FPM),cumulative_FPM/1000000,
     type="s",
     xlab="fraction top peaks",
     ylab= "fraction total FPM")
dev.off()


mean_FPM=mean(coverage_ownpeaks$FPM)

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
  out_info[paste0("fraction_peaks_over_",x,"_FPM")] <- round(sum(coverage_ownpeaks$FPM>=x)/N_own_peaks,2)
  
}

pdf(paste0(sample,".ecdf_FPM.pdf"))
plot(ecdf(log10(FPM + 0.001)))
dev.off()


#Coverage of consensus peaks
#"N_cons_peaks_with_overlap", #exclude unique peaks from here
#"mean_fraction_overlap", #from those calculate on average what fraction of the consensus peak is covered by this sample's peaks
#"N_cons_peaks_min_coverage", # >= X TPM (calculate TPM considering only covered length for each peak)
#"mean_length_covered_cons_peaks_min_coverage",
intersect_conspeaks$cons_peak_length=intersect_conspeaks$V3 - intersect_conspeaks$V2 +1
intersect_conspeaks$fraction_peak_covered=intersect_conspeaks$V17/intersect_conspeaks$cons_peak_length
intersect_conspeaks_nonunique=intersect_conspeaks[intersect_conspeaks$V5!=sample,]
N_cons_peaks_with_overlap <- length(unique(intersect_conspeaks_nonunique$id))
out_info["N_cons_peaks_with_overlap"] <- N_cons_peaks_with_overlap
fraction_covered=aggregate(intersect_conspeaks_nonunique$fraction_peak_covered,by=list(id=intersect_conspeaks_nonunique$id),sum)
out_info["mean_fraction_overlap"] <- round(mean(fraction_covered$x),2)


#coverage_conspeaks_min_unique <- coverage_conspeaks_min[coverage_conspeaks_min$V5==sample,]
peaks_unique_for_sample=unique(intersect_conspeaks$id[intersect_conspeaks$V5==sample])
out_info["N_unique_peaks"] <- length(peaks_unique_for_sample)
out_info["Total_counts_unique_peaks"] <- sum(coverage_conspeaks[coverage_conspeaks$id%in%peaks_unique_for_sample,sampleID])


#N_own_peaks_intersect_alt_treat 
out_info["N_own_peaks_intersect_alt_treat"] <- length(unique(intersect_alt_treat$id))
out_info["fraction_own_peaks_intersect_alt_treat"] <- length(unique(intersect_alt_treat$id))/N_own_peaks
outdf <- data.frame(info=names(out_info),sample=out_info)
colnames(outdf)[2] <- sample
write.csv(outdf,paste0(sample,".sample_peaks_info_-log10qval.",cutoff,".csv"),row.names = F)



#}

