samps <- list.files("/scratch/llorenzi/MDS/ATAC-seq",pattern="AT-")

file1 <- read.table(paste0("/scratch/llorenzi/MDS/ATAC-seq/", samps[1],"/",samps[1],".consensusPeaks.fragment.counts"))

file1$id <- paste0(file1$V1,":",file1$V2,"-",file1$V3)
all_counts <- file1[, c("id","V7")]
colnames(all_counts)[2] <- samps[1]
all_cov <- file1[,c("id","V8")]
colnames(all_cov)[2] <- samps[1]

for (samp in samps[-1]) {
  tmpfile <- read.table(paste0("/scratch/llorenzi/MDS/ATAC-seq/", samp,"/",samp,".consensusPeaks.fragment.counts"))
  tmpfile$id <- paste0(tmpfile$V1,":",tmpfile$V2,"-",tmpfile$V3)
  
  tmp_counts <- tmpfile[, c("id","V7")]
  all_counts <- cbind(all_counts,tmp_counts[,2])
  colnames(all_counts)[ncol(all_counts)] <- samp
  
  tmp_cov <- tmpfile[,c("id","V8")]
  all_cov <- cbind(all_cov,tmp_cov[,2])
  colnames(all_cov)[ncol(all_cov)] <- samp
  
}

#write.table(all_counts,"/scratch/llorenzi/MDS/ATAC-seq/merged_peaks/coverage_all_samples/all_samples.counts.tsv",quote = F,sep = "\t",row.names = F)
write.table(all_cov,"/scratch/llorenzi/MDS/ATAC-seq/merged_peaks/coverage_all_samples/all_samples.covered_peak_length.tsv",quote = F,sep = "\t",row.names = F)
