dir="/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/QC/per_sample_peak_info/per_sample_peak_info_-qval_9/"
options(scipen = 999)
setwd(dir)
files <- list.files(pattern = "sample_peaks_info_-log10qval.9.csv")

all_info <- read.csv( files[1])


for (fi in files[-1]) {
  tmpfile <- read.csv(fi)
  all_info <- cbind(all_info,tmpfile[,2])
  colnames(all_info)[ncol(all_info)] <- colnames(tmpfile)[2]
  

}

write.csv(all_info,"all_samples_sample_peaks_info.-log10qval.9.csv",row.names = F)
