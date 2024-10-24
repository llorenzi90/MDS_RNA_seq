dir="/scratch/llorenzi/MDS/ATAC-seq/per_sample_peak_info"
options(scipen = 999)
setwd(dir)
files <- list.files(pattern = "sample_peaks_info.csv")

all_info <- read.csv( files[1])


for (fi in files[-1]) {
  tmpfile <- read.csv(fi)
  all_info <- cbind(all_info,tmpfile[,2])
  colnames(all_info)[ncol(all_info)] <- colnames(tmpfile)[2]
  

}

write.csv(all_info,"all_samples_sample_peaks_info.csv",row.names = F)
