options(scipen = 999)
setwd("/Users/llorenzi/MDS/RNA-seq/analyses/hisat2_summary")
li=list.files(pattern="summary")
#total_sequenced_reads=read.table("total_sequenced_reads.txt")
summary_alignment_QC_all_samples=c()
#library(ggplot2)

for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".hisat2.summary","",fi)
  pat=strsplit(name,split = "-")[[1]][1]
  LPS=strsplit(name,split = "-")[[1]][2]

  da=read.table(fi,sep = "#")
  da_split=strsplit(trimws(da$V1), split=": ")
  trimws(da_split[[grep("overall alignment rate",da$V1,ignore.case = T)]][2])
  tmpdf <- as.data.frame(da_split)
  colnames(tmpdf) <- tmpdf[1,]
  tmpdf <- tmpdf[-1,-1]
  tmpdf <- cbind(data.frame(sample=name,patient=pat,LPS=LPS),tmpdf)
  summary_alignment_QC_all_samples <- rbind(summary_alignment_QC_all_samples,tmpdf)  
  
}
write.csv(summary_alignment_QC_all_samples,"../QC/alignment_QC_summary.csv",row.names = F)
