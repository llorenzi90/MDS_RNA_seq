setwd("~/Work_local/MDS/count_files/")

files=list.files()
counts=read.table(files[1])

colnames(counts) <- gsub("AT.","",gsub(".sorted.markedDups.primary_chr.proper_pairs.minq2.bam","",colnames(counts)))

for (fi in files[2:length(files)]) {
  c2=read.table(fi)
  colnames(c2) <- gsub("AT.","",gsub(".sorted.markedDups.primary_chr.proper_pairs.minq2.bam","",colnames(c2)))
  
  counts=cbind(counts,c2)
}

write.table(counts,"../merged_macs2peaksNoDups_MDS_all_samples.counts",sep = "\t",quote = F)
effective_lib_sizes=colSums(counts)
sort(effective_lib_sizes/1000000)

