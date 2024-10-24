#setwd("~/MDS/RNA-seq/analyses/htseq_count_data/")
setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/unstranded_htseq_count/")
library(data.table)
countfiles=list.files(pattern = ".htseq.count" )
htseq_counts_all=read.table(countfiles[1],header = F,row.names = 1)
for (fi in countfiles[-1]) {
  htseq_counts_all <- cbind(htseq_counts_all,read.table(fi,header = F,row.names = 1))
}
colnames(htseq_counts_all) <- gsub(".htseq.count","",countfiles)
non_mapped_stats <- htseq_counts_all[grep("ENSG",rownames(htseq_counts_all),invert = T),]
htseq_counts <-  htseq_counts_all[grep("ENSG",rownames(htseq_counts_all)),]
total_sequenced_reads <- apply(htseq_counts_all,2,sum)
total_gene_counts <- apply(htseq_counts,2,sum)
total_gene_counts/total_sequenced_reads*100
total_non_mapped_counts=apply(non_mapped_stats,2, sum)
all_stats=rbind("__total_non_mapped_counts"=total_non_mapped_counts,non_mapped_stats)
all_stats=rbind("__total_gene_counts"=total_gene_counts,all_stats)
all_stats=rbind("__total_sequenced_counts"=total_sequenced_reads,all_stats)
write.csv(cbind(gene_id=rownames(htseq_counts),htseq_counts),
          "second_batch_gencodev38.counts.csv",quote = F,row.names = F)
write.csv(t(all_stats),"all_stats_htseq.second_batch.csv")
all_stats["__total_gene_counts",]/all_stats["__no_feature",]
#dir.create("../htseq_count_data/count_matrix")

#merge with previous data
new_all_stats <- cbind(colnames(all_stats),t(all_stats))
new_all_stats <- as.data.frame(new_all_stats)
new_all_stats$V1 <- gsub("-",".",new_all_stats$V1)
prev_all_stats <- read.csv("all_stats_htseq.csv")
colnames(prev_all_stats) <- colnames(new_all_stats)
both_all_stats <- rbind(prev_all_stats,new_all_stats)
write.csv(both_all_stats,"all_stats_htseq.all_samples.csv")

#do the same for count matrix
prev_count_data <- read.csv("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv",header = T)
write.table(prev_count_data,"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/first_batch_gencodev38.counts.csv",quote = F,row.names = F)
colnames(prev_count_data) <- gsub("X","",colnames(prev_count_data))
colnames(htseq_counts) <- gsub("-",".",colnames(htseq_counts))

all_count_data=cbind(prev_count_data,htseq_counts)
write.csv(all_count_data,"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv",
          row.names = F,quote = F)
