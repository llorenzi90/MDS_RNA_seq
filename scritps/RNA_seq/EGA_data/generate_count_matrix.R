setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/htseq_counts/")
library(data.table)
countfiles=list.files(pattern = "htseq.count" )
htseq_counts_all=read.table(countfiles[1],row.names = 1)
for (fi in countfiles[-1]) {
  htseq_counts_all <- cbind(htseq_counts_all,read.table(fi,row.names = 1))
}

non_mapped_stats <- htseq_counts_all[grep("ENSG",rownames(htseq_counts_all),invert = T),]
htseq_counts <-  htseq_counts_all[grep("ENSG",rownames(htseq_counts_all)),]
total_sequenced_reads <- apply(htseq_counts_all,2,sum)
total_gene_counts <- apply(htseq_counts,2,sum)
total_gene_counts/total_sequenced_reads*100
total_non_mapped_counts=apply(non_mapped_stats,2, sum)
all_stats=rbind("__total_non_mapped_counts"=total_non_mapped_counts,non_mapped_stats)
all_stats=rbind("__total_gene_counts"=total_gene_counts,all_stats)
all_stats=rbind("__total_sequenced_counts"=total_sequenced_reads,all_stats)

dir.create("stats_plots")
setwd("stats_plots/")
pdf("all_stats.boxplot.pdf")
par(las=2)
boxplot(t(all_stats/1000000),ylab="# reads (M)")
dev.off()

library(tidyverse)
ggpdata <- gather(as.data.frame(t(all_stats)),key = "read_type",value = "read_number")
library(ggplot2)
pdf("all_stats.boxplot.ggplot.pdf")
ggplot(ggpdata,aes(y=read_number/1000000,x=read_type))+geom_boxplot()+
  ylab("# reads (M)")+ 
  theme_classic() +theme(axis.text.x  = element_text(angle = 45,hjust = 1))
dev.off()

all_stats_prop <- apply(as.data.frame(all_stats), 1,function(x)x/as.numeric(all_stats["__total_sequenced_counts",]))
ggpdata <- gather(as.data.frame(all_stats_prop),key = "read_type",value = "read_number")
pdf("all_stats.boxplot.ggplot.fraction.pdf")
ggplot(ggpdata,aes(y=read_number,x=read_type))+geom_boxplot()+
  ylab("fraction of sequenced reads")+ 
  theme_classic() +theme(axis.text.x  = element_text(angle = 45,hjust = 1))
dev.off()

write.csv(all_stats,"../htseq_count_stats.csv")

colnames(htseq_counts) <- gsub(".bam.unstranded.htseq.count","",countfiles)
write.csv(cbind(gene_id=rownames(htseq_counts),htseq_counts),"../all_samples_gencodev19.counts.csv",quote = F,row.names = F)


###Convert gene ids to gene names#####
library(biomaRt)

listEnsemblArchives()
ensGRCh37 <- useEnsembl(biomart = "genes",version = "GRCh37")
mart=useDataset(dataset = "hsapiens_gene_ensembl",mart = ensGRCh37)
searchAttributes(mart = mart,pattern = "gene_name")
searchAttributes(mart = mart,pattern = "symbol")
searchAttributes(mart = mart,pattern = "version")

test=getBM(attributes = "ensembl_gene_id_version",mart = mart)
table(rownames(htseq_counts)%in%test$ensembl_gene_id_version)
rownames(htseq_counts)[!rownames(htseq_counts)%in%test$ensembl_gene_id_version]

G_list <- getBM(filters= "ensembl_gene_id_version",
                attributes= c("ensembl_gene_id_version","ensembl_gene_id",
                              "external_gene_name", "hgnc_symbol"),
                values=rownames(htseq_counts),mart= mart)
gene_names=G_list$external_gene_name[match(rownames(htseq_counts),G_list$ensembl_gene_id_version)]

#########keep gene_ids for those that have no translation#######
table(is.na(gene_names))
gene_names[is.na(gene_names)] <- rownames(htseq_counts)[is.na(gene_names)]

###deduplicate duplicated gene names##########
table(duplicated(gene_names))
duplicated_gene_names <- gene_names[duplicated(gene_names)]


ord=order(duplicated_gene_names)
ndups=table(duplicated_gene_names)
table(unique(duplicated_gene_names[ord])==names(ndups))

dup_vect <- sapply(ndups, function(x)return(paste0("_",seq(1:x))))
modified_names <- duplicated_gene_names
modified_names[ord] <- paste0(duplicated_gene_names[ord],unlist(dup_vect))

gene_names[duplicated(gene_names)] <- modified_names
anyNA(gene_names)
any(duplicated(gene_names))

write.csv(cbind(gene_name=gene_names,htseq_counts),"../all_samples_gencodev19.counts.gene_name.csv",quote = F,row.names = F)
