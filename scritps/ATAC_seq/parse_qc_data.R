setwd("~/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/QC/sample_stats/")
options(scipen = 999)
li=list.files(pattern="bowtie2.log")
total_sequenced_reads=read.table("total_sequenced_reads.txt")
summary_QC_all_samples=c()
library(ggplot2)
number_of_macs2_peaks <- read.table("MDS_number_macs2_peaks_per_sample.txt")
number_of_macs2_peaks$sample_id=sapply(number_of_macs2_peaks$V2,function(x)gsub(".primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak","",basename(x)))
cols=c("sample_id","patient","LPS",
       "total_sequenced_read_pairs",
       "read_pairs_after_trimming",
       "overall_alignment_percent",
       "properly_aligned_pairs",
       "pairs_aligned_exactly_once",
       "pairs_aligned_more_than_once",
       "single_aligned_mates",
       "percent_properly_aligned_pairs",
       "number_of_MT_reads",
       "percent_MT_raw_reads",
       "MT_proper_pairs",
       "percent_MT_proper_pairs",
       "primarychr_proper_pairs",
       "percent_primarychr_proper_pairs",
       "primarychr_MAPQ20_proper_pairs",
       "percent_primarychr_MAPQ20_proper_pairs",
       "percent_proper_pairs_pass_MAPQ",
       "percent_duplicates",
       "mean_insert_size",
       "insert_size_sd",
	"number_of_macs2_peaks"
)
for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".bowtie2.log","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi,sep = "#")
  chrM=read.table(paste0(name,".chrM_SN.txt"),sep="\t", fill = T)
  primarychr <- read.table(paste0(name,".primary_chr_SN.txt"),sep="\t", fill = T)
  primarychr.quality <- read.table(paste0(name,".primary_chr.quality_SN.txt"),sep="\t", fill = T)
  primarychr.quality.concordant <- read.table(paste0(name,".primary_chr.quality.concordant.sorted.markedDups_SN.txt"),sep="\t", fill = T)
  da_split=strsplit(trimws(da$V1), split=" +")
  da_split[[grep("overall alignment rate",da$V1)]][1]
  
  
  reads_after_trimming=as.numeric(da_split[[grep("reads; of these:",da$V1)]][1])
  non_concordant_pairs=as.numeric(da_split[[grep("aligned concordantly 0 times",da$V1)[1]]][1])
  pairs_aligned_exactly_once=as.numeric(da_split[[grep("aligned concordantly exactly 1",da$V1)[1]]][1])
  pairs_aligned_more_than_once=as.numeric(da_split[[grep("aligned concordantly >1 times",da$V1)[1]]][1])
  concordantly_aligned_pairs=as.numeric(reads_after_trimming) - as.numeric(non_concordant_pairs)
  # mates_from_unproper_pairs=as.numeric(da_split[[grep("mates make up the pairs; of these:",da$V1)]][1])
  unpaired_aligned_exactly_once=as.numeric(da_split[[grep("aligned exactly 1 time",da$V1)]][1])
  unpaired_aligned_more_than_once=as.numeric(da_split[[grep("aligned >1 times",da$V1)]][1])
  single_aligned_mates =unpaired_aligned_exactly_once+ unpaired_aligned_more_than_once
  MTreads <- chrM$V2[chrM$V1=="raw total sequences:"]
  MT_reads_percent <- round((MTreads/(reads_after_trimming*2))*100,2)
  MT_proper_pairs <- chrM$V2[chrM$V1=="reads properly paired:"]/2
  MT_proper_pairs_percent <- round(MT_proper_pairs/concordantly_aligned_pairs*100,2)
  primarychr_proper_pairs <- primarychr$V2[primarychr$V1=="reads properly paired:"]/2
  primarychr_proper_pairs_percent <- round(primarychr_proper_pairs/concordantly_aligned_pairs*100,2)
  primarychr_quality_proper_pairs <- primarychr.quality$V2[primarychr.quality$V1=="reads properly paired:"]/2
  primarychr_quality_proper_pairs_percent <- round(primarychr_quality_proper_pairs/concordantly_aligned_pairs*100,2)
  percent_proper_pairs_pass_qual <- round(primarychr_quality_proper_pairs/primarychr_proper_pairs*100,2)
  percent_dups <- round(primarychr.quality.concordant$V2[primarychr.quality.concordant$V1=="reads duplicated:"]/primarychr.quality.concordant$V2[primarychr.quality.concordant$V1=="raw total sequences:"]*100,2)
  mean_insert_size <- round(primarychr.quality.concordant$V2[primarychr.quality.concordant$V1=="insert size average:"],2)
  insert_size_sd <- round(primarychr.quality.concordant$V2[primarychr.quality.concordant$V1=="insert size standard deviation:"],2)
  samp_number_of_macs2_peaks <- number_of_macs2_peaks$V1[number_of_macs2_peaks$sample_id==name]

  vec=c(name,patient, LPS,
        as.numeric(total_sequenced_reads$V2[total_sequenced_reads$V1==name]),
        reads_after_trimming,
        gsub("%","",da_split[[grep("overall alignment rate",da$V1)]][1]),
        concordantly_aligned_pairs,
        pairs_aligned_exactly_once,
        pairs_aligned_more_than_once,
        single_aligned_mates,
        round(concordantly_aligned_pairs/reads_after_trimming*100,2),
        MTreads,
        MT_reads_percent,
        MT_proper_pairs,
        MT_proper_pairs_percent,
        primarychr_proper_pairs,
        primarychr_proper_pairs_percent,
        primarychr_quality_proper_pairs, 
        primarychr_quality_proper_pairs_percent,
        percent_proper_pairs_pass_qual,
        percent_dups,
        mean_insert_size,
        insert_size_sd,
	samp_number_of_macs2_peaks
  )

  
  names(vec) <- cols
  summary_QC_all_samples=rbind(summary_QC_all_samples,vec)
}
summary_QC_all_samples <- as.data.frame(summary_QC_all_samples)
summary_QC_all_samples$total_sequenced_read_pairs=as.numeric(summary_QC_all_samples$total_sequenced_read_pairs)

write.csv(summary_QC_all_samples,"summary_QC_after_alignment.csv",row.names = F,quote = F)
#test <- read.table("test.txt",header = T)
library(ggplot2)
colnames(summary_QC_all_samples)
#plot stats and convert fields to numeric
for (f in colnames(summary_QC_all_samples)[4:23]) {
  summary_QC_all_samples[,f] <- as.numeric(summary_QC_all_samples[,f] )
  g <- ggplot(summary_QC_all_samples,
              aes_string(x="patient",fill="LPS",y=f)) + geom_col(position = "dodge2")+
    theme_classic()+ ggtitle(f)
  print(g)
}

pdf("QC_plots.pdf")
for (f in colnames(summary_QC_all_samples)[4:23]) {
  g <- ggplot(summary_QC_all_samples,
              aes_string(x="patient",fill="LPS",y=f)) + geom_col(position = "dodge2")+
    theme_classic()+ ggtitle(f)
  print(g)
}
dev.off()


for (f in colnames(summary_QC_all_samples)[4:23]) {
  png(paste0(f,"_barplot.png"))
  g <- ggplot(summary_QC_all_samples,
              aes_string(x="patient",fill="LPS",y=f)) + geom_col(position = "dodge2")+
    theme_classic()+ ggtitle(f)
  print(g)
  dev.off()
}




#fragment sizes plots
library(ggplot2)
li=list.files(pattern=".*chrM_incl.fragment_length_count.txt")
for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".chrM_incl.fragment_length_count.txt","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi)
  vect <- rep(da$V2,da$V1)
  tm=as.data.frame(vect)
  g1 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name) +xlab("Fragment size (bp)")
  print(g1)
  g2 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name)+
    xlab("Fragment size (bp)")+
    scale_y_log10()
  print(g2)
  
}


pdf("fragment_size_distributions.incl_chrM.pdf")
for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".chrM_incl.fragment_length_count.txt","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi)
  vect <- rep(da$V2,da$V1)
  tm=as.data.frame(vect)
  g1 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name) +xlab("Fragment size (bp)")
  print(g1)
  g2 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name)+
    xlab("Fragment size (bp)")+
    scale_y_log10()
  print(g2)
  
}
dev.off()

#generate png plots
for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".chrM_incl.fragment_length_count.txt","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi)
  vect <- rep(da$V2,da$V1)
  tm=as.data.frame(vect)
  g1 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name) +xlab("Fragment size (bp)")
  png(paste0(name,"_fragment_size_dist_inclchrM.png"))
  print(g1)
  dev.off()
  
  g2 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name)+
    xlab("Fragment size (bp)")+
    scale_y_log10()
  png(paste0(name,"_fragment_size_dist_inclchrM.logscale.png"))
  
  print(g2)
  dev.off()
}

##excluding chrM:
li=list.files(pattern=".*chrM_excl.fragment_length_count.txt")
for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".chrM_excl.fragment_length_count.txt","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi)
  vect <- rep(da$V2,da$V1)
  tm=as.data.frame(vect)
  g1 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name) +xlab("Fragment size (bp)")
  print(g1)
  g2 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name)+
    xlab("Fragment size (bp)")+
    scale_y_log10()
  print(g2)
  
}

li=list.files(pattern=".*chrM_excl.fragment_length_count.txt")
pdf("fragment_size_distributions.excl_chrM.pdf")
for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".chrM_excl.fragment_length_count.txt","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi)
  vect <- rep(da$V2,da$V1)
  tm=as.data.frame(vect)
  g1 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name) +xlab("Fragment size (bp)")
  print(g1)
  g2 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name)+
    xlab("Fragment size (bp)")+
    scale_y_log10()
  print(g2)
  
}
dev.off()

#generate png plots
li=list.files(pattern=".*chrM_excl.fragment_length_count.txt")

for (i in 1:length(li)) {
  fi=li[i]
  name=gsub(".chrM_excl.fragment_length_count.txt","",fi)
  LPS=strsplit(name,split = "\\.")[[1]][2]
  patient=gsub("AT-", "",strsplit(name,split = "\\.")[[1]][1])
  da=read.table(fi)
  vect <- rep(da$V2,da$V1)
  tm=as.data.frame(vect)
  g1 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name) +xlab("Fragment size (bp)")
  png(paste0(name,"_fragment_size_dist_exclchrM.png"))
  print(g1)
  dev.off()
  
  g2 <- ggplot(tm, aes(x=vect))+
    geom_density()+
    theme_classic()+
    ggtitle(name)+
    xlab("Fragment size (bp)")+
    scale_y_log10()
  png(paste0(name,"_fragment_size_dist_exclchrM.logscale.png"))
  
  print(g2)
  dev.off()
}

