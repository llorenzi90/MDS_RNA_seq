#mywdir="~/Work_local/MDS/"
mywdir="~/Nextcloud/MDS/ATAC-seq/analyses/"
setwd(mywdir)
data_files=list.files(mywdir,pattern = "^MDS")
all_data=c()
for (fi in data_files) {
  td=read.table(fi)
  tn=unlist(strsplit(gsub("MDS_","",gsub(".txt","",fi)),split = "_"))
  td=td$V1
  td=gsub("AT-","",td)
  td=gsub("MDS","",td)
  data_release_code=paste(tn[2:3],collapse = "_")
  td_sp=sapply(td, function(x) strsplit(x,split = "\\.|-"))
  patient=sapply(td_sp,function(x)x[1])
  LPS=sapply(td_sp,function(x)x[2])
  
  tmpda=data.frame(patient=patient,LPS=LPS,type=tn[1],data_release_code=data_release_code)
  all_data=rbind(tmpda,all_data)
  }

patients_with_data=sort(as.numeric(unique(all_data$patient)))

patients_with_ATAC_seq_data=sort(as.numeric(unique(all_data$patient[all_data$type=="ATAC-seq"])))

table(patients_with_ATAC_seq_data%in%patients_with_data)

unique(all_data$data_release_code)

date_data_realease=c("28/09/2021","08/09/2021",
                     "14/04/2021","08/11/2021")
names(date_data_realease)=unique(all_data$data_release_code)

all_data$date_data_realease=date_data_realease[match(all_data$data_release_code,names(date_data_realease))]

table(all_data$type)

#dups=read.table("/home/llorenzi/Work_local/multiqc_all/multiqc_data/multiqc_picard_dups.txt")
dups=read.table("QC/QC/multiqc_all/multiqc_data/multiqc_picard_dups.txt")

dups$PERCENT_DUPLICATION=gsub(",",".",dups$PERCENT_DUPLICATION)
#bowtie2=read.table("/home/llorenzi/Work_local/multiqc_all/multiqc_data/multiqc_bowtie2.txt",header = T)
bowtie2=read.table("QC/QC/multiqc_all/multiqc_data/multiqc_bowtie2.txt",header = T)

# samtools_raw=read.table("/home/llorenzi/Work_local/MDS/analyses/QC/sample_stats/stats_files/raw_markedDups/multiqc_data/multiqc_samtools_stats.txt",header = T)
# samtools_primary_chr=read.table("/home/llorenzi/Work_local/MDS/analyses/QC/sample_stats/stats_files/primary_chr_markedDups/multiqc_data/multiqc_samtools_stats.txt",header = T)
# samtools_primary_chr_proper_pairs_mq2=read.table("/home/llorenzi/Work_local/MDS/analyses/QC/sample_stats/stats_files/markedDups/multiqc_data/multiqc_samtools_stats.txt",header = T)

samtools_raw=read.table("QC/QC/sample_stats/stats_files/raw_markedDups/multiqc_data/multiqc_samtools_stats.txt",header = T)
samtools_primary_chr=read.table("QC/QC/sample_stats/stats_files/primary_chr_markedDups/multiqc_data/multiqc_samtools_stats.txt",header = T)
samtools_primary_chr_proper_pairs_mq2=read.table("QC/QC/sample_stats/stats_files/markedDups/multiqc_data/multiqc_samtools_stats.txt",header = T)

#%mapped
samtools_raw$reads_mapped/samtools_raw$raw_total_sequences*100
stperc=round(samtools_raw$reads_mapped/samtools_raw$raw_total_sequences*100,2)

btperc=bowtie2$overall_alignment_rate
table(stperc==btperc)

#% reads that map to mitochondria from the total reads
100*(samtools_raw$raw_total_sequences - 
       samtools_raw$reads_unmapped - 
       samtools_primary_chr$raw_total_sequences)/samtools_raw$raw_total_sequences

#% reads that map to mitochondria from the mapped reads
100*(samtools_raw$raw_total_sequences - 
       samtools_raw$reads_unmapped - 
       samtools_primary_chr$raw_total_sequences)/samtools_raw$reads_mapped

#% raw reads that are duplicated
samtools_raw$reads_duplicated/samtools_raw$raw_total_sequences*100
samtools_raw$reads_duplicated_percent
#% primary chr reads that are duplicated
samtools_primary_chr$reads_duplicated/samtools_primary_chr$raw_total_sequences*100
samtools_primary_chr$reads_duplicated_percent
#% primary chr and good quality that are duplicated
samtools_primary_chr_proper_pairs_mq2$reads_duplicated/samtools_primary_chr_proper_pairs_mq2$raw_total_sequences*100
samtools_primary_chr_proper_pairs_mq2$reads_duplicated_percent

#
raw_reads=samtools_raw$raw_total_sequences
unmapped_reads=samtools_raw$reads_unmapped
MT_reads=(samtools_raw$raw_total_sequences - 
            samtools_raw$reads_unmapped - 
            samtools_primary_chr$raw_total_sequences)
primary_chr_proper_mq2=samtools_primary_chr_proper_pairs_mq2$raw_total_sequences
low_qual=raw_reads- unmapped_reads - MT_reads - primary_chr_proper_mq2
duplicated_primary_chr_proper_mq2=samtools_primary_chr_proper_pairs_mq2$reads_duplicated
unique_primary_chr_proper_mq2=primary_chr_proper_mq2 - duplicated_primary_chr_proper_mq2
unique_primary_chr_proper_mq2/raw_reads*100

#mapped to peaks
#countdata=read.table("~/Work_local/MDS/merged_macs2peaksNoDups_MDS_all_samples.counts")
countdata=read.table("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/merged_macs2peaksNoDups_MDS_all_samples.counts")
reads_mapped_to_peaks=colSums(countdata)
rm(countdata)
samples=gsub("AT-","",samtools_raw$Sample)
table(names(reads_mapped_to_peaks)%in%samples)
reads_mapped_to_peaks=reads_mapped_to_peaks[match(samples,names(reads_mapped_to_peaks))]
reads_mapped_to_peaks/raw_reads*100


unique_primary_chr_proper_mq2_non_peaks=unique_primary_chr_proper_mq2 - reads_mapped_to_peaks

#gather data of interest:
rawdf=data.frame(samples=samples,
                 unmapped=unmapped_reads,
                 low_qual=low_qual,
                 MT_reads=MT_reads,
                 good_duplicated=duplicated_primary_chr_proper_mq2,
                 good_uniq_nonpeaks=unique_primary_chr_proper_mq2_non_peaks,
                 good_uniq_in_peaks=reads_mapped_to_peaks)

library(tidyverse)
library(ggplot2)
transformedDF=gather(rawdf,key = "read_type",value = "number(M)",2:7)
transformedDF$samples=gsub("MDS","",transformedDF$samples)
transformedDF$`number(M)`=transformedDF$`number(M)`/1000000
source("~/Rfunct/cbPalette.R")
deflevs=levels(as.factor(transformedDF$read_type))
transformedDF$read_type=factor(transformedDF$read_type,levels = deflevs[c(2,3,1,5,4,6)])
g=ggplot(transformedDF,aes(x=samples,fill=read_type,y=`number(M)`)) +
  geom_bar(position = "stack",stat = "identity") +scale_fill_manual(values = cbPalette[c(1,4,7,8,3,6)])

p=g + theme(axis.text.x = element_text(    size=4, angle=45,hjust = 1))

pdf("QC/all_samples_read_stats.pdf")
print(p)
dev.off()
# general_stats=fread("/home/llorenzi/Work_local/MDS/analyses/QC/multiqc_all/multiqc_data/multiqc_general_stats.txt",header = T)



#write summary sequenced data
rawdf=data.frame(samples=samples,
                 total_read_pairs=bowtie2$total_reads,
                 overall_alignment_percent=bowtie2$overall_alignment_rate,
                 unmapped=unmapped_reads,
                 low_qual=low_qual,
                 MT_reads=MT_reads,
                 good_duplicated=duplicated_primary_chr_proper_mq2,
                 good_uniq_nonpeaks=unique_primary_chr_proper_mq2_non_peaks,
                 good_uniq_in_peaks=reads_mapped_to_peaks)
write.csv(rawdf,"~/Nextcloud/MDS/ATAC-seq/analyses/QC/MDS_summary_QC_after_alignment.JAN2022.csv",row.names = F)
