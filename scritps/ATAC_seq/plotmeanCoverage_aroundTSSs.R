library(rtracklayer
        )
# gtf=readGFF("/home/llorenzi/inflammatory_genes_gencode.v38.annotation.genes.gtf")
# table(gtf$type
#       )
# list_inf_genes <- read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/list_Hallmark_inflammatory_response.txt")
# table(gtf$gene_name%in%list_inf_genes$V1)

# tss.inf_genes=gtf[gtf$gene_name%in%list_inf_genes$V1
#                     ,c("seqid","start")]
# 
#tss.inf_genes=read.table("~/Nextcloud/references/TSS_inflammatory_genes_hg38.txt",header = T)
tss.inf_genes=read.table("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_inflammatory_genes_hg38.txt",header = T)

#generate granges object with -3000 - +3000 bps around the 200 TSSs

seqs <- lapply(tss.inf_genes$seqid,function(x)rep(x,6000))
starts <- lapply(tss.inf_genes$start, function(x)return(c((x - 2999):x,(x+1):(x + 3000))))

length(starts)

extended_granges <- GRanges(unlist(seqs),IRanges(start=unlist(starts),width = 1))

#for each sample, import the regions of the bigwig files matching the regions around TSSs 
datadir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"

bwfiles=list.files(datadir,full.names = T,pattern = "bam.CPM.bw")

basename(bwfiles)
mut_data=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/mut_data_all_patients.txt",sep = "\t",header = T)

samples=gsub(".sorted.markedDups.primary_chr.proper_pairs.minq2.bam.CPM.bw","",basename(bwfiles))
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=gsub("AT-","",sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1]))
Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]
Cohesin=patient%in%Cohesin_patients


list_meanCov=list()
i=0
for (bw in bwfiles) {
  i=i+1
  tmp_tssCPM=import(bw,which=extended_granges)
  
  #extend the GRanges object so each position around -3000 - +3000 bp of each TSS has a CPM value
  all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
  
  #convert the extended ranges into a list of length equal to the number of TSSs
  splitted_ranges=split(all_TSSposCPM, ceiling(seq_along(all_TSSposCPM)/6000))
  
  #convert the list into a dataframe
  TSSposCPMdf=as.data.frame(splitted_ranges)
  
  #Take the average CPM for each position around the TSS across the 200 genes
  mean_Cov=apply(TSSposCPMdf,1,mean)
  list_meanCov[[i]] <- mean_Cov
  
}

#calculate the maximum y value across all samples, so all samples fit in the plot
maxyval=max(sapply(list_meanCov, max))


setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/")
dir.create("meanCoverage_aroundTSSs")
setwd("meanCoverage_aroundTSSs/")
dir.create("inflammatory_genes")
setwd("inflammatory_genes/")

#Plot 1 all samples Cohesin red, non-cohesin black
pdf("all_samples_CohesinRED_noCohesinBLACK.pdf")
par(las=1)
plot(-3000:2999,list_meanCov[[1]],type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "Cohesin mut",bty = "n",fill = c("red","black"))
for (sa in 1:(length(list_meanCov) - 1)) {
  color="black"
  if(Cohesin[sa])color="red"
  points(-3000:2999,list_meanCov[[sa]],type="s",col=color)
}
dev.off()

#Plot 2 all LPS1 samples Cohesin red, non-cohesin black
#set boolean to filter mean_cov list
LPS_cond <- LPS==1

maxyval=max(sapply(list_meanCov[LPS_cond], max))

pdf("LPS1_samples_CohesinRED_noCohesinBLACK.pdf")
par(las=1)
plot(-3000:2999,list_meanCov[LPS_cond][[1]],type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "Cohesin mut",bty = "n",fill = c("red","black"))
for (sa in 1:(length(list_meanCov[LPS_cond]) - 1)) {
  color="black"
  if(Cohesin[LPS_cond][sa])color="red"
  points(-3000:2999,list_meanCov[LPS_cond][[sa]],type="s",col=color)
}
dev.off()


#Plot 3 all LPS2 samples Cohesin red, non-cohesin black
#set boolean to filter mean_cov list
LPS_cond <- LPS==2

maxyval=max(sapply(list_meanCov[LPS_cond], max))

pdf("LPS2_samples_CohesinRED_noCohesinBLACK.pdf")
par(las=1)
plot(-3000:2999,list_meanCov[LPS_cond][[1]],type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "Cohesin mut",bty = "n",fill = c("red","black"))
for (sa in 1:(length(list_meanCov[LPS_cond]) - 1)) {
  color="black"
  if(Cohesin[LPS_cond][sa])color="red"
  points(-3000:2999,list_meanCov[LPS_cond][[sa]],type="s",col=color)
}
dev.off()


#plot 4 take average of LPS and non-LPS
tmp=as.data.frame(list_meanCov[LPS==1])
mean_nonLPS_all_samples=apply(tmp, 1, mean)

tmp=as.data.frame(list_meanCov[LPS==2])
mean_LPS_all_samples=apply(tmp, 1, mean)

pdf("sample_average_LPSvsNOLPS.pdf")
maxyval=max(c(max(mean_nonLPS_all_samples),max(mean_LPS_all_samples)))
plot(-3000:2999,mean_nonLPS_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "LPS",bty = "n",fill = c("red","black"))

points(-3000:2999,mean_LPS_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS",col="red")
dev.off()


##plot 5 take average of Cohesin vs non-Cohesin
tmp=as.data.frame(list_meanCov[!Cohesin])
mean_nonCohesin_all_samples=apply(tmp, 1, mean)

tmp=as.data.frame(list_meanCov[Cohesin])
mean_Cohesin_all_samples=apply(tmp, 1, mean)

maxyval=max(c(max(mean_nonCohesin_all_samples),max(mean_Cohesin_all_samples)))

pdf("sample_average_CohesinvsNOCohesin.pdf")
plot(-3000:2999,mean_nonCohesin_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "Cohesin mut",bty = "n",fill = c("red","black"))

points(-3000:2999,mean_Cohesin_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS",col="red")
dev.off()

##plot 6 same as 5 but for LPS1 only
tmp=as.data.frame(list_meanCov[LPS==1&!Cohesin])
mean_nonCohesin_all_samples=apply(tmp, 1, mean)

tmp=as.data.frame(list_meanCov[LPS==1&Cohesin])
mean_Cohesin_all_samples=apply(tmp, 1, mean)

maxyval=max(c(max(mean_nonCohesin_all_samples),max(mean_Cohesin_all_samples)))

pdf("nonLPS_sample_average_CohesinvsNOCohesin.pdf")
plot(-3000:2999,mean_nonCohesin_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "Cohesin mut",bty = "n",fill = c("red","black"))

points(-3000:2999,mean_Cohesin_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS",col="red")
dev.off()

##plot 7 same as 5 but for LPS2 only
tmp=as.data.frame(list_meanCov[LPS==2&!Cohesin])
mean_nonCohesin_all_samples=apply(tmp, 1, mean)

tmp=as.data.frame(list_meanCov[LPS==2&Cohesin])
mean_Cohesin_all_samples=apply(tmp, 1, mean)

maxyval=max(c(max(mean_nonCohesin_all_samples),max(mean_Cohesin_all_samples)))

pdf("LPS_sample_average_CohesinvsNOCohesin.pdf")
plot(-3000:2999,mean_nonCohesin_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS")
legend("topright",legend = c("YES","NO"),title = "Cohesin mut",bty = "n",fill = c("red","black"))

points(-3000:2999,mean_Cohesin_all_samples,type="s",ylim=c(0,maxyval+0.02),ylab="mean CPM",xlab="distance from TSS",col="red")
dev.off()

