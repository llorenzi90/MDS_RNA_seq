wdir=("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/results/csaw/")
setwd(wdir)
window_data <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/win.data.MDS.RDS")
library(csaw)
library(edgeR)
library("ExploreModelMatrix")
#generate the design matrix:

bam.paths=window_data$bam.files
samples=gsub("(/scratch/llorenzi/MDS/ATAC-seq/AT-MDS[0-9]*.[12]/AT-)(.*)(.primary_chr.quality.concordant.sorted.markedDups.bam)",
             "\\2",
             bam.paths)

LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1])
coldata=data.frame(samples,patient,LPS)

coldata$LPS <- as.factor(coldata$LPS)

design <- model.matrix(~ coldata$LPS+ coldata$patient )
colnames(design) <- gsub("coldata\\$","",colnames(design))
design

filter.stat <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/windows.filter.stat.globalBackground.RDS")

hist(filter.stat$filter, xlab="Log-fold change from global background", 
     breaks=100, main="", col="grey80", xlim=c(0, 5))
abline(v=log2(3), col="green", lwd=2)
abline(v=log2(5), col="red", lwd=2)
fc=5
table(filter.stat$filter > log2(fc))
fc=3
fil <- filter.stat$filter > log2(fc)
table(filter.stat$filter > log2(fc))
table(filter.stat$filter > log2(fc))/length(filter.stat$filter)*100
fc
filtered.data.fc3 <- window_data[fil,]

options(scipen = 999) 

ranges <- unlist(strsplit(as.character(ranges(rowRanges(filtered.data.fc3))),split = "-"))
ranges <- cbind(ranges[seq(1,length(ranges),by=2)],
                ranges[seq(2,length(ranges),by=2)])
bed_filtered_peaks <- cbind(as.character(seqnames(rowRanges(filtered.data.fc3))),
                            ranges)
bed_filtered_peaks <- as.data.frame(bed_filtered_peaks)
bed_filtered_peaks$V2 <- as.integer(bed_filtered_peaks$V2) -1 

#write.table(bed_filtered_peaks,"../../data/filtered_windows_log2FC_3.bed",sep = "\t",quote = F,row.names = F,col.names = F)

fc=2
fil <- filter.stat$filter > log2(fc)
table(fil)
table(fil)/length(fil)*100


