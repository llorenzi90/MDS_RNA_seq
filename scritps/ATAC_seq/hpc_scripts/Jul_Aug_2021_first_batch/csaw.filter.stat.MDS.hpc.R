wdir=("/home/llorenzi/csaw_MDS")
setwd(wdir)
window_data <- readRDS("data/win.data.MDS.RDS")
library(csaw)
library(edgeR)
#generate the design matrix:

bam.paths=window_data$bam.files
samples=gsub("(/scratch/llorenzi/MDS/ATAC-seq/AT-MDS[0-9]*.[12]/AT-)(.*)(.primary_chr.quality.concordant.sorted.markedDups.bam)",
             "\\2",
             bam.paths)

LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1])
coldata=data.frame(samples,patient,LPS)

coldata$LPS <- as.factor(coldata$LPS)

design <- model.matrix(~0+coldata$patient + coldata$LPS)
colnames(design) <- gsub("coldata\\$","",colnames(design))
design


#filter windows based on global background
#fc <- 5
#binned <- readRDS("binned_counts.RDS")
#bin.size=10000L
#binned <- windowCounts(bam.paths, bin=TRUE, width=bin.size, param=param)
binned <- readRDS("data/binned_counts_MDS.10K.RDS") #
filter.stat <- filterWindowsGlobal(window_data, background=binned)#this crashes locally, run in hpc
#filtered.data <- window_data[filter.stat$filter > log2(fc),]
#saveRDS(filtered.data,"filtered.windows.globalBackground.RDS")
saveRDS(filter.stat,"windows.filter.stat.globalBackground.RDS")
