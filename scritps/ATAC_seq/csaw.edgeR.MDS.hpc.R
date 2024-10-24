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
fc <- 5
#binned <- readRDS("binned_counts.RDS")
#bin.size=10000L
#binned <- windowCounts(bam.paths, bin=TRUE, width=bin.size, param=param)
binned <- readRDS("data/binned_counts_MDS.10K.RDS") #
filter.stat <- filterWindowsGlobal(window_data, background=binned)#this crashes locally, run in hpc
filtered.data <- window_data[filter.stat$filter > log2(fc),]
saveRDS(filtered.data,"filtered.windows.globalBackground.RDS")

#introduce a loop to try out the two alternative normalization methods
norm_method_list <- list(composition=binned,
                         efficiency=filtered.data)
for (norm_method in names(norm_method_list)) {
  setwd(wdir)
  dir.create(norm_method)
  setwd(norm_method)
  sink(paste0(norm_method,".log"))
  filtered.data <- window_data[filter.stat$filter > log2(fc),]
  filtered.data <- normFactors(norm_method_list[[norm_method]], se.out=filtered.data)
  filtered.data$norm.factors
  
  #construct the DGE list
  y <- asDGEList(filtered.data)
  rownames(y) <- paste(seqnames(rowRanges(filtered.data)),ranges(rowRanges(filtered.data)),sep = ":")
  #create the design matrix
  
  rownames(design) <- colnames(y)
  #estimate dispersion
  y <- estimateDisp(y, design)
  summary(y$trended.dispersion)
  fit <- glmQLFit(y, design, robust=TRUE)
  summary(fit$var.post)
  par(mfrow=c(1,2))
  o <- order(y$AveLogCPM)
  png(paste0("BiolCoeffVar.vs.AveLogCPM.",norm_method,".png"))
  plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
       ylim=c(0, 1), xlab=expression("Ave. Log2 CPM"),
       ylab=("Biological coefficient of variation"))
  dev.off()
  png(paste0("QLDisp.",norm_method,".png"))
  plotQLDisp(fit)
  dev.off()
  
  saveRDS(fit,paste0(norm_method,"_fit.RDS"))
  
  
  sink()
}
