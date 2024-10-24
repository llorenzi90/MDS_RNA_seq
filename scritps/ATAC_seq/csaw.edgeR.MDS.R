wdir=("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/results/csaw/")
setwd(wdir)
window_data <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/win.data.MDS.RDS")
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

#function to write result tables
write.res.edgeR <- function(x,text=deparse(substitute(x))){
  tm <- as.data.frame(x)
  tw <- cbind(peak_ID=rownames(tm),tm)
  write.csv(tw,paste0(text,".csv"),quote = F,row.names = F)
  
}


#filter windows based on global background
fc <- 5
#binned <- readRDS("binned_counts.RDS")
#bin.size=10000L
#binned <- windowCounts(bam.paths, bin=TRUE, width=bin.size, param=param)
binned <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/binned_counts_MDS.10K.RDS") #
filter.stat <- filterWindowsGlobal(window_data, background=binned)#this crashes locally, run in hpc
filtered.data <- window_data[filter.stat$filter > log2(fc),]

filtered.data <- readRDS("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/csaw_rds/filtered.windows.globalBackground.RDS")
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
  #filtered.data <- normFactors(binned, se.out=filtered.data)
  
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
  png(paste0("coeffvar.vs.avelogCPM.",norm_method,".png"))
  plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
       ylim=c(0, 1), xlab=expression("Ave. Log2 CPM"),
       ylab=("Biological coefficient of variation"))
  dev.off()
  png(paste0("QLDisp.",norm_method,".png"))
  plotQLDisp(fit)
  dev.off()
  
  
  fit120vs0 <- glmQLFTest(fit, contrast=c(0,0,0,0, 1))
  summary(decideTests(fit120vs0))
  results120vs0 <- topTags(fit120vs0,n = Inf)
  
  fit15vs0 <- glmQLFTest(fit, contrast=c(0,0,0,1,0))
  summary(decideTests(fit15vs0))
  results15vs0 <- topTags(fit15vs0,n = Inf)
  
  fit120vs15 <- glmQLFTest(fit, contrast=c(0,0,0,-1, 1))
  summary(decideTests(fit120vs15))
  results120vs15 <- topTags(fit120vs15,n = Inf)
  
  write.res.edgeR(results120vs0)
  write.res.edgeR(results15vs0)
  write.res.edgeR(results120vs15)
  
  sink()
}

#write filtered peaks to show in IGV
options(scipen = 999) 

ranges <- unlist(strsplit(as.character(ranges(rowRanges(filtered.data))),split = "-"))
ranges <- cbind(ranges[seq(1,length(ranges),by=2)],
                ranges[seq(2,length(ranges),by=2)])
bed_filtered_peaks <- cbind(as.character(seqnames(rowRanges(filtered.data))),
                            ranges)
bed_filtered_peaks <- as.data.frame(bed_filtered_peaks)
bed_filtered_peaks$V2 <- as.integer(bed_filtered_peaks$V2) -1 
setwd(wdir)

write.table(bed_filtered_peaks,"filtered_windows.bed",sep = "\t",quote = F,row.names = F,col.names = F)
plotMD(fit120vs0)
