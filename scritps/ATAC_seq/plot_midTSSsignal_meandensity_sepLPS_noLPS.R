## ---------------------------
##
##
## Purpose of script: Plot profile signal plots (average coverage)
##                    around TSS of different set of genes 
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-23
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: I use bigwig files normalized by CPM as input signal
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
#require(tidyverse)
#require(data.table)
require(rtracklayer)
require(ComplexHeatmap)
## ---------------------------

#run: Rscript plot_midTSSsignal_heatmap.R TSS_file datadir outdir sample_data_file Cohesin
#define input data and paths:
cArgs=commandArgs(trailingOnly = T)

#we need: folder with bw files
# TSS file with hg38 coordinates 
TSS_file <- cArgs[1] #example "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_files/HALLMARK_INFLAMMATORY_RESPONSE.TSSs.txt"
#TSS_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_files/HALLMARK_INFLAMMATORY_RESPONSE.TSSs.txt"
setname <- gsub("\\.[a-zA-Z0-9]*$","",gsub("\\.txt$","",basename(TSS_file)))

#coverage datadir
datadir <- cArgs[2]
#datadir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"

# outdir
outdir <- cArgs[3]
#outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/average_plots/"
if(!dir.exists(outdir)) dir.create(outdir,recursive = T)

# sample data table
sample_data_file <- cArgs[4] 
#sample_data_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"

#feature to group data: #has to match a colname in the sample data table
#group_feat <- cArgs[5]
#group_feat <- "Cohesin"

#file_type 
#filextension <- cArgs[6]
filextension <- ".bw" #currently only bigwig files supported


sort_by="CPM"
region_to_plot=c(-1000,1000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)

#Color scale parameters
range_col_values=c(-6,0)
cada=0.02
mi=range_col_values[1]
ma=range_col_values[2]



######Read and process TSS positions####
##Read TSSs coords 
TSSc <- read.table(TSS_file,header = T)
colnames(TSSc)[1:2] <- c("seqid","start") 
  
extended_granges <- GRanges(TSSc$seqid,
                            IRanges(start = TSSc$start + llimit,
                                    width = rwidth))

#Keep only peaks in canonical chrs
canonical_chrs <- paste0("chr",c(seq(1:22),"X","Y"))
print("peaks in non-canonical chromosomes")
print(table(!seqnames(extended_granges)%in%canonical_chrs))
print("removing peaks in following non-canonical chrs:")
print(unique(as.character(seqnames(extended_granges)[!seqnames(extended_granges)%in%canonical_chrs])))
extended_granges <- extended_granges[seqnames(extended_granges)%in%canonical_chrs]

###Remove possible overlapping regions 
tfo=findOverlaps(extended_granges,extended_granges)
overlaps=data.frame(q=queryHits(tfo),s=subjectHits(tfo))
overlaps <- overlaps[overlaps$q!=overlaps$s,]
#first remove the repeated rows, keep only those in which q is smaller than s 
overlapsf <- overlaps[overlaps$q<overlaps$s,]
rows_to_remove <- unique(overlapsf$s)
if(length(rows_to_remove)>0)extended_granges <- extended_granges[-rows_to_remove]


######Read and process coverage data####
##Read in coverage data and generate the variables of interest
bwfiles <- list.files(datadir,full.names = T,pattern = filextension)

sample_table <- read.csv(sample_data_file) #this have to have a column "samples"
sample_table$file <- bwfiles[unlist(lapply(sample_table$samples,function(x)grep(x,bwfiles)))]

list_all_PeaksposCPM <- list()
i=0
for (bw in sample_table$file) {
#for (bw in sample_table$file[1:10]) { #to test
  i=i+1
  tmp_tssCPM=import(bw,which=extended_granges)
  #extend the GRanges object so each position around -Xkb - +X kb of each TSS has a CPM value
  all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
  list_all_PeaksposCPM[[i]] <- all_TSSposCPM
}
print("Finished data import\nChecking if incomplete peaks")
#####Check that each feature have all positions. If incomplete, remove!####
# generate some indexes that will be useful to remove incomplete 
#peaks in case they happen:
#1) generate single-position ids in the peaks object:
posIDorigin <- paste0(as.character(rep(seqnames(extended_granges),each=rwidth)),"_",
                      unlist(apply(as.data.frame(ranges(extended_granges)),1,function(x)return(x[1]:x[2]))))
#2)generate corresponding peakIDs for each position in the original peaks:
chrs_origin <- as.character(rep(seqnames(extended_granges),each=rwidth))
ranges_origin <- rep(unlist(apply(as.data.frame(ranges(extended_granges)),
                                  1,function(x)return(paste(x[1:2],collapse = "-")))),
                     each=rwidth)
peakIDorigin <- paste0(chrs_origin,":",ranges_origin)
#3) In the same way, generate single-position ids in the coverage object
posIDCovobj <- paste0(as.character(rep(seqnames(tmp_tssCPM),width(ranges(tmp_tssCPM)))),"_",
                      unlist(apply(as.data.frame(ranges(tmp_tssCPM)),1,function(x)return(x[1]:x[2]))))
#4) Use both to assign peak of origin for each position in the coverage object
peakIDCovobj <- peakIDorigin[match(posIDCovobj,posIDorigin)]
#5) use rle to check if any incomplete peakIDs
rle_peakIDCovobj <- rle(peakIDCovobj)
if(any(rle_peakIDCovobj$lengths!=rwidth)){
  peaks_to_remove <- rle_peakIDCovobj$values[rle_peakIDCovobj$lengths!=rwidth]
  print("peaks to remove because incomplete:")
  print(peaks_to_remove)
  print("These correspond to")
  positions_to_remove <- peakIDCovobj%in%peaks_to_remove
  table(positions_to_remove)
  print(as.numeric(table(positions_to_remove)[2]))
  print("positions")
  #remove positions to remove in each sample position-coverage vector 
  list_all_PeaksposCPM <- lapply(list_all_PeaksposCPM,function(x)return(x[!positions_to_remove]))
  print("After removing incomplete peaks, there are left")
  print(lapply(list_all_PeaksposCPM, function(x)length(x)/rwidth))
  print("TSSs to plot")
  
}else{
  print("There are")
  print(lapply(list_all_PeaksposCPM, function(x)length(x)/rwidth))
  print("TSSs to plot")
}
#####

####Compute mean CPM for each position*gene for each sample:
mean_posCPM <- lapply(list_all_PeaksposCPM, function(x)aggregate(x,by=list(pos=rep(1:rwidth,length(x)/rwidth)),mean))
mean_posCPM <- as.data.frame(lapply(mean_posCPM,function(s)return(s$x)))

######################PLOTS######################
#####1) Plot density for each sample showing each group in a different color
#define color palette (color-blind friendly)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#define function to plot
plot_densities <- function(dfin,lwd=2,main="",
         xlab="distance from TSS",
         ylab="mean CPM",xvals=llimit:ulimit,
         cols=cbPalette[6:7],
         sampleTable=sample_table,
         groupCol=group_feat){
  #Define sample sets: has to be a factor of two groups:
  #group_feat="LPS" #to test
  sampleTable[,groupCol] <- as.factor(sampleTable[,groupCol])
  sgr <- paste0(groupCol,"_",levels(sampleTable[,groupCol]))
  vcols <- cols[as.numeric(as.factor(sampleTable[,groupCol]))]
  
  #we use it to generate average vectors for each position*peak
  #across all samples from each group

  maxyval=max(apply(dfin, 2,max))
  minyval=min(apply(dfin, 2,min))
  plot(xvals, dfin[,1],col=vcols[1],
       type="s",ylim=c(minyval,maxyval+0.02),ylab=ylab,
       xlab=xlab,lwd=lwd,main=main
  )
  
  legend("topright",lwd = lwd,legend = sgr,
         title = groupCol,bty = "n",col = cols)
  for (sa in 2:ncol(dfin)) {
    points(xvals, dfin[,sa],
           type="s",col=vcols[sa],lwd=2)
  }
  
}

#test:
#plot_densities(mean_posCPM,sampleTable = sample_table[1:10,],groupCol = "LPS",main="LPS_vs_noLPS")

plot_base_name <- paste0(setname,".meanCoverage")
gf="Cohesin"
comp <- paste(paste0(gf,levels(as.factor(sample_table[,gf]))),
              collapse = "vs")
for (lps in levels(as.factor(sample_table[,"LPS"]))) {
  setwd(outdir)
  od=paste0("LPS_",lps)
  dir.create(od)
  setwd(od)
  
  mean_posCPM_subset <- mean_posCPM[,sample_table$LPS==lps]
  sample_table_subset <- sample_table[sample_table$LPS==lps,] 
  #A) plot all samples individually
  graphics.off()
  pdf(paste(od,comp,plot_base_name,"pdf",sep = "."))
  plot_densities(mean_posCPM_subset,
                 sampleTable = sample_table_subset,
                 groupCol = gf)
  dev.off()
  
  #B) compute average for each group
  mean_posCPM_subset_grpAvg <- t(aggregate(t(mean_posCPM_subset),
                                           by=list(grp=sample_table_subset[,gf]),
                                           mean)[,-1])
  sample_tablex <- data.frame(levels(as.factor(sample_table_subset[,gf])))
  colnames(sample_tablex) <- gf          
  graphics.off()
  pdf(paste("Avg_samples",od,comp,plot_base_name,"pdf",sep = "."))
  plot_densities(mean_posCPM_subset_grpAvg,sampleTable = sample_tablex,groupCol = gf,ylab = "mean CPM (samples avg)")
  dev.off()
  
  }

