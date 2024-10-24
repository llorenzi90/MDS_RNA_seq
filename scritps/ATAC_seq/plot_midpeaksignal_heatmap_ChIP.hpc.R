## ---------------------------
##
##
## Purpose of script: Plot profile signal plots (heatmaps)
##                    around ENCODE ChIP-seq peaks 
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-22
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

#run: Rscript plot_midpeaksignal_heatmap_ChIP.R peaks_file chain_file datadir outdir sample_data_file Cohesin
#define input data and paths:
cArgs=commandArgs(trailingOnly = T)

#we need: folder with bw files
# peaks file
peaks_file <- cArgs[1] #example "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/ENCODE_datasets/ATF3_K562_ENCFF950AHH.bigBed"
                       #example 2   "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/ENCODE_datasets/IRF_K562_ENCFF943JSO.bed.gz"  
#peaks_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/ENCODE_datasets/SMC3_K562_ENCFF339LVM.bigBed"
setname <- gsub("\\.[a-zA-Z0-9]*$","",gsub("\\.gz$","",basename(peaks_file)))

# chain file
chain_file <- cArgs[2] #
#chain_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/liftover_chains/hg19ToHg38.over.chain"

#datadir
datadir <- cArgs[3]
#datadir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"

# outdir
outdir <- cArgs[4]
#outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/heatmaps/"


# sample data table
sample_data_file <- cArgs[5] 
#sample_data_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"

#feature to group data: #has to match a colname in the sample data table
group_feat <- cArgs[6]
#group_feat <- "Cohesin"

#file_type 
#filextension <- cArgs[5]
filextension <- ".bw" #currently only bigwig files supported

#some parameters that for now I will keep fixed but that we may 
#want to play around with in the future:
Ntop=cArgs[7] #top peaks to pick based on macs2 qval
#
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



######Read and process ENCODE peaks####
##Read ENCODE peaks 
if(grepl(".bed",basename(peaks_file))) {
  encodepeaks <- import(peaks_file,format = "narrowPeak")
  } else { if(grepl(".bigBed",basename(peaks_file))){
    encodepeaks <- import(peaks_file)
  }else warning("Currently non-supported peak input file format")
  }


#Keep only peaks in canonical chrs
canonical_chrs <- paste0("chr",c(seq(1:22),"X","Y"))
print("peaks in non-canonical chromosomes")
print(table(!seqnames(encodepeaks)%in%canonical_chrs))
print("removing peaks in following non-canonical chrs:")
print(unique(as.character(seqnames(encodepeaks)[!seqnames(encodepeaks)%in%canonical_chrs])))
encodepeaks <- encodepeaks[seqnames(encodepeaks)%in%canonical_chrs]

###Sort ENCODE peaks in order of q-val first, then score, then signalValue
sortedpeaks <- encodepeaks[order(encodepeaks$qValue,encodepeaks$score,encodepeaks$signalValue,decreasing = T)]

###Remove overlapping peaks prioritizing the peaks with higher rank (that's why we sorted in previous step)
tfo=findOverlaps(sortedpeaks,sortedpeaks)
overlaps=data.frame(q=queryHits(tfo),s=subjectHits(tfo))
overlaps <- overlaps[overlaps$q!=overlaps$s,]
#first remove the repeated rows, keep only those in which q is smaller than s 
overlapsf <- overlaps[overlaps$q<overlaps$s,]
#now the filtered overlaps table
#has unique comparisons between overlapping peaks
#in which the first value is always the best ranked
#because the encodepeaks were previously sorted
#then peaks in filtered column 2 (s) are the ones to remove
rows_to_remove <- unique(overlapsf$s)
if(length(rows_to_remove)>0){
  sortedpeaks_noOverlaps <- sortedpeaks[-rows_to_remove]}else{
    sortedpeaks_noOverlaps <- sortedpeaks
  }

###Extract mid point of peaks and preform liftOver to hg38 coordinates (coverage data is hg38) 
chrs=seqnames(sortedpeaks_noOverlaps)
coords=ranges(sortedpeaks_noOverlaps)
#midpoints <- start(coords) + round( width(coords)/2) #wrong
midpoints <- start(coords) + round( width(coords)/2) -1

midpointsGR <- GRanges(seqnames = chrs,
                       ranges = IRanges(start = midpoints,width = 1))

names(midpointsGR) <- paste0(seqnames(midpointsGR),"_",
                             start(midpointsGR))
##lifOver:
chain=import.chain(chain_file)
midpoints_sortedpeaks_noOVerlaps_liftover <- liftOver(midpointsGR,chain)

names(midpoints_sortedpeaks_noOVerlaps_liftover) <- names(midpointsGR)
midpointsGR_liftOver <- unlist(midpoints_sortedpeaks_noOVerlaps_liftover)

lost_peaks <- names(midpoints_sortedpeaks_noOVerlaps_liftover)[!names(midpoints_sortedpeaks_noOVerlaps_liftover)%in%names(midpointsGR_liftOver)]
print("lost peaks:")
print(lost_peaks)#some sequences are lost after the liftover

print("rank of lost peaks:")
print(which(names(midpoints_sortedpeaks_noOVerlaps_liftover)%in%lost_peaks))
#I need to filter out those peaks that have no possible liftover:

midpointsGR_noliftOver <- midpointsGR[names(midpointsGR)%in%names(midpointsGR_liftOver)]
print("length midpointsGR")
print(length(midpointsGR))
print("length midpointsGR_liftOver")
print(length(midpointsGR_liftOver))
print("length midpointsGR_noliftOver")
print(length(midpointsGR_noliftOver))

###Generate the extended ranges with the liftovered coordinates
chrs=seqnames(midpointsGR_liftOver)
midpoints <- start(midpointsGR_liftOver)
extended_granges <- GRanges(chrs,
                            IRanges(start = midpoints + llimit,
                                    width = rwidth))

##Need to filter overlaps again:
tfo=findOverlaps(extended_granges,extended_granges)
overlaps=data.frame(q=queryHits(tfo),s=subjectHits(tfo))
overlaps <- overlaps[overlaps$q!=overlaps$s,]
overlapsf <- overlaps[overlaps$q<overlaps$s,]
rows_to_remove <- unique(overlapsf$s)
if(length(rows_to_remove)>0){
  extended_granges_noOverlaps <- extended_granges[-rows_to_remove]}else{
    extended_granges_noOverlaps <- extended_granges
  }

#take the top Ntop extended granges or all peaks if indicated or if total peaks less than Ntop
if(Ntop!="all"){
  if(Ntop>length(extended_granges_noOverlaps)) {print(paste0(Ntop," is larger than total peaks. Using all ",length(extended_granges_noOverlaps)," peaks with no overlaps")) 
  Ntoppeaks <- extended_granges_noOverlaps 
  }else Ntoppeaks <- extended_granges_noOverlaps[1:Ntop]
}else Ntoppeaks <- extended_granges_noOverlaps


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
  print(paste0("processing sample ",i))
  tmp_tssCPM=import(bw,which=Ntoppeaks)
  #extend the GRanges object so each position around -Xkb - +X kb of each peak has a CPM value
  all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
  list_all_PeaksposCPM[[i]] <- all_TSSposCPM
  print(object.size(list_all_PeaksposCPM))
}
print("Finished data import\nChecking if incomplete peaks")
#####Check that each feature have all positions. If incomplete, remove!####
# generate some indexes that will be useful to remove incomplete 
#peaks in case they happen:
#1) generate single-position ids in the peaks object:
posIDorigin <- paste0(as.character(rep(seqnames(Ntoppeaks),each=rwidth)),"_",
                      unlist(apply(as.data.frame(ranges(Ntoppeaks)),1,function(x)return(x[1]:x[2]))))
#2)generate corresponding peakIDs for each position in the original peaks:
chrs_origin <- as.character(rep(seqnames(Ntoppeaks),each=rwidth))
ranges_origin <- rep(unlist(apply(as.data.frame(ranges(Ntoppeaks)),
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
  print("peaks to plot")
  
}else{
  print("There are")
  print(lapply(list_all_PeaksposCPM, function(x)length(x)/rwidth))
  print("peaks to plot")
}
#####

####A) For average densities: 
###Compute mean CPM for each position*gene for each sample:
print("Computing mean CPM for each position-feature for each sample")
mean_posCPM <- lapply(list_all_PeaksposCPM, function(x)aggregate(x,by=list(pos=rep(1:rwidth,length(x)/rwidth)),mean))
mean_posCPM <- as.data.frame(lapply(mean_posCPM,function(s)return(s$x)))
if(!dir.exists(outdir)) dir.create(outdir,recursive = T)
setwd(outdir)
rds_base_name <- paste0(setname,".top",Ntop,"peaks.mean_sampleCPM.rds")
saveRDS(mean_posCPM,rds_base_name)

###Compute mean CPM for each position*peak for each sample set:
#Define sample sets: has to be a factor of two groups:
print("Computing mean CPM for each position-peak for each group")
sample_table[,group_feat] <- as.factor(sample_table[,group_feat])
sgr <- paste0(group_feat,"_",levels(sample_table[,group_feat]))

#we use it to generate average vectors for each position*peak
#across all samples from each group
grl <- as.list(levels(sample_table[,group_feat]))
names(grl) <- sgr
list_Peak_pos_meanCPM <- lapply(grl, function(x){
  apply(as.data.frame(list_all_PeaksposCPM[sample_table[,group_feat]==x]),1,mean)
 # apply(as.data.frame(list_all_PeaksposCPM[sample_table[1:20,group_feat]==x]),1,mean)#to test
  
})

rm(list_all_PeaksposCPM)
gc()
###Convert the average vectors just created into matrices of genes vs position 
#for this we use the functions split "split divides the data in the vector x into the groups defined by f."
data_to_plot <- lapply(list_Peak_pos_meanCPM,function(x){
  t(as.data.frame(split(x, 
                        ceiling(seq_along(x)/rwidth))))})


##Next we want to transform our data in 2 different ways taking into
#account that both matrices will be plotted together, so they need to have
#the same gene order and the same color scale:
#   1) sort rows in both matrices by mean total CPM signal between Cohesin and non-Cohesin
#   2) convert values into log2 scale to plot
#   3) generate color palette

#   1) sort data to plot by mean total CPM
total_CPMs <- lapply(data_to_plot, function(x)apply(x, 1,sum)) #for each element of our list,
# for each peak (1 indicates rows) takes the sum of all CPMs (total CPM per row)
#now, for each peak take the average of total CPMs per row between both groups
mean_totalCPM <- apply(as.data.frame(total_CPMs),1,mean)
#Finally, we use this vector to sort the rows in both matrices in the same way
rowOrder=order(mean_totalCPM,decreasing = T)# 
data_to_plot <- lapply(data_to_plot, function(x){
  return(x[rowOrder,])
  
})

#2)   Convert values into log2 scale to plot
data_to_plot <- lapply(data_to_plot,function(x)return(log2(x+0.01)))
rds_base_name <- paste0(paste(names(data_to_plot),collapse = "_vs_"),".",setname,".top",Ntop,"peaks.log2.datatoplot.rds")
# print("data_to plot size:")
# print(object.size(data_to_plot))
saveRDS(data_to_plot, file = rds_base_name)

#3)   Create palette to apply to both plots using the absolute minimum and max values across both matrices
#a) vectorize both matrices and check values. This part is a bit manual
# The idea is to select a range of values where most of the CPMs fall into
# and use them to scale the colors of the palette

#plot(density(as.matrix(as.data.frame(data_to_plot))))
#summary(as.vector(as.matrix(as.data.frame(data_to_plot))))
#For now I have not implemented the automatization of this, but would like to do in the future
# for now the range of values to vary the color scale is an input parameter:

# range_col_values=c(-6,0)
# cada=0.02
# mi=range_col_values[1]
# ma=range_col_values[2]
# 
# #palette.breaks=seq(-6,-2,0.02)
# palette.breaks=seq(mi,ma,cada) #palette breaks for changing colors
# 
# palette.breaks.tmp=seq(round(min(as.data.frame(data_to_plot))-1),
#                        round(max(as.data.frame(data_to_plot))+1),0.02) #palette breaks for all values
# palette.breaks.tmp <-  round(palette.breaks.tmp,2)
# color.palette.tmp=vector(length = (length(palette.breaks.tmp)-1))
# 
# #I will try 3 types of palette:
# palette_list <- list(white_blue_orange=c("white","#0072B2","#D55E00"),
#                      blue_white_orange=c("#0072B2","white","#D55E00"),
#                      white_blue=c("white","#0072B2"))
# 
# 
# #Now that we have our data ready to plot we need to adjust some plotting parameters:
# 
# ##Generate column labels = -1.5, 1, 0.5 , 0 (TSS), 0.5, 1, 1.5
# labcol <- rep("",rwidth)
# names(labcol)=as.character((llimit:ulimit)/1000)
# labcol[1+c(0,1,2,3,4)*round(rwidth/4)] <- names(labcol[1+c(0,1,2,3,4)*round(rwidth/4)])
# #points2label=c("-1","-0.5","0","0.5","1")
# #labcol[points2label] <- points2label
# 
# if(!dir.exists(outdir)) dir.create(outdir,recursive = T)
# setwd(outdir)
# 
# plot_base_name <- paste0(paste(names(data_to_plot),collapse = "_vs_"),".",setname,".top",Ntop,"peaks.log2.heatmap")
# 
# for (pname in names(palette_list)) {
#   color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
#   #Set the palette and palette breaks so all plots have the same range of colors
#   color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
#   color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
#   color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
# 
#   #################################HEATMAPS#####################
#   heats <- lapply(1:2,function(n){
#     if(n==1) showleg=T else showleg=F
#     m=data_to_plot[[n]]
#     return(ComplexHeatmap::pheatmap(m,cluster_rows = F,
#                                     cluster_cols = F,
#                                     show_colnames = T,
#                                     labels_col = labcol,
#                                     show_rownames = F,
#                                     color = color.palette.tmp,
#                                     breaks = palette.breaks.tmp,
#                                     use_raster=F,
#                                     legend = showleg,
#                                     fontsize_number = 9, column_names_side = c("top"),
#                                     angle_col = c("0"),
#                                     heatmap_legend_param = list(title="log2(CPM+0.01)"),
#                                     main = names(data_to_plot)[n]))
#   } )
# 
#   #graphics.off()
#   print("heatmap list object size:")
#   print(object.size(heats))
#   saveRDS(heats,paste0(plot_base_name,".colorscale_",pname,".heatmap.rds"))
#    # pdf(paste0(plot_base_name,".colorscale_",pname,".pdf"))
#    # print(heats[[1]] + heats[[2]])
#    # dev.off()
#   # 
#   # graphics.off()
#   # tiff(paste0(plot_base_name,".colorscale_",pname,".tiff"), width = 8, height = 8, units = 'in', res = 300)
#   # print(heats[[1]] + heats[[2]])
#   # dev.off()
# 
# 
# }

