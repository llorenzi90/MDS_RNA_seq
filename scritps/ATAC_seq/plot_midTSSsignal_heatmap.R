## ---------------------------
##
##
## Purpose of script: Plot profile signal plots (heatmaps)
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
#outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/heatmaps/"
if(!dir.exists(outdir)) dir.create(outdir,recursive = T)

# sample data table
sample_data_file <- cArgs[4] 
#sample_data_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"

#feature to group data: #has to match a colname in the sample data table
group_feat <- cArgs[5]
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

###Compute mean CPM for each position*peak for each sample set:
#Define sample sets: has to be a factor of two groups:
#group_feat="LPS" #to test
sample_table[,group_feat] <- as.factor(sample_table[,group_feat])
sgr <- paste0(group_feat,"_",levels(sample_table[,group_feat]))

#we use it to generate average vectors for each position*peak
#across all samples from each group
grl <- as.list(levels(sample_table[,group_feat]))
names(grl) <- sgr
print("Generating average vectors for each position across all samples from each group")
list_Peak_pos_meanCPM <- lapply(grl, function(x){
  apply(as.data.frame(list_all_PeaksposCPM[sample_table[,group_feat]==x]),1,mean)
  #apply(as.data.frame(list_all_PeaksposCPM[sample_table[1:10,group_feat]==x]),1,mean)#to test
  
})

rm(list_all_PeaksposCPM)
gc()
###Convert the average vectors just created into matrices of genes vs position 
#for this we use the functions split "split divides the data in the vector x into the groups defined by f."
print("Converting average vectors into matrices of feature vs position")
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
print("sorting data to plot by mean total CPM (between both groups")
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
print("Converting values to log2 scale")
data_to_plot <- lapply(data_to_plot,function(x)return(log2(x+0.01)))


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

#palette.breaks=seq(-6,-2,0.02)
print("Generating color palette and column labels")
palette.breaks=seq(mi,ma,cada) #palette breaks for changing colors

palette.breaks.tmp=seq(round(min(as.data.frame(data_to_plot))-1),
                       round(max(as.data.frame(data_to_plot))+1),0.02) #palette breaks for all values
palette.breaks.tmp <-  round(palette.breaks.tmp,2)
color.palette.tmp=vector(length = (length(palette.breaks.tmp)-1))

#I will try 3 types of palette:
palette_list <- list(white_blue_orange=c("white","#0072B2","#D55E00"),
                     blue_white_orange=c("#0072B2","white","#D55E00"),
                     white_blue=c("white","#0072B2"))


#Now that we have our data ready to plot we need to adjust some plotting parameters:

##Generate column labels = -1.5, 1, 0.5 , 0 (TSS), 0.5, 1, 1.5 
labcol <- rep("",rwidth)
names(labcol)=as.character((llimit:ulimit)/1000)
labcol[1+c(0,1,2,3,4)*round(rwidth/4)] <- names(labcol[1+c(0,1,2,3,4)*round(rwidth/4)])
#points2label=c("-1","-0.5","0","0.5","1")
#labcol[points2label] <- points2label

setwd(outdir)

plot_base_name <- paste0(paste(names(data_to_plot),collapse = "_vs_"),".",setname,".TSSs.log2.heatmap")

print("Running loop to make plots with different color scales")
for (pname in names(palette_list)) {
  print(paste("Color palette:",pname))
  color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
  #Set the palette and palette breaks so all plots have the same range of colors
  color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
  color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
  color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
  
  #################################HEATMAPS#####################
  heats <- lapply(1:2,function(n){
    if(n==1) showleg=T else showleg=F
    m=data_to_plot[[n]]
    return(ComplexHeatmap::pheatmap(m,cluster_rows = F,
                                    cluster_cols = F,
                                    show_colnames = T,
                                    labels_col = labcol, 
                                    show_rownames = F,
                                    color = color.palette.tmp,
                                    breaks = palette.breaks.tmp,
                                    use_raster=F,
                                    legend = showleg,
                                    fontsize_number = 9, column_names_side = c("top"),
                                    angle_col = c("0"),
                                    heatmap_legend_param = list(title="log2(CPM+0.01)"),
                                    main = names(data_to_plot)[n]))
  } )
  
  graphics.off()
  pdf(paste0(plot_base_name,".colorscale_",pname,".pdf"))
  print(heats[[1]] + heats[[2]])
  dev.off()
  
  graphics.off()
  tiff(paste0(plot_base_name,".colorscale_",pname,".tiff"),
  width = 8, height = 8, units = 'in', res = 300)
  print(heats[[1]] + heats[[2]])
  dev.off()
  
  
}

