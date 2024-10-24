## ---------------------------
##
##
## Purpose of script: plot heatmap ATAC-seq profile for differential peaks
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-02-10
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
library(ComplexHeatmap)
library(rtracklayer)
## ---------------------------


#Input parameters
#
sort_by="CPM"
region_to_plot=c(-1500,1500)


##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)

#load DESeq2 results
deseq2res=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/DESeq2/all_samples/cohesin_vs_nocohesin/all_samples_cohesin_vs_nocohesin_results.csv",header = T)
padjcutoff=0.05
upPeaks <- deseq2res %>% filter(padj<=padjcutoff , !is.na(log2FoldChange) ,log2FoldChange>0)

#take mid point of peaks
chrs=sapply(strsplit(upPeaks$X,split = ":"),function(x)return(x[1]))
coords=sapply(strsplit(upPeaks$X,split = ":"),function(x)return(x[2]))
coords=strsplit(coords,split = "-")
midpoints=sapply(coords, function(x){
  x=as.numeric(x) 
  return(x[1] + round((x[2] - x[1] +1)/2))})

#select coordinates from differential peaks 

#generate granges object with -3000 - +3000 bps around the 200 TSSs
extended_granges <- GRanges(chrs,
                            IRanges(start = midpoints + llimit,
                                    width = rwidth))

#for each sample, import the regions of the bigwig files matching the regions around TSSs 
datadir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"

bwfiles=list.files(datadir,full.names = T,pattern = "bam.CPM.bw")

mut_data=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/mut_data_all_patients.txt",sep = "\t",header = T)

samples=gsub(".sorted.markedDups.primary_chr.proper_pairs.minq2.bam.CPM.bw","",basename(bwfiles))
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=gsub("AT-","",sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1]))
Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]
Cohesin=patient%in%Cohesin_patients

list_all_upPeaksposCPM <- list()
i=0
for (bw in bwfiles) {
  i=i+1
  tmp_tssCPM=import(bw,which=extended_granges)
  #extend the GRanges object so each position around -Xkb - +X kb of each TSS has a CPM value
  all_TSSposCPM <- rep(tmp_tssCPM$score,width(ranges(tmp_tssCPM)))
  list_all_upPeaksposCPM[[i]] <- all_TSSposCPM
}

#save this, as it takes time to run
dir2savedata="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/Rdata/"

dir.create(dir2savedata)
write_rds(list_all_upPeaksposCPM,paste0(dir2savedata,"upPeaks.1.5kbAroundmidpoint_allSamples.CPM.rds"))

#checking output
sapply(list_all_upPeaksposCPM,length)
length(list_all_upPeaksposCPM[[1]])/rwidth
#933.8727
#the length should be equal to Npeaks*rwidth, 936*3001=2808936
#however the number of total position is lower in 6384 bases
#Thinking about it I realised this is probably because some peaks, when extended 1.5 to each side
#from the middle, end up overlapping adjacent peaks
library()
test=findOverlaps(extended_granges,extended_granges)
overlaps=data.frame(q=queryHits(test),s=subjectHits(test))
overlaps <- overlaps[overlaps$q!=overlaps$s,]
#there are 14 overlapping peaks. What to do with them?
# first check pval of these peaks

View(upPeaks[overlaps$q,])
# check overlap if we take shorter region around peak center
for (newllim in c(-1000,-500)) {
  newrwidth=length(newllim:-(newllim))
  tmp_gr=extended_granges <- GRanges(chrs,
                                     IRanges(start = midpoints + newllim,
                                             width = newrwidth))
  tmpov=findOverlaps(tmp_gr,tmp_gr)
  tmpovdf=data.frame(q=queryHits(tmpov),s=subjectHits(tmpov))
  tmp_filt <- tmpovdf[tmpovdf$q!=tmpovdf$s,]
  print(length(unique(tmp_filt$q)))
}

#Option 1: duplicate the values of positions that are duplicated
#wont do this for now
#Option 2: remove the 14 overlapping peaks
#to avoid re-doing all the processing that takes time and memory
upPeaksID=upPeaks$X
upPeaksID_to_remove=upPeaksID[overlaps$q]
all_pos=apply(as.data.frame(extended_granges)[,2:3],1,function(x)seq(from=x[1],to=x[2],by=1))
test=matrix(data = c(1,2,3,4,5,6,7,8),nrow = 4,ncol = 2)
test
as.vector(test)
as.vector(all_pos)[2999:3003]
all_chr_position=paste(rep(chrs,each=rwidth),as.vector(all_pos),sep = ":")
all_chr_position[2999:3003]
table(duplicated(all_chr_position))
peakID_per_pos=rep(upPeaksID,each=rwidth)
undup_chr_pos=all_chr_position[!duplicated(all_chr_position)]
undup_peakIDs=peakID_per_pos[!duplicated(all_chr_position)]
filtercond=!undup_peakIDs%in%upPeaksID_to_remove
#use this filter to filter the list of positions. 
filtered_list_all_upPeaksposCPM = lapply(list_all_upPeaksposCPM,function(x)return(x[filtercond]))
list_all_upPeaksposCPM=filtered_list_all_upPeaksposCPM
rm(filtered_list_all_upPeaksposCPM)
#Option 3: 
#use shorter distance around peak center
#i can try it later

#list_all_upPeaksposCPM <- read_rds("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/Rdata/200inflammatoryGenes.1.5kbAroundTSS_allSamples.CPM.rds")
#this list has one element for each sample
#each element has the CPM value for each position across the 200 genes

###Compute mean CPM for each position*gene for each sample set:
#the logical vector Cohesin, generated before indicates if the patient 
#has Cohesin mutations or not
#we use it to generate average vectors for each position*gene 
#across all samples from each group
#dim(as.data.frame(list_all_TSSposCPM[Cohesin]))
#dim(as.data.frame(list_all_TSSposCPM[!Cohesin]))
Cohesin_TSSpos_meanCPM <- apply(as.data.frame(list_all_upPeaksposCPM[Cohesin]),1,mean)
NONCohesin_TSSpos_meanCPM <- apply(as.data.frame(list_all_upPeaksposCPM[!Cohesin]),1,mean)
rm(list_all_upPeaksposCPM)
gc()
###Convert the average vectors just created into matrices of genes vs position 
#for this we use the functions split "split divides the data in the vector x into the groups defined by f."
Cohesin_meanCPM_all_genes <- t(as.data.frame(split(Cohesin_TSSpos_meanCPM, 
                                                   ceiling(seq_along(Cohesin_TSSpos_meanCPM)/rwidth))))
NONCohesin_meanCPM_all_genes <- t(as.data.frame(split(NONCohesin_TSSpos_meanCPM, 
                                                      ceiling(seq_along(NONCohesin_TSSpos_meanCPM)/rwidth))))
dim(Cohesin_meanCPM_all_genes)
dim(NONCohesin_meanCPM_all_genes)

###Place both matrices into a list, to make following computations easier
data_to_plot=list(Non_Cohesin_patients=NONCohesin_meanCPM_all_genes,
                  Cohesin_patients=Cohesin_meanCPM_all_genes)

###Next we want to transform our data in 2 different ways taking into
#account that both matrices will be plotted together, so they need to have
#the same gene order and the same color scale:
#   1) sort rows in both matrices by mean total CPM signal between Cohesin and non-Cohesin
#   2) convert values into log2 scale to plot


#   1) sort data to plot by mean total CPM
total_CPMs <- lapply(data_to_plot, function(x)apply(x, 1,sum)) #for each element of our list,
# for each gene (1 indicates rows) takes the sum of all CPMs
total_CPMs$Non_Cohesin_patients
total_CPMs$Cohesin_patients
#now, for each gene take the average of total CPMs between both groups
mean_totalCPM <- apply(as.data.frame(total_CPMs),1,mean)
#Finally, we use this vector to sort the rows in both matrices in the same way
rowOrder=order(mean_totalCPM,decreasing = T)# function "order" retrieves 
#a numeric vector of rearrenged indexes following increasing(default) 
#or decreasing order (we want genes with highest signal on top, so we choose decresing)
data_to_plot <- lapply(data_to_plot, function(x){
  return(x[rowOrder,])
  
})


#let's make a quick test before converting to log2
hmcols<-colorRampPalette(c("#0072B2","white","#D55E00"))(256)
#Generate column labels = -1.5, 1, 0.5 , 0 (TSS), 0.5, 1, 1.5 
labcol <- rep("",rwidth)
points2label=c("-1.5","-1","-0.5","0","0.5","1","1.5")
names(labcol)=as.character((llimit:ulimit)/1000)
labcol[points2label] <- points2label


heats <- lapply(data_to_plot,function(m) ComplexHeatmap::pheatmap(m,cluster_rows = F,
                                                                  cluster_cols = F,
                                                                  show_colnames = T,
                                                                  labels_col = labcol,
                                                                  color = hmcols,
                                                                  show_rownames = F,
                                                                  use_raster=F,
                                                                  fontsize_number = 9, 
                                                                  column_names_side = c("top"),
                                                                  angle_col = c("0")))

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/inflammatory_genes/heatmaps/")

pdf("test.ATAC_seq_UPpeaks.nolog2.pdf")
heats[[1]] +heats[[2]]
dev.off()

#2)   Convert values into log2 scale to plot
data_to_plot <- lapply(data_to_plot,function(x)return(log2(x+0.01)))

#Now that we have our data ready to plot we need to adjust some plotting parameters:

#let's make a quick test
hmcols<-colorRampPalette(c("#0072B2","white","#D55E00"))(256)

heats <- lapply(data_to_plot,function(m) ComplexHeatmap::pheatmap(m,cluster_rows = F,
                                                                  cluster_cols = F,
                                                                  show_colnames = T,
                                                                  labels_col = labcol,
                                                                  color = hmcols,
                                                                  show_rownames = F,
                                                                  use_raster=F,
                                                                  fontsize_number = 9, 
                                                                  column_names_side = c("top"),
                                                                  angle_col = c("0")))

setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/inflammatory_genes/heatmaps/")

# pdf("test.log2.pdf")
# heats[[1]] +heats[[2]]
# dev.off()


# heats <- lapply(data_to_plot,function(m) ComplexHeatmap::pheatmap(m,cluster_rows = F,
#                                                                   cluster_cols = F,
#                                                                   show_colnames = T,
#                                                                   labels_col = labcol, 
#                                                                   show_rownames = F,
#                                                                   color = color.palette.tmp,
#                                                                   breaks = palette.breaks.tmp,
#                                                                   use_raster=F,
#                                                                   fontsize_number = 9, column_names_side = c("top"),
#                                                                   angle_col = c("0")))

### Adjust color scale
#3)   Create palette to apply to both plots using the absolute minimum and max values across both matrices
#a) vectorize both matrices and check values. This part is a bit manual
# The idea is to select a range of values where most of the CPMs fall into
# and use them to scale the colors of the palette

plot(density(as.matrix(as.data.frame(data_to_plot))))
summary(as.vector(as.matrix(as.data.frame(data_to_plot))))
#let's try with the quartiles first
# range_col_values=c(round(quantile(as.vector(as.matrix(as.data.frame(data_to_plot))),0.25)),
#                    round(quantile(as.vector(as.matrix(as.data.frame(data_to_plot))),0.75)))
# #looks ugly...


range_col_values=c(-6,0)
cada=0.02
mi=range_col_values[1]
ma=range_col_values[2]
#palette.breaks=seq(-6,-2,0.02)
palette.breaks=seq(mi,ma,cada) #palette breaks for changing colors

palette.breaks.tmp=seq(round(min(as.data.frame(data_to_plot))-1),
                       round(max(as.data.frame(data_to_plot))+1),0.02) #palette breaks for all values
palette.breaks.tmp <-  round(palette.breaks.tmp,2)
color.palette.tmp=vector(length = (length(palette.breaks.tmp)-1))
which(palette.breaks.tmp==as.character(ma - cada))
#try 3 types of palette:
palette_list <- list(white_blue_orange=c("white","#0072B2","#D55E00"),
                     blue_white_orange=c("#0072B2","white","#D55E00"),
                     white_blue=c("white","#0072B2"))

for (pname in names(palette_list)) {
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
  pdf(paste0("NONCohesin_vs_Cohesin.log2.heatmap.colorscale_",pname,".ATAC_seqUPpeaks.pdf"))
  print(heats[[1]] + heats[[2]])
  dev.off()
  
}



png(paste0("NONCohesin_vs_Cohesin.log2.heatmap.colorscale_",pname,".ATAC_seqUPpeaks.png"))
print(heats[[1]] + heats[[2]])
dev.off()
