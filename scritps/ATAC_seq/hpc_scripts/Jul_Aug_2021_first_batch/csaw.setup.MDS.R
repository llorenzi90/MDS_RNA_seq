
#create bam.paths vector
library(csaw)
bam.paths <- c()
samples <- read.table("/home/llorenzi/MDS_files.txt")

bam.paths <- paste0("/scratch/llorenzi/MDS/ATAC-seq/",samples$V1,"/",samples$V1,".primary_chr.quality.concordant.sorted.markedDups.bam")

#define filtering parameters
#ENCODE black listed regions
#BL <- read.table("~/share/Cuartero Group/CUARTERO GROUP/references/human/ENCODE_blacklist/ENCFF356LFX.bed")
BL <- read.table("/home/llorenzi/references/human/ENCFF356LFX.bed")
blacklist <- GRanges(BL$V1,IRanges(BL$V2,BL$V3))

param <- readParam(minq = 20, pe = "both", dedup=TRUE,discard = blacklist)
 #another parameter that could be set is max.frag = e.g 200...


#Obtaining window-level counts 
win.data <- windowCounts(bam.paths, param=param, width=170)
#win.data
saveRDS(win.data,"win.data.MDS.RDS")
#win.data=readRDS("win.data.MDS.RDS")

bin.size=2000L
binned <- windowCounts(bam.paths, bin=TRUE, width=bin.size, param=param)

saveRDS(binned,"binned_counts_MDS.2K.RDS")

surrounds <- 2000
neighbor <- suppressWarnings(resize(rowRanges(win.data), surrounds, fix="center"))
wider <- regionCounts(bam.paths, regions=neighbor, param=param)

saveRDS(wider,"neighbour_counts_MDS.2K.RDS")


#add larger bin counts
bin.size=10000L
binned <- windowCounts(bam.paths, bin=TRUE, width=bin.size, param=param)
saveRDS(binned,"binned_counts_MDS.10K.RDS")



