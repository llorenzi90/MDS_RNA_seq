##########################################################################################
########### 6. performing in silico QC using the ATACseqQC package

## Getting chromosome length information
#run this on terminal:
#samtools view -H /Users/llorenzi/MDS/ATAC-seq/data/bam_files/AT-MDS35.1.primary_chr.quality.concordant.sorted.markedDups.bam | grep -P 'chr\d+|chrX' | cut -f2,3  | perl -n -e 's/[SL]N://g' | \
#awk 'BEGIN{FS=OFS="\t"} {print $1, 1, $2, "hg38"}' > auto.x.chrom.human.txt
#NOTE: this doesn't work in mac because grep -P option is not available, I ran it in linux

##########################################################################################
#################### 6.1 Creating ATACseqQC R scripts: ATACseqQC.R 


## This R script accept a single command line argument - the name of a BAM file for QC analysis. 
## An index file for the BAM file must be included in the same directory.

## Loading all required packages

library(motifStack)
library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
library(ChIPpeakAnno)


## Getting the BAM file name and sample ID
args<- commandArgs(trailingOnly=TRUE)
bamfile.sample.ID <- args[1]
setwd("/scratch/llorenzi/MDS/ATAC-seq/")
setwd(bamfile.sample.ID)
bamfile <- paste0(bamfile.sample.ID,".primary_chr.quality.concordant.sorted.markedDups.bam")
#bamfile.sample.ID <- gsub(".bam", "", basename(bamfile))

## Plotting size distribution of fragments (Figure 1G)
pdf(paste0(bamfile.sample.ID, ".fragment.size.distribution.pdf"), width =10, height=8) 
fragSize <- fragSizeDist(bamFiles=bamfile, bamFiles.labels=bamfile.sample.ID)
dev.off()


## BAM file tags to be included when read in by the readBamFile
tags <- c("AS", "NM", "MD")

## Create a output directory where all subsequent BAM files will be output
outPath <- paste0(bamfile.sample.ID, ".splited.bam")
if (!dir.exists(outPath)){
  dir.create(outPath)
}

## Build GRanges for the human genome hg38 excluding unplaced scaffoldsand chrY
human.genome <- read.delim("/home/llorenzi/auto.x.chrom.human.txt", header=F)

seqlev <- as.character(human.genome[,1])
gr <- GRanges(as.character(human.genome[,1]), IRanges(human.genome[,2], human.genome[,3]))

## For QC, use read alignments from chromosomes 1 and 2 as representatives
which <- gr[seqnames(gr) %in% c("chr1", "chr2")]

## Reading in paired end read alignment
gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)

## Shifting the coordinates of 5' ends of the aligned reads in the bam file, +4 for reads mapping
## to the positive strand, -5 for reads mapping to the negative strand
gal1 <- shiftGAlignmentsList(gal)
shiftedBamfile <- file.path(outPath, paste0(bamfile,".shifted.bam"))
## Outputting the shifted BAM file
export(gal1, shiftedBamfile)

# gal1 <- readBamFile(shiftedBamfile,tag = tags,which = which,asMates = FALSE)
# names(gal1) <- mcols(gal1)$qname
# 
# 
# BiocManager::install("phastCons100way.UCSC.hg38")
# 
# library(phastCons100way.UCSC.hg38)

## Getting information of known transcripts from UCSC genome browser database
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)


## Classifying reads into nucleosome-free, mono-, di- and tri-nucleosome bins based on their fragment sizes.
genome <- Hsapiens
saveRDS(genome,"Hsapiens.RDS")
genome2 <- readRDS("Hsapiens.RDS")
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome)


## outputting split BAM files
null <- writeListOfGAlignments(objs, outPath)

## Generating heatmaps and coverage curves for nucleosome-free, mono-, di- and tri-nucleosome occupied regions

bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))

## extracting TSSs coordinates
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

## estimating the library size for normalization
librarySize <- estLibSize(bamfiles)


## calculating the signals around TSSs.
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(bamfiles, TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)


## log2-transforming signals
names(sigs) <- gsub(".bam", "", basename(names(sigs)))
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

## plotting heatmap showing signals  for nucleosome-free and oligonucleosome-bound regions around TSSs. (Figure 1 H and 1I)
pdf(paste0(bamfile.sample.ID, ".heatmap and averaged coverage.pdf"))
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)

out <- featureAlignedDistribution(sigs, 
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l", 
                                  yalb="Averaged coverage", 
                                  ylab="Averaged coverage")
dev.off()


## Plotting CTCF footprints. (Figure 2C)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
print(CTCF[[1]], digits=2)
seqlev <- c("chr1", "chr2")

nucleosome_free_bamfile <- file.path(outPath, "NucleosomeFree.bam")

pdf(paste0(bamfile.sample.ID, ".CTCF.footprint.pdf"))
factorFootprints(nucleosome_free_bamfile, pfm=CTCF[[1]], 
                 genome=genome,
                 min.score="95%", seqlev=seqlev,
                 upstream=100, downstream=100)
dev.off()

