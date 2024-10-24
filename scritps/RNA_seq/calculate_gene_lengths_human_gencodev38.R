## ---------------------------
##
##
## Purpose of script: calculate gene lenght from a GTF file
##                    in this case the gencode mouse m39 annotation
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-05-28
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes: script adapted from https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
## suggested by Devon Ryan in this post: https://bioinformatics.stackexchange.com/questions/2567/how-can-i-calculate-gene-length-for-rpkm-calculation-from-counts-data
##   
##
## ---------------------------

options(scipen = 999) 

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)

## ---------------------------


#!/usr/bin/Rscript
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

GTFfile = "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/annotation/gencode/gencode.v38.annotation.gtf"
#FASTAfile = "/home/ryand/Documents/Misc/Mus_musculus/Ensembl/GRCm38.71/Sequence/Mus_musculus.GRCm38.71.dna.toplevel.fa"

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="GRCh38.38", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

#Open the fasta file
#FASTA <- FaFile(FASTAfile)
#open(FASTA)

#Add the GC numbers
#elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
# calc_GC_length <- function(x) {
#   nGCs = sum(elementMetadata(x)$nGCs)
#   width = sum(elementMetadata(x)$widths)
#   c(width, nGCs/width)
# }

calc_length <- function(x) {
  
  width = sum(elementMetadata(x)$widths)
  return(width)
}


output <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length)
#colnames(output) <- c("Length", "GC")
outfile=gsub("annotation.gtf","gene_lengths.tsv",GTFfile)
write.table(output, outfile, sep="\t", col.names = F,quote = F)
