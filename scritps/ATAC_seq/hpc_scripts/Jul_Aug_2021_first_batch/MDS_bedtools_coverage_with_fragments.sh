#!/bin/bash
##SBATCH -J 
##SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mem=32G
##SBATCH -p std
#SBATCH --ntasks=10
#SBATCH --time=10:05:30

#######		Script to process info on per-sample peaks		#######
samp=$1
echo $samp

PPN=10
#move to sample folder and get file names
cd /scratch/llorenzi/MDS/ATAC-seq/$samp


#sample peaks
peaks=$samp.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak
#sample bam
bam=$samp.primary_chr.quality.concordant.sorted.markedDups.bam
#consensus peak file
conspeaks=/scratch/llorenzi/MDS/ATAC-seq/merged_peaks/merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed

#load modules
#a) samtools to remove duplicates and sort by name
module load samtools/1.9
#b) bedtools to convert alignments to fragment bed and to compute coverage
module load bedtools2/2.29.0


#1)generate fragment bed
samtools view -@ $PPN -b -F1024 $bam|samtools sort -@ $PPN -n |bedtools bamtobed -bedpe |cut -f1,2,6-10 > $samp.primary_chr.quality.concordant.sorted.noDups.fragments.bed

bed=$samp.primary_chr.quality.concordant.sorted.noDups.fragments.bed

#2) compute coverage with fragments bed file
#1.1)coverage over own sample peak file
bedtools coverage -a $peaks -b $bed > $samp.narrowPeak.fragment.counts
#1.2)run bedtools coverage over consensus peak file
bedtools coverage -a $conspeaks -b $bed > $samp.consensusPeaks.fragment.counts

	


