#!/bin/bash
##SBATCH -J 
##SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mem=32G
##SBATCH -p std
##SBATCH --ntasks=10
#SBATCH --time=01:25:30

while read s
do
	cut -f1-4,9,6 $s >> /scratch/llorenzi/MDS/ATAC-seq/merged_peaks/all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak 
done < narrowPeak_files_paths.txt



#module load samtools/1.9
module load bedtools2/2.29.0

cd /scratch/llorenzi/MDS/ATAC-seq/merged_peaks

sort -k1,1 -k2,2n all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.narrowPeak > all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.narrowPeak
sed "s/.primary_chr.quality.concordant.sorted.markedDups.nodups_peak_[0-9a-z]*//"  all_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.narrowPeak | bedtools merge -c 4,4,6 -o count_distinct,distinct,mean -i - > merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed


