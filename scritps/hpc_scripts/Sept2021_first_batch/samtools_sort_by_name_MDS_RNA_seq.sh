#!/bin/bash

##SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
##SBATCH --ntasks=8
#SBATCH --time=0-5
PPN=$2

sample=$1

echo -e "$PPN\t$sample"

gtf="/home/llorenzi/references/human/annotation/gencode/gencode.v38.annotation.gtf"

cd /scratch/llorenzi/MDS_RNA_seq_SEPT2021/$sample 


#perform htseq count 
#first sort the sam file by name
module load samtools/1.9 

samtools sort -@ $PPN -n -o $sample.namesorted.bam $sample.hisat2.sam

#cat list_MDS_samples_failed.txt |while read samp; do sbatch -p mem -n 1 -J $samp.samtools_sortbyname samtools_sort_by_name_MDS_RNA_seq.sh $samp 1;done

