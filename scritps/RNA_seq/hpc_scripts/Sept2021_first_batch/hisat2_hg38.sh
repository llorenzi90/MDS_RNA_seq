#!/bin/bash

##SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
#SBATCH --ntasks=16
#SBATCH --time=0-5
PPN=16

#usage: give complete path example: sbatch /scratch/llorenzi/MDS_RNA_seq_SEPT2021/1-2

module load hisat2/2.1.0

path=$1
sample=$(basename $path)

hisat2Index="/home/llorenzi/references/human/ncbi/hisat2_index/GRCh38_no_alt_as"

cd $path 

hisat2 -p $PPN -x $hisat2Index -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz -S $sample.hisat2.sam --novel-splicesite-outfile $sample.new_splice_sites.txt --rna-strandness RF --summary-file $sample.hisat2.summary --new-summary
