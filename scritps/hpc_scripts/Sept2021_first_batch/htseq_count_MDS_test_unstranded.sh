#!/bin/bash

#SBATCH -e %x%j.err
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
#module load samtools/1.9 

#samtools sort -@ $PPN -n -o $sample.namesorted.bam $sample.hisat2.sam

module load conda/current
conda activate htseq

htseq-count -n $PPN --stranded=no --mode intersection-nonempty --supplementary-alignments ignore $sample.namesorted.bam $gtf

#submit: 
#PPN=8 (or another value)
#for f in /scratch/llorenzi/MDS_RNA_seq_SEPT2021/*;do s=`echo $f|cut -d "/" -f5`;echo $s; sbatch -J $s.htseq_MDS -n $PPN htseq_count_MDS_test_unstranded.sh $s $PPN; done
