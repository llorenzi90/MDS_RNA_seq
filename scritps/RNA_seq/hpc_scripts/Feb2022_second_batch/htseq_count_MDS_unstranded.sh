#!/bin/bash

#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem
##SBATCH --ntasks=8
#SBATCH --time=0-5
PPN=8

sample=${1:-"34-1"}

echo -e "$PPN\t$sample"

gtf="/scratch/llorenzi/references/human/annotation/gencode/gencode.v38.annotation.gtf"

cd /scratch/llorenzi/MDS_RNA_seq_Feb2022/$sample 


# #perform htseq count 
# #first sort the sam file by name
# module load samtools/1.9 
# echo -e "################################\n\tsamtools sort started\n################################\n"
# 
# samtools sort -@ $PPN -n -o $sample.namesorted.bam $sample.hisat2.sam
# echo -e "################################\n\tsamtools sort done\n################################\n"

module load conda/current
conda activate htseq
echo -e "################################\n\thtseq-count started\n################################\n"
htseq-count -n $PPN --stranded=no --mode intersection-nonempty --supplementary-alignments ignore $sample.namesorted.bam $gtf > $sample.htseq.count
echo -e "################################\n\thtseq-count done\n################################\n"
#submit: 
#PPN=8 (or another value)
#tail -n +2 list_of_files.txt|while read s; do sbatch -J $s.htseq-count htseq_count_MDS_unstranded.sh $s; done
