#!/bin/bash

#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem
##SBATCH --ntasks=8
#SBATCH --time=0-5
PPN=$2

sample=${1:-"34-1"}

echo -e "$PPN\t$sample"

gtf="/scratch/llorenzi/references/human/annotation/gencode/gencode.v38.annotation.gtf"

cd /scratch/llorenzi/MDS_RNA_seq_Feb2022/$sample 


#perform htseq count 
#first sort the sam file by name
module load samtools/1.9 
echo -e "################################\n\tsamtools sort started\n################################\n"

samtools sort -@ $PPN -n -o $sample.namesorted.bam $sample.hisat2.sam
echo -e "################################\n\tsamtools sort done\n################################\n"

module load conda/current
conda activate htseq
echo -e "################################\n\thtseq-count started\n################################\n"
htseq-count -n $PPN --stranded=reverse --mode intersection-nonempty --supplementary-alignments ignore $sample.namesorted.bam $gtf > $sample.htseq.count
htseq-count -n $PPN --stranded=no --mode intersection-nonempty --supplementary-alignments ignore $sample.namesorted.bam $gtf > $sample.htseq.unstranded.count
echo -e "################################\n\thtseq-count done\n################################\n"
#submit: 
#PPN=8 (or another value)
#for f in /scratch/llorenzi/MDS_RNA_seq_SEPT2021/*;do s=`echo $f|cut -d "/" -f5`;echo $s; sbatch -J $s.htseq_MDS -n $PPN htseq_count_MDS.sh $s $PPN; done
