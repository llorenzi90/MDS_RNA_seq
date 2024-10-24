#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem #default is std
#SBATCH --ntasks=10
#SBATCH --time=00:35:30
##SBATCH --mem=16G

myd=/scratch/llorenzi/mytest/MDS25.1
mkdir -p $myd

cp /scratch/llorenzi/MDS_new_ATAC_seq/MDS25.1/MDS25.1.sorted.markedDups.primary_chr.proper_pairs.minq2.bam $myd
cp /scratch/llorenzi/MDS_new_ATAC_seq/MDS25.1/MDS25.1.sorted.markedDups.primary_chr.proper_pairs.minq2.bam.bai $myd

/home/llorenzi/scripts/ATAC_seq_by_steps.sh -ppn 10 -rs bigWig -only $myd
