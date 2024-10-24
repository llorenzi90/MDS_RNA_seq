#!/bin/bash
##SBATCH -J 
##SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90àgmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH -p std
##SBATCH --ntasks=10
#SBATCH --time=10:05:30


echo $1
module load conda/current
conda activate macs2


cd /scratch/llorenzi/MDS/ATAC-seq/$1


macs2 callpeak -t $1.primary_chr.quality.concordant.sorted.markedDups.bam -n $1.primary_chr.quality.concordant.sorted.markedDups.nodups -f BAMPE -g hs --nomodel --nolambda --keep-dup 1 --call-summits --verbose 3

macs2 callpeak -t $1.primary_chr.quality.concordant.sorted.markedDups.bam -n $1.primary_chr.quality.concordant.sorted.markedDups.nodups -f BAMPE -g hs --nomodel --nolambda --keep-dup 1 --call-summits --verbose 3 --broad

#run: for samp in `ls /scratch/llorenzi/MDS/ATAC-seq/`; do sbatch -J $samp.macs2 MDS_macs2.sh $samp; done
