#!/bin/bash
##SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH -p mem #default is std
##SBATCH --ntasks=10
#SBATCH --time=00:35:30
##SBATCH --mem=16G

#set threads
#PPN=10 #this has to match ntasks 

sample=$1

cd $wdir
#mkdir $samplename
cd $sample

#5) filter out mitochondrial reads
echo `date  +'%r'`
echo -e "### $sample SAMTOOLS no MT fragment length extraction###" 
#calculate fragment lenght distribution this time excluding chrM reads
samtools view $sample.primary_chr.quality.concordant.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $sample.chrM_excl.fragment_length_count.txt
echo -e "### $sample SAMTOOLS fragment length DONE###" 
echo `date  +'%r'`
