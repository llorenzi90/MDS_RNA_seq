#!/bin/bash
##SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH -p mem
##SBATCH --ntasks=16
#SBATCH --time=01:05:30

module load fastqc/0.11.8
samp=$1

fi=/scratch/llorenzi/MDS_RNA_seq_SEPT2021/$samp
echo $fi

	
cd $fi
mkdir fastqc
fastqc ${samp}_1.fq.gz -o fastqc
fastqc ${samp}_2.fq.gz -o fastqc

