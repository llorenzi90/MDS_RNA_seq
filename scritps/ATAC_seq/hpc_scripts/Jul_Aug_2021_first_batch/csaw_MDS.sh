#!/bin/bash
##SBATCH -J 
#SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p mem
#SBATCH --ntasks=10
#SBATCH --time=02:05:30

module load conda/current

conda activate csaw

#module load R/4.0.2

#Rscript /home/llorenzi/scripts/csaw.setup.MDS.R
#Rscript /home/llorenzi/scripts/csaw.edgeR.MDS.hpc.R
Rscript /home/llorenzi/scripts/csaw.filter.stat.MDS.hpc.R
