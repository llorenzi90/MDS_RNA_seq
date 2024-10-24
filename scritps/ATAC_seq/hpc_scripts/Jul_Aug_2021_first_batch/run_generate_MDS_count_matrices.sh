#!/bin/bash

#SBATCH -J generate_MDS_count_matrices
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p std
##SBATCH --ntasks=8
#SBATCH --time=0-2

module load R/3.6.3

Rscript /home/llorenzi/scripts/generate_count_covlength_matrices.R
