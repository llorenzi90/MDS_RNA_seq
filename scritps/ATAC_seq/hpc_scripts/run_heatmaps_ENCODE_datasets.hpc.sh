#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem #default is std
##SBATCH --ntasks=10
#SBATCH --time=3:35:30
##SBATCH --mem=16G

s=$1

peaks_file="/scratch/llorenzi/ENCODE_datasets/$s"
chain_file="/scratch/llorenzi/references/human/liftover_chains/hg19ToHg38.over.chain"
datadir="/scratch/llorenzi/MDS_ATAC_seq_bigWig_files"
outdir="/scratch/llorenzi/MDS_ATAC_seq_plots/meanCoverage_aroundENCODEPeaks/heatmaps/"
sample_data_file="/scratch/llorenzi/MDS_ATAC_seq_plots/processed_samples.metadata.csv"
group_feat="Cohesin"
Ntop="all"


module load conda/current
conda activate atacseqqc

echo $peaks_file
echo -e "Rscript '/home/llorenzi/scripts/plot_midpeaksignal_heatmap_ChIP.hpc.R' "${peaks_file}" "${chain_file}" "${datadir}" "${outdir}" "${sample_data_file}" $group_feat $Ntop"
Rscript '/home/llorenzi/scripts/plot_midpeaksignal_heatmap_ChIP.hpc.R' "${peaks_file}" "${chain_file}" "${datadir}" "${outdir}" "${sample_data_file}" $group_feat $Ntop


