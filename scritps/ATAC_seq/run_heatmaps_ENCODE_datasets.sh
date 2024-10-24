#!/bin/bash

chain_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/liftover_chains/hg19ToHg38.over.chain"
datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"
outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/heatmaps/"
sample_data_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"
group_feat="Cohesin"

while read s 
do 
	peaks_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/ENCODE_datasets/$s"
	echo $peaks_file
	echo -e "Rscript '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/scripts/plot_midpeaksignal_heatmap_ChIP.R' "${peaks_file}" "${chain_file}" "${datadir}" "${outdir}" "${sample_data_file}" $group_feat"
	Rscript '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/scripts/plot_midpeaksignal_heatmap_ChIP.R' "${peaks_file}" "${chain_file}" "${datadir}" "${outdir}" "${sample_data_file}" $group_feat

done < '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/ENCODE_datasets/list_of_ENCODE_datasets.txt'
