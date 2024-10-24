#!/bin/bash


datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"
outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/average_plots"
sample_data_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"
#group_feat="Cohesin"

while read s 
do 
	TSS_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_files/$s"
	echo $TSS_file
	echo -e "Rscript '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/scripts/plot_midTSSsignal_meandensity.R' "${TSS_file}" "${datadir}" "${outdir}" "${sample_data_file}""
	Rscript '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/scripts/plot_midTSSsignal_meandensity.R' "${TSS_file}" "${datadir}" "${outdir}" "${sample_data_file}" 

done < '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/list_gene_sets_to_plot.txt'

