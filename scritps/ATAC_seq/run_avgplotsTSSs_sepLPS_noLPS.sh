#!/bin/bash


datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"
outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundTSSs/average_plots"
sample_data_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"
#group_feat="Cohesin"

for s in HALLMARK_INFLAMMATORY_RESPONSE.TSSs.txt HALLMARK_INTERFERON_GAMMA_RESPONSE.TSSs.txt HALLMARK_INTERFERON_ALPHA_RESPONSE.TSSs.txt HALLMARK_TNFA_SIGNALING_VIA_NFKB.TSSs.txt
do 
	TSS_file="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_files/$s"
	echo $TSS_file
	echo -e "Rscript '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/scripts/plot_midTSSsignal_meandensity_sepLPS_noLPS.R' "${TSS_file}" "${datadir}" "${outdir}" "${sample_data_file}""
	Rscript '/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/scripts/plot_midTSSsignal_meandensity_sepLPS_noLPS.R' "${TSS_file}" "${datadir}" "${outdir}" "${sample_data_file}" 

done 

#~/shares/INVESTIGACIO/Cuartero\ Group/CUARTERO\ GROUP/MDS/ATAC-seq/scripts/run_avgplotsTSSs_sepLPS_noLPS.sh 2>&1 | tee Avgplots_sepLPS_noLPS.log
