#!/bin/bash
#while read s
#do 
#	Rscript ~/shares/INVESTIGACIO/Cuartero\ Group/CUARTERO\ GROUP/MDS/ATAC-seq/scripts/plot_tiff_fromRDS.R $s 2>&1 | tee $s.plots.log
#done < ENCODEdatasets_to_plot.txt

#Run again ATF3 that failed for heatmaps (run only heatmaps)
#s=ATF3_K562_ENCFF950AHH.topallpeaks
#Rscript ~/shares/INVESTIGACIO/Cuartero\ Group/CUARTERO\ GROUP/MDS/ATAC-seq/scripts/plot_tiff_fromRDS.R $s hmap 2>&1 | tee $s.hmpas_only.log

outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/New_datasets/"
datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/RDS_files/New_datasets/"

cd "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks"
while read s
do 
	Rscript ~/shares/INVESTIGACIO/Cuartero\ Group/CUARTERO\ GROUP/MDS/ATAC-seq/scripts/plot_tiff_fromRDS.R $s -o "${outdir}" -d "${datadir}" 2>&1 | tee $s.New_dataset.plots.log
done < NewENCODEdatasets_to_plot.txt
