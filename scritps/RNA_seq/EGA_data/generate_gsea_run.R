## ---------------------------
##
##
## Purpose of script: generate bash to run GSEA MDS
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-15
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------
source("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/get_relative_path_Rfunction.R")

#command example: gsea-cli.sh GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/gene_sets/h.all.v7.3.symbols.gmt -collapse Remap_Only -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/llorenzi/Nextcloud/TCGA/LL_analyses/TCGA_data_analysis/results/limma_voom/rank_files/ASXL1_limma_voom.rnk -scoring_scheme weighted -rpt_label my_analysis -chip /Users/llorenzi/Nextcloud/TCGA/LL_analyses/TCGA_data_analysis/data/GSEA/match_gene_sets.chip -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /Users/llorenzi/gsea_home/output/mar31
#cmm1 <- "./gsea-cli.sh GSEAPreranked -gmx /Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results/Mm.h.all.v7.1.entrez.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk /Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results/rank_files/"
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/")
dir.create("GSEA")
setwd("GSEA")
gsea_rel_path=get_relative_path("/home/llorenzi/software/GSEA_4.2.2/")

rank_file='/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/DESEq2/CohesinMut_vs_nonCohesinMut/rank_files/DESeq2.rnk'
relrank=get_relative_path(rank_file)


cmm1 <- paste0(gsea_rel_path,"/gsea-cli.sh GSEAPreranked -gmx ")
               
cmm2 <- " -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk "
cmm3 <- " -scoring_scheme weighted -rpt_label "
cmm4 <- " -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out "

commands <- c()
#for(f in list.files("/Users/llorenzi/Nextcloud/CBPa_RNA_seq/DESeq2_results/rank_files/")){
genesets=list.files("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/",pattern = ".gmt",full.names = T)
genesets <- genesets[c(4,1,2,3)]

nam="EGA_data"
for(gset in genesets){
  #rank_file='/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/DESeq2/all_samples/cohesin_vs_nocohesin/rank_files/all_samples_CohesinvsNoCohesin_DESeq2.rnk'
  #relrank=get_relative_path('/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/DESeq2/all_samples/cohesin_vs_nocohesin/rank_files/all_samples_CohesinvsNoCohesin_DESeq2.rnk')
  dir.create(basename(gset))
  relgsetpath <- get_relative_path(gset)
  #nam <- gsub("_DESeq2.rnk","",basename(relrank))
  
  commands <- c(commands, paste0(cmm1,relgsetpath,
                                 cmm2,relrank,
                                 cmm3, nam,
                                 cmm4, basename(gset)))
}

#write(commands,"/Users/llorenzi/Nextcloud/CBPa_RNA_seq/run_GseaPreranked_CBPa_DESeq2.sh")
write(commands,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/scripts/run_GseaPreranked_EGA_data_DESeq2.sh")
getwd()
#run this script from "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA_seq_EGA_blood2017/GSEA"

