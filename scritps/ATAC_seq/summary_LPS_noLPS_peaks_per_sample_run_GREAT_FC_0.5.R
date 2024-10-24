## ---------------------------
##
##
## Purpose of script: Run GREAT with peaks based on fold change
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-09-15
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

#load merged peaks data
merged_peaks <- fread("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/peak_files/merged_peaks.primary_chr.quality.concordant.sorted.markedDups.nodups_peaks.sorted.bed")
samples <- unique(grep(",",merged_peaks$V5,invert = T,value = T))
patient=sort(unique(gsub("(AT-MDS[0-9]*)(.*)","\\1",samples)))
#patient <- gsub("\\.","-",patient)
LPS=gsub("(AT-MDS[0-9]*.)(.*)","\\2",samples)
merged_peaks$id <- paste0(merged_peaks$V1,":",merged_peaks$V2,"-",merged_peaks$V3)



#add count data
datadir <- "~/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/count_matrices/merged_macs2peaks/"
countdata <- read.table(list.files(datadir,pattern = "counts",full.names = T),header = T)

rownames(countdata) <- countdata$id
countdata <- countdata[,-1]

covered_length <- read.table(list.files(datadir,pattern = "covered_peak_length",full.names = T),header = T)
rownames(covered_length) <- covered_length$id
#covered_length <- covered_length[,-1]

library(DESeq2)
name=colnames(countdata)
patient=gsub("(AT.MDS[0-9]*)(.*)","\\1",name)
#patient <- gsub("\\.","-",patient)
LPS=gsub("(AT.MDS[0-9]*.)(.*)","\\2",name)
coldata=data.frame(name,patient,LPS)
rownames(coldata) <- coldata$name
table(rownames(coldata)==colnames(countdata))

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~0+LPS + patient)

rl_data <- rlog(ddsMat, blind = TRUE, fitType = "parametric")
rlc <- assay(rl_data)
#generate matrix with FCs per patient
FCs <- sapply(unique(patient),function(x)rlc[,paste0(x,".2")] - rlc[,paste0(x,".1")])

head(FCs)
apply(FCs, 2, summary)
min(apply(FCs, 2, summary)[5,])
apply(FCs, 2, function(x)sum(x>=0.5))
#generate bedfiles with those peaks that have a rlog FC of 0.5 or more
FC=0.5
getwd()
setwd("~/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/LPS_vs_noLPS_analysis/")


setwd("patients_peaks_bed_files/")
dir.create(paste0("patients_peaks_LPSvsNOLPS_minlogFC_",FC))
setwd(paste0("patients_peaks_LPSvsNOLPS_minlogFC_",FC))
for(pat in patient){
  write.table(merged_peaks[merged_peaks$id%in%rownames(FCs)[FCs[,pat]>=FC],c(1:3,7)],
              paste0(pat,"_minlogFC_",FC,".txt"),sep = "\t",quote = F,row.names = F,col.names = F)
  
}

#Now upload these files to DropBox
#Get the dropbox links and copy them to file

#Run GREAT with all samples
#generate bed files with all patient's peaks that are exclusive for LPS+:

setwd("~/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/LPS_vs_noLPS_analysis/")

dplinks <- read.table("Dropbox_links_to_bed_files_FC_0.5.txt")
patient <- gsub("(.*)(AT.MDS[0-9]*)(.*)","\\2",dplinks$V1)
#pat=dplinks$V1[1]
#link=dplinks$V2[1]
cmm1 <- "wget -O "
cmm2 <- '"http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=hg38&requestName='
cmm3 <- "&requestURL="
comms <- c()
samp_counts <-0
for (i in 1:nrow(dplinks)) {
  samp_counts=samp_counts+1
  pat=patient[i]
  link <- dplinks$V1[i]
  comms <- c(comms,
             paste0(cmm1,
                    pat,
                    "_GREAT_batch.tsv ",
                    cmm2,
                    pat,
                    cmm3,
                    gsub("/","%2F",gsub(":","%3A",link)),
                    '"'))
  if(samp_counts==4){
    write(comms,paste0("run_GREAT_HTTP_request_samples_FC_0.5",i -4,"-",i,".sh"))
    samp_counts=0
    comms=c()
  }
}


#generate urls with requests to get web outputs

cmm4 <- 'http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?requestSpecies=hg38&requestName='
cmm5 <- "&requestURL="
comms <- c()
for (i in 1:nrow(dplinks)) {
  pat=patient[i]
  link <- dplinks$V1[i]
  comms <- c(comms,
             paste0(cmm4,
                    pat,
                    cmm5,
                    
                    gsub("/","%2F",gsub(":","%3A",link))))
  
}
write(comms,"HTTP_request_links_all_samples_FC_0.5.txt")

# test_res_from_web <- read.table("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/LPS_vs_noLPS_analysis/greatExportAll_AT-MDS35_from_web.tsv",fill = T,quote = "",sep = "\t")
# colnames(test_res_from_web) <- clns
# View(test_res_from_web[grep("GO Bio",test_res_from_web$Ontology),])
# 
# #check example result:
# #res <- read.table("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/LPS_vs_noLPS_analysis/AT-MDS35_GREAT_batch.tsv",fill = T,quote = "",sep = "\t",skip = 1)
# clns <- c("Ontology",
#           "ID",
#           "Desc",
#           "BinomRank",
#           "BinomP",
#           "BinomBonfP",
#           "BinomFdrQ",
#           "RegionFoldEnrich",
#           "ExpRegions",
#           "ObsRegions",
#           "GenomeFrac",
#           "SetCov",
#           "HyperRank",
#           "HyperP",
#           "HyperBonfP",
#           "HyperFdrQ",
#           "GeneFoldEnrich",
#           "ExpGenes",
#           "ObsGenes",
#           "TotalGenes",
#           "GeneSetCov",
#           "TermCov",
#           "Regions",
#           "Genes")
# colnames(res) <- clns
# View(res[grep("GO Bio",res$Ontology),])
# 
