# MDS_RNA_seq
## Gene expression and chromatin accessibility in response to inflammation of Myelodysplastic Syndromes (MDS) patient samples.

 

While it is clear that inflammatory deregulation is a hallmark of Myelodysplastic Syndromes (MDS), it is yet unclear how MDS mutations promote changes in inflammatory response or inflammatory gene expression. Knowing it would be a huge step towards understanding the biology of this disease. Moreover, it could help to identify specific molecular signatures for different mutations and thus help guide the startification of patients for treatment choice.  

Here, we will analyze gene expression (RNA-seq), chromatin accessibility (ATAC-seq), mutations and clinical data of ~40 MDS bone marrow samples. For RNA-seq and ATAC-seq, each sample has been treated with LPS to examine their capacity of mounting an inflammatory response (so 2 libraries/patient=control and LPS). Accessibility data will be used to understand the enhancers and gene regulatory modules that are associated to each mutation,  while gene expression data is used to understand the genes that are expressed by each mutant.  

 

 

 

## Targeted DNA Analysis 

Identify mutations in all patients 

Perform Oncoprint plots

Heatmap of significant overlap and exclusion between mutations 

Investigate clinical features of cohesin mutations and the rest: are ther any significant difference? 

 

 

 

## ATAC-seq analysis 

*samples ending in .1 are non-treated, and .2 are LPS-treated 

Basic processing 

QC 

Trim adapters 

Alignment 

Fragment size distribution; % of mapped reads; % of duplicates 

BIGWIG files 

Peak calling (MACS2) 

Summary of number of peaks per sample and FRIP 

BED files 

PCA of LPS and PBS treated samples (separately and together) 

Heatmap of sample clustering: again, LPS and PBS separately and together 

 

Identify gained peaks in LPS-treated samples 

PCA 

Peak annotation 

Motif enrichment 

 
