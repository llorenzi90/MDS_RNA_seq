TPM <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.gene_name.TPM.csv")
rownames(TPM)=TPM$gene_name
TPM <- TPM[,-1]
colnames(TPM) <- gsub("X","",colnames(TPM))



metadata <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv")

metadata$ATAC_seq="yes"
metadata$RNA_seq="no"
metadata$RNA_seq[metadata$patient.LPS%in%colnames(TPM)]="yes"
table(metadata$Cohesin[metadata$ATAC_seq=="yes"])
table(metadata$Cohesin[metadata$RNA_seq=="yes"])
write.csv(metadata,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv",row.names = F)
