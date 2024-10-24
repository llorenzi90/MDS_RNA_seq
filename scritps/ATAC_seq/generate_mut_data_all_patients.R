#setwd("/Users/llorenzi/Nextcloud/MDS/OncoPrint/Analisis_Pame/")
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/Analisis_Pame/")
#note that all data is in sheet 3
mut_data=readxl::read_xlsx("Resumen_analisis paneles Pame.xlsx",sheet = 3)
extra_samples=readxl::read_xlsx("Resumen_analisis paneles Pame.xlsx",sheet =4)
#extra_samples=extra_samples[c(1:3,8),]
MDS3=extra_samples[1:3,]

#for now I am only interested in genes mutated
MDS3[,!tolower(colnames(MDS3))%in%tolower(colnames(mut_data))] <- NA
matching_cols=match(tolower(colnames(mut_data)),tolower(colnames(MDS3)))
restof_cols=(1:ncol(MDS3))[!(1:ncol(MDS3))%in%matching_cols]
matching_cols[is.na(matching_cols)] <- restof_cols
MDS3 <- MDS3[,matching_cols]
colnames(MDS3) <- colnames(mut_data)



MDS35=extra_samples[7:8,]
colnames(MDS35)=MDS35[1,]
MDS35[,!tolower(colnames(MDS35))%in%tolower(colnames(mut_data))] <- NA
matching_cols=match(tolower(colnames(mut_data)),tolower(colnames(MDS35)))
restof_cols=(1:ncol(MDS35))[!(1:ncol(MDS35))%in%matching_cols]
matching_cols[is.na(matching_cols)] <- restof_cols
MDS35 <- MDS35[,matching_cols]
colnames(MDS35) <- colnames(mut_data)

mut_data=rbind(mut_data,MDS3,MDS35[-1,])




length(unique(mut_data$Sample))
length(unique(mut_data$Gene))
#I just need the genes mutated per patient
mat=as.matrix(unclass(table(mut_data$Gene,mut_data$Sample)))
#Classify genes by pathwaty using the Leukemia 2014 paper as reference
S3table <- read.csv("../OncoPrint/Suppl_TableS3_Leukemia_2014.csv")
table(rownames(mat)%in%S3table$Gene)
rownames(mat)[!rownames(mat)%in%S3table$Gene]
#there are 5 genes that are not in the paper supp table
# but U2AF1/1L5 are just U2AF1 in the paper
Pathway_annot <- data.frame(Gene=rownames(mat),Pathway=S3table$Pathway[match(rownames(mat),S3table$Gene)])
table(S3table$Pathway)
# add the missing annotations manually
#CSF3R is a receptor 
manual_annot <- c(CSF3R="Receptors/Kinases")
#CUX1 is a transcription factor "Transcription" group
manual_annot <- c(manual_annot,CUX1="Transcription")
#KMT2A is a methyl-tranferase
manual_annot <- c(manual_annot,KMT2A="Chromatin modification")
# STBP1 according to Wikipedia:
#The SETBP1 gene provides instructions for making a protein known as the SET binding protein 1, which is widely distributed throughout somatic cells. The protein is known to bind to another protein called SET. SETBP1 is a DNA-binding protein that forms part of a group of proteins that act together on histone methylation to make chromatin more accessible and regulate gene expression.[6] There is still more to learn about the overall function of the SETBP1 protein and the effect of SET binding.
#So it is a chromatin modifier
manual_annot <- c(manual_annot,SETBP1="Chromatin modification")
manual_annot <- c(manual_annot,`U2AF1;U2AF1L5`="RNA splicing")

Pathway_annot$Pathway[Pathway_annot$Gene%in%names(manual_annot)] <- manual_annot[match(Pathway_annot$Gene[Pathway_annot$Gene%in%names(manual_annot)],
                                                                                       names(manual_annot))]

grp=as.factor(Pathway_annot$Pathway)

mut_data$Pathway=Pathway_annot$Pathway[match(mut_data$Gene,Pathway_annot$Gene)]

#write.table(mut_data,"~/Nextcloud/MDS/mut_data_all_patients.txt",sep = "\t",row.names = F,quote = F)
write.table(mut_data,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/mut_data_all_patients.txt",sep = "\t",row.names = F,quote = F)
