datadir="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/data/bigWig_files/"

bwfiles=list.files(datadir,full.names = T,pattern = "bam.CPM.bw")

mut_data=read.table("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/mut_data_all_patients.txt",sep = "\t",header = T)

samples=gsub(".sorted.markedDups.primary_chr.proper_pairs.minq2.bam.CPM.bw","",basename(bwfiles))
LPS=sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][2])
patient=gsub("AT-","",sapply(samples,function(x)strsplit(x,split = "\\.")[[1]][1]))
patient <- gsub("MDS","",patient)
Cohesin_patients=unique(mut_data$Sample[mut_data$Pathway=="Cohesin"])
Cohesin_patients <- Cohesin_patients[!is.na(Cohesin_patients)]
Cohesin=patient%in%Cohesin_patients
coldata=data.frame(samples,patient,LPS,Cohesin)                                                                                                   
#coldata$patient.LPS <- paste0(coldata$patient,coldate$LPS,sep=".")
write.csv(coldata,"/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv",quote = F,row.names = F)
