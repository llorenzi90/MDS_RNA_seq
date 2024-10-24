## ---------------------------
##
##
## Purpose of script: Make complex OncoPrint plot
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-11-05
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

library(ComplexHeatmap)
## ---------------------------
#setwd("/home/llorenzi@CARRERASRESEARCH.ORG/share/Cuartero Group/CUARTERO GROUP/MDS/mutations/OncoPrint")
setwd("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/OncoPrint")

#note that all data is in sheet 3
raw_data=readxl::read_xlsx("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/mutations/Analisis_Pame/Resumen_analisis paneles Pame.xlsx",sheet = 3)

length(unique(raw_data$Sample))
length(unique(raw_data$Gene))
#for the type of plot I want to make, similar to the one in doi:10.1038/leu.2013.336, I don't need to use the type of mutation for now
#I just need the genes mutated per patient
mat=as.matrix(unclass(table(raw_data$Gene,raw_data$Sample)))

mat[mat!=0] <- "MUT"
mat[mat==0] <- "" 
mat <- mat[rownames(mat)!="SIN VARIANTES",]
mat <- mat[,colSums(mat!="")!=0]
#quick test
oncoPrint(mat)

#group genes
write.csv(data.frame(Gene=rownames(mat)),"mutated_genes.csv",row.names = F)
#add the groups based on supplementary table 3 from the paper: 
S3table <- read.csv("Suppl_TableS3_Leukemia_2014.csv")
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

oncoPrint(mat,row_split=grp)

# sort the groups by frequency of the most frequent gene 
Pathway_annot$freq <- rowSums(mat=="MUT")

max_freq_per_group=Pathway_annot %>%
  group_by(Pathway) %>%
  group_map(~ max(.x$freq))

max_freq_per_group <- unlist(max_freq_per_group)
names(max_freq_per_group) <- levels(grp)

grp <- factor(grp,levels = names(max_freq_per_group)[order(max_freq_per_group,decreasing = T)])


OP=oncoPrint(mat,row_split=grp, row_title_gp = gpar( fontsize = 12),row_title_rot = 0)
OP 


#test remove_empty_rows
matt=as.matrix(unclass(table(raw_data$Gene,raw_data$Sample)))

matt[mat!=0] <- "MUT"
matt[mat==0] <- "" 
#quick test
oncoPrint(matt,remove_empty_columns = T,remove_empty_rows = T,row_split=grp, row_title_gp = gpar( fontsize = 12),row_title_rot = 0)


##add clinical data as annotations

#test
oncoPrint(mat,row_split=grp, row_title_gp = gpar( fontsize = 12),row_title_rot = 0,
bottom_annotation = HeatmapAnnotation(
                                   foo1 = 1:ncol(mat),
                                   bar1 = anno_points(1:ncol(mat))
))

clinical_data <- readxl::read_xlsx("~/2104_gene_panels_FusterCuartero/Analisis_Pame/Datos clinicos_pacientes SC.xlsx")
matched_clinical_data <- clinical_data[match(colnames(mat),paste0("MDS",clinical_data$`Case ID`)),]

pdf("test2.pdf",)
oncoPrint(mat,row_split=grp, row_title_gp = gpar( fontsize = 8),row_title_rot = 0,
          bottom_annotation = HeatmapAnnotation(
            WHO = matched_clinical_data$`WHO 2017 subtype`,
            gender = matched_clinical_data$Gender,
            age = anno_points(matched_clinical_data$`Age at diagnosis`)
          ))  
dev.off()



pdf("test3.pdf",)
oncoPrint(mat,row_split=grp, row_title_gp = gpar( fontsize = 8),row_title_rot = 0,
          row_names_gp = gpar(fontsize = 8),
          pct_gp = gpar(fontsize = 8),
          bottom_annotation = HeatmapAnnotation(
            WHO = matched_clinical_data$`WHO 2017 subtype`,
            gender = matched_clinical_data$Gender,
            age = anno_points(matched_clinical_data$`Age at diagnosis`),annotation_name_gp =  gpar(fontsize = 8), annotation_legend_param = list( title_gp = gpar(fontsize = 8),labels_gp=gpar(fontsize=8))
          ),show_heatmap_legend = F)  
dev.off()


#font size:
fz=6
pdf("test3.pdf",)
oncoPrint(mat,row_split=grp, row_title_gp = gpar( fontsize = fz),row_title_rot = 0,
          row_names_gp = gpar(fontsize = fz),
          pct_gp = gpar(fontsize = fz),
          bottom_annotation = HeatmapAnnotation(
            WHO = matched_clinical_data$`WHO 2017 subtype`,
            gender = matched_clinical_data$Gender,
            age = anno_points(matched_clinical_data$`Age at diagnosis`),annotation_name_gp =  gpar(fontsize =fz), annotation_legend_param = list( title_gp = gpar(fontsize = fz),labels_gp=gpar(fontsize=fz))
          ),show_heatmap_legend = F)  
dev.off()


#size is fixed once text size is adjusted 

#Now control the colours
fz=6
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

alterfun= list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = cbp1[1], col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "black", col = NA))
  })

WHO_types=matched_clinical_data$`WHO 2017 subtype`[matched_clinical_data$`WHO 2017 subtype`!="NA"&!is.na(matched_clinical_data$`WHO 2017 subtype`)]
WHO_cols=colorRampPalette(colors = cbp1[-1])(length(unique(WHO_types)))
names(WHO_cols) <- unique(WHO_types)

pdf("test4.pdf")
oncoPrint(mat,
          alter_fun = alterfun,
          col=c(MUT="black"), row_split=grp, row_title_gp = gpar( fontsize = fz),row_title_rot = 0,
          row_names_gp = gpar(fontsize = fz),
          pct_gp = gpar(fontsize = fz),
          bottom_annotation = HeatmapAnnotation(
            WHO = matched_clinical_data$`WHO 2017 subtype`,
            gender = matched_clinical_data$Gender,
            age = anno_points(matched_clinical_data$`Age at diagnosis`), 
            col = list(gender = c(Male=cbp1[7],Female=cbp1[6]),
                       WHO = WHO_cols),na_col = "white",
            annotation_name_gp =  gpar(fontsize =fz), annotation_legend_param = list( title_gp = gpar(fontsize = fz),labels_gp=gpar(fontsize=fz))
          ),show_heatmap_legend = F)  
dev.off()
