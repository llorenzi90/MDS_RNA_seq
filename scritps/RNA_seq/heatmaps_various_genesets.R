## ---------------------------
##
##
## Purpose of script: make heatmap of cohesin patients
##                    vs non-cohesin for interferon alpha genes
##
## Author: Lucia Lorenzi
##
## Date Created: 2022-03-11
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
library(ComplexHeatmap)

## ---------------------------
#
#read TPM data
TPM <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.gene_name.TPM.csv")
rownames(TPM)=TPM$gene_name
TPM <- TPM[,-1]
colnames(TPM) <- gsub("X","",colnames(TPM))

#read in metadata
metadata <- read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv")
#metadata <- cbind(metadata, unite(metadata[,2:3],sep = ".",col = "patient.LPS"))
#table(colnames(TPM)%in%metadata$patient.LPS)
#write.csv(metadata,"~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv",row.names = F)

#load DESeq2 results 
deseq2res=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/DESeq2/all_samples/cohesin_vs_nocohesin/all_samples_CohesinvsNoCohesinDESeq2_results.gene_names.csv")
padjco=0.05
deseq2res <- deseq2res%>%mutate(DIFF=ifelse(padj<=padjco&!is.na(padj),
                                            ifelse(log2FoldChange>0,"UP","DOWN"),
                                            "NO"))
source("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/scripts/cbPalette.R")

genesetfiles=list.files("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/TSS_files/",full.names = T)

for (gs in genesetfiles) {
  geneset=gsub("_"," ",gsub(".TSSs.txt","",basename(gs)))
  
  intfa <- read.table(gs,header = T)
  table(intfa$gene_name%in%rownames(TPM))
  
  data_to_plot <- TPM[rownames(TPM)%in%intfa$gene_name,]
  data_to_plot <- log2(data_to_plot+0.01)
  
  plot(density(as.matrix(as.data.frame(data_to_plot))))
  summary(as.vector(as.matrix(as.data.frame(data_to_plot))))
  #For now I have not implemented the automatization of this, but would like to do in the future
  # for now the range of values to vary the color scale is an input parameter:
  
  range_col_values=c(-6,10)
  cada=0.02
  mi=range_col_values[1]
  ma=range_col_values[2]
  
  #palette.breaks=seq(-6,-2,0.02)
  palette.breaks=seq(mi,ma,cada) #palette breaks for changing colors
  
  palette.breaks.tmp=seq(round(min(as.data.frame(data_to_plot))-1),
                         round(max(as.data.frame(data_to_plot))+1),0.02) #palette breaks for all values
  palette.breaks.tmp <-  round(palette.breaks.tmp,2)
  color.palette.tmp=vector(length = (length(palette.breaks.tmp)-1))
  
  #I will try 3 types of palette:
  palette_list <- list(white_blue_orange=c("white","#0072B2","#D55E00"),
                       blue_white_orange=c("#0072B2","white","#D55E00"),
                       white_blue=c("white","#0072B2"))
  
  
  #Now that we have our data ready to plot we need to adjust some plotting parameters:
  
  ##Generate column labels = -1.5, 1, 0.5 , 0 (TSS), 0.5, 1, 1.5 
  #points2label=c("-1","-0.5","0","0.5","1")
  #labcol[points2label] <- points2label
  
  #sort columns
  rownames(metadata) <- metadata$patient.LPS
  metadata$LPS[metadata$LPS==2] <- "yes"
  metadata$LPS[metadata$LPS==1] <- "no"
  
  data_to_plot <- data_to_plot[,order(colnames(data_to_plot)%in%metadata$patient.LPS[metadata$Cohesin])]
  colnames(data_to_plot)%in%metadata$patient.LPS[metadata$Cohesin==TRUE]
  annocoldf=metadata[match(colnames(data_to_plot),rownames(metadata)),c("Cohesin","LPS")]
  colnames(annocoldf) <- c("Cohesin_mutant", "LPS_treatment")
  
  ann_colors = list(Cohesin_mutant =c("TRUE"=cbPalette[4],"FALSE"=cbPalette[1]),
                    LPS_treatment = c(yes=cbPalette[7], no = cbPalette[3]))
  
  #without setting the palette:
  setwd("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/plots/")
  dir.create("heatmaps")
  setwd("heatmaps/")
  dir.create(gsub(" ","_",geneset))
  setwd(gsub(" ","_",geneset))
  
  #first save the  plot with default colors
  showleg=T 
  graphics.off()
  
  pdf("default_pheatmap_Colors.pdf")
  hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot),cluster_rows = T,
                           cluster_cols = F,
                           show_colnames = T,
                           fontsize_col = 4,
                           angle_col=c("45"),
                           annotation_col =  annocoldf,
                           annotation_colors = ann_colors,
                           show_rownames = F,
                           # color = color.palette.tmp,
                           # breaks = palette.breaks.tmp,
                           use_raster=F,
                           legend = showleg,
                           fontsize_number = 9, column_names_side = c("top"),
                           
                           heatmap_legend_param = list(title="log2(TPM+0.01)"),
                           main = geneset)
  print(hm)
  dev.off()
  
  for (pname in names(palette_list)) {
    color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
    #Set the palette and palette breaks so all plots have the same range of colors
    color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
    color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
    color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
    
    #################################HEATMAPS#####################
    #default scale:
    graphics.off()
    
    pdf(paste0("default_scaling_",pname,"_colors.pdf"))
    hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot),cluster_rows = T,
                                cluster_cols = F,
                                show_colnames = T,
                                fontsize_col = 4,
                                angle_col=c("45"),
                                annotation_col =  annocoldf,
                                annotation_colors = ann_colors,
                                show_rownames = F,
                                color = colorRampPalette(palette_list[[pname]])(100),
                                # breaks = palette.breaks.tmp,
                                use_raster=F,
                                legend = showleg,
                                fontsize_number = 9, column_names_side = c("top"),
                                heatmap_legend_param = list(title="log2(TPM+0.01)"),
                                main = geneset)
    print(hm)
    dev.off()
    
    #manual scale
    graphics.off()
    
    pdf(paste0("manual_scaling_",pname,"_colors.pdf"))
    
    hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot),cluster_rows = T,
                                cluster_cols = F,
                                show_colnames = T,
                                fontsize_col = 4,
                                angle_col=c("45"),
                                annotation_col =  annocoldf,
                                annotation_colors = ann_colors,
                                show_rownames = F,
                                color = color.palette.tmp,
                                breaks = palette.breaks.tmp,
                                use_raster=F,
                                legend = showleg,
                                fontsize_number = 9, column_names_side = c("top"),
                                heatmap_legend_param = list(title="log2(TPM+0.01)"),
                                main = geneset)
    print(hm)
    dev.off()
  } 
  dir.create("plots_with_row_annot")
  setwd("plots_with_row_annot/")
  ifrnagenes=deseq2res[match(rownames(data_to_plot),deseq2res$gene_name),]
  annorowdf <- data.frame(rank=ifrnagenes$stat,DIFF=ifrnagenes$DIFF)
  rownames(annorowdf) <- rownames(data_to_plot)
  #order must match the order of the input matrix:

  ann_colors$DIFF=c("UP"=cbPalette[7],"NO"=cbPalette[1],"DOWN"=cbPalette[6])
  ann_colors$rank=colorRampPalette(c("white",cbPalette[8]))(100)
  annorowdf$rank[is.na(annorowdf$rank)] <- 0
  for (pname in names(palette_list)) {
    color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
    #Set the palette and palette breaks so all plots have the same range of colors
    color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
    color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
    color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
    
    #################################HEATMAPS#####################
    #default scale:
    graphics.off()
    
    pdf(paste0("default_scaling_",pname,"_colors.pdf"))
    hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot),cluster_rows = T,
                                cluster_cols = F,
                                show_colnames = T,
                                fontsize_col = 4,
                                angle_col=c("45"),
                                annotation_col =  annocoldf,
                                annotation_row = annorowdf,
                                annotation_colors = ann_colors,
                                show_rownames = F,
                                color = colorRampPalette(palette_list[[pname]])(100),
                                # breaks = palette.breaks.tmp,
                                use_raster=F,
                                legend = showleg,
                                fontsize_number = 9, column_names_side = c("top"),
                                heatmap_legend_param = list(title="log2(TPM+0.01)"),
                                main = geneset)
    print(hm)
    dev.off()
    
    #manual scale
    graphics.off()
    
    pdf(paste0("manual_scaling_",pname,"_colors.pdf"))
    
    hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot),cluster_rows = T,
                                cluster_cols = F,
                                show_colnames = T,
                                fontsize_col = 4,
                                angle_col=c("45"),
                                annotation_col =  annocoldf,
                                annotation_row = annorowdf,
                                annotation_colors = ann_colors,
                                show_rownames = F,
                                color = color.palette.tmp,
                                breaks = palette.breaks.tmp,
                                use_raster=F,
                                legend = showleg,
                                fontsize_number = 9, column_names_side = c("top"),
                                heatmap_legend_param = list(title="log2(TPM+0.01)"),
                                main = geneset)
    print(hm)
    dev.off()
  } 
  
  
  
}

# #lo que continúa sin cuadrarme es que el GSEA de IFNa salía como muy upregulado en los cohesin mutants
# 18:20
# por el heatmap, uno diría lo contrario
# 18:24
# pero bueno, no tenemos que resolverlo ahora
# 18:24
# con el heatmap me vale para el report!

#let's check this


#volcano plot
library(ggplot2)

tt=ggplot(deseq2res, aes(x=log2FoldChange,y=-log2(padj))) +geom_point()+
  theme_classic() 
tt + geom_point(data =ifrnagenes ,
                aes(x=log2FoldChange,y=-log2(padj),col="red"))+
  geom_hline(yintercept = -log2(0.05))

View(deseq2res[match(gene_order_after_Clustering,deseq2res$gene_name),])

padjco=0.05


table(deseq2res$DIFF[deseq2res$gene_name%in%rownames(data_to_plot)])


#what if we plot the deseq2 normalized counts?
countdata=read.csv("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/RNA-seq/analyses/htseq_count_data/count_matrix/all_samples_gencodev38.counts.csv") 
rownames(countdata) <- rownames(TPM)

countdata <- countdata[,-1]
colnames(countdata) <- gsub("X","",colnames(countdata))
metadata <- metadata[match(colnames(countdata),metadata$patient.LPS),]
library(DESeq2)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,                                                                                   
                                 colData = metadata,                                                                                       
                                 design = ~Cohesin )  
keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]

ddsMat <- estimateSizeFactors(ddsMat)
dds <- DESeq(ddsMat)
resultsNames(dds)
res <- results(dds,alpha = 0.05)
resdf=as.data.frame(res)
normcounts <- counts(ddsMat,normalized=T)
log2_norm_counts <- log2(normcounts+0.01)

data_to_plot_log2norm <- log2_norm_counts[match(gene_order_after_Clustering,rownames(log2_norm_counts)),
                                          match(colnames(data_to_plot),colnames(log2_norm_counts))]

annorowdf <- annorowdf[match(rownames(data_to_plot_log2norm),rownames(annorowdf)),]
pdf(paste0("default_scaling_",pname,"_colors.log2DESeq2_normalized.pdf"))
hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot_log2norm), 
                            cluster_cols = F,
                            cluster_rows=F,
                            show_colnames = T,
                            fontsize_col = 4,
                            angle_col=c("45"),
                            annotation_col =  annocoldf,
                            annotation_row = annorowdf,
                            annotation_colors = ann_colors,
                            show_rownames = F,
                            color = colorRampPalette(palette_list[[pname]])(100),
                            # breaks = palette.breaks.tmp,
                            use_raster=F,
                            legend = showleg,
                            fontsize_number = 9, column_names_side = c("top"),
                            heatmap_legend_param = list(title="log2(deseq2norm+0.01)"),
                            main = geneset)
print(hm)
dev.off()

data_to_plot_mean <- aggregate(t(data_to_plot),by=list(Cohesin=metadata$Cohesin[match(colnames(data_to_plot),
                                                                          metadata$patient.LPS)]),mean)
rownames(data_to_plot_mean) <- data_to_plot_mean$Cohesin   

data_to_plot_mean <- t(data_to_plot_mean[,-1])
data_to_plot_mean <- data_to_plot_mean[match(gene_order_after_Clustering,rownames(data_to_plot_mean)),]


pdf(paste0("default_scaling_",pname,"_colors.meanlog2TPM.pdf"))

hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot_mean), 
                            cluster_cols = F,
                            cluster_rows=F,
                            show_colnames = T,
                            fontsize_col = 4,
                            angle_col=c("45"),
                            #annotation_col =  annocoldf,
                            annotation_row = annorowdf,
                            annotation_colors = ann_colors,
                            show_rownames = F,
                            color = colorRampPalette(palette_list[[pname]])(100),
                            # breaks = palette.breaks.tmp,
                            use_raster=F,
                            legend = showleg,
                            fontsize_number = 9, column_names_side = c("top"),
                            heatmap_legend_param = list(title="mean(log2(TPM+0.01))"),
                            main = geneset)
print(hm)
dev.off()


data_to_plot_deseq2mean <- aggregate(t(data_to_plot_log2norm),by=list(Cohesin=metadata$Cohesin[match(colnames(data_to_plot_log2norm),
                                                                                      metadata$patient.LPS)]),mean)
rownames(data_to_plot_deseq2mean) <- data_to_plot_deseq2mean$Cohesin   

data_to_plot_deseq2mean <- t(data_to_plot_deseq2mean[,-1])
data_to_plot_deseq2mean <- data_to_plot_deseq2mean[match(gene_order_after_Clustering,rownames(data_to_plot_deseq2mean)),]

pdf(paste0("default_scaling_",pname,"_colors.meanlog2deseq2norm.pdf"))
hm=ComplexHeatmap::pheatmap(as.matrix(data_to_plot_deseq2mean), 
                            cluster_cols = F,
                            cluster_rows=F,
                            show_colnames = T,
                            fontsize_col = 4,
                            angle_col=c("45"),
                            #annotation_col =  annocoldf,
                            annotation_row = annorowdf,
                            annotation_colors = ann_colors,
                            show_rownames = F,
                            color = colorRampPalette(palette_list[[pname]])(100),
                            # breaks = palette.breaks.tmp,
                            use_raster=F,
                            legend = showleg,
                            fontsize_number = 9, column_names_side = c("top"),
                            heatmap_legend_param = list(title="mean(log2(deseq2norm+0.01))"),
                            main = geneset)
print(hm)
dev.off()

