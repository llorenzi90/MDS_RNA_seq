cArgs=commandArgs(trailingOnly = T)
samp=cArgs[1]
if(length(cArgs)==1) what_to_plot="both" else{
  if(cArgs[2]==1|cArgs[2]=="avg") what_to_plot="avg" else if(cArgs[2]==2|cArgs[2]=="hmap"){
    what_to_plot="hmap"
  } else what_to_plot="both"

}
if(any(cArgs=="-o")|any(cArgs=="-outdir")){
  outdir=cArgs[which(cArgs%in%c("-o","-outdir"))+1]
}else outdir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/"

if(any(cArgs=="-d")|any(cArgs=="-datadir")){
  datadir=cArgs[which(cArgs%in%c("-d","-datadir"))+1]
}else datadir="/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/meanCoverage_aroundENCODEPeaks/RDS_files/"



print(paste0(what_to_plot," plots will be done"))
#datadir=cArgs[2]
#outdir=cArgs[3]
#sample_data_file=cArgs[4]
#group_feat=cArgs[5]

group_feat="Cohesin"
outdir_hm <- paste0(outdir,"/heatmaps")
outdir_ap <- paste0(outdir,"/average_plots")
sample_data_file="~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/processed_samples.metadata.csv"
sample_table <- read.csv(sample_data_file) #this have to have a column "samples"

#samp="STAT1_K562_ENCFF884LDR.topallpeaks"

sort_by="CPM"
region_to_plot=c(-1000,1000)
##
#lower limit
llimit=region_to_plot[1]
#upper limit
ulimit=region_to_plot[2]
#region width
rwidth=length(llimit:ulimit)

#Color scale parameters
range_col_values=c(-6,0)
cada=0.02
mi=range_col_values[1]
ma=range_col_values[2]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if(what_to_plot=="both"|what_to_plot=="avg"){
  ######################PLOT average plots######################
  #####1) Plot density for each sample showing each group in a different color
  #define color palette (color-blind friendly)
  #define function to plot
  plot_densities <- function(dfin,lwd=2,main="",
                             xlab="distance from peak center",
                             ylab="mean CPM",xvals=llimit:ulimit,
                             cols=cbPalette[6:7],
                             sampleTable=sample_table,
                             groupCol=group_feat){
    #Define sample sets: has to be a factor of two groups:
    #group_feat="LPS" #to test
    sampleTable[,groupCol] <- as.factor(sampleTable[,groupCol])
    sgr <- paste0(groupCol,"_",levels(sampleTable[,groupCol]))
    vcols <- cols[as.numeric(as.factor(sampleTable[,groupCol]))]
    
    #we use it to generate average vectors for each position*peak
    #across all samples from each group
    
    maxyval=max(apply(dfin, 2,max))
    minyval=min(apply(dfin, 2,min))
    plot(xvals, dfin[,1],col=vcols[1],
         type="s",ylim=c(minyval,maxyval+0.02),ylab=ylab,
         xlab=xlab,lwd=lwd,main=main
    )
    
    legend("topright",lwd = lwd,legend = sgr,
           title = groupCol,bty = "n",col = cols)
    for (sa in 2:ncol(dfin)) {
      points(xvals, dfin[,sa],
             type="s",col=vcols[sa],lwd=2)
    }
    
  }
  
  #test:
  #plot_densities(mean_posCPM,sampleTable = sample_table[1:10,],groupCol = "LPS",main="LPS_vs_noLPS")
  mean_posCPM_file=paste0(samp,".mean_sampleCPM.rds")
  mean_posCPM <- readRDS(paste0(datadir,mean_posCPM_file))
  
  if(!dir.exists(outdir_ap))dir.create(outdir_ap,recursive = T)
  setwd(outdir_ap)
  plot_base_name <- paste0(samp,".meanCoverage")
  
  for (gf in c("Cohesin","LPS")) {
    comp <- paste(paste0(gf,levels(as.factor(sample_table[,gf]))),
                  collapse = "vs")
    
    #A) plot all samples individually
    graphics.off()
    pdf(paste("all_samples",comp,plot_base_name,"pdf",sep = "."))
    plot_densities(mean_posCPM,sampleTable = sample_table,groupCol = gf)
    dev.off()
    
    #B) compute average for each group
    mean_posCPM_grpAvg <- t(aggregate(t(mean_posCPM),by=list(grp=sample_table[,gf]),mean)[,-1])
    sample_tablex <- data.frame(levels(as.factor(sample_table[,gf])))
    colnames(sample_tablex) <- gf          
    graphics.off()
    pdf(paste("Avg_samples",comp,plot_base_name,"pdf",sep = "."))
    plot_densities(mean_posCPM_grpAvg,sampleTable = sample_tablex,groupCol = gf,ylab = "mean CPM (samples avg)")
    dev.off()
    
  }
  
}


if(what_to_plot=="both"|what_to_plot=="hmap"){
  ######################PLOT heatmaps######################
  library(ComplexHeatmap)
  
  data_to_plot_file=paste0("Cohesin_FALSE_vs_Cohesin_TRUE.",samp,".log2.datatoplot.rds")
  data_to_plot <- readRDS(paste0(datadir,data_to_plot_file))
  
  #2)   Create palette to apply to both plots using the absolute minimum and max values across both matrices
  #a) vectorize both matrices and check values. This part is a bit manual
  # The idea is to select a range of values where most of the CPMs fall into
  # and use them to scale the colors of the palette
  
  #plot(density(as.matrix(as.data.frame(data_to_plot))))
  #summary(as.vector(as.matrix(as.data.frame(data_to_plot))))
  #For now I have not implemented the automatization of this, but would like to do in the future
  # for now the range of values to vary the color scale is an input parameter
  
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
  labcol <- rep("",rwidth)
  names(labcol)=as.character((llimit:ulimit)/1000)
  labcol[1+c(0,1,2,3,4)*round(rwidth/4)] <- names(labcol[1+c(0,1,2,3,4)*round(rwidth/4)])
  #points2label=c("-1","-0.5","0","0.5","1")
  #labcol[points2label] <- points2label
  
  if(!dir.exists(outdir_hm)) dir.create(outdir_hm,recursive = T)
  setwd(outdir_hm)
  
  plot_base_name <- paste0(paste(names(data_to_plot),collapse = "_vs_"),".",samp,".log2.heatmap")
  
  for (pname in names(palette_list)) {
    color.palette  <- colorRampPalette(palette_list[[pname]])(length(palette.breaks) - 1)
    #Set the palette and palette breaks so all plots have the same range of colors
    color.palette.tmp[(which(palette.breaks.tmp==as.character(mi)):which(palette.breaks.tmp==as.character(ma - cada)))]=color.palette
    color.palette.tmp[1:(which(palette.breaks.tmp==as.character(mi))-1)]=color.palette.tmp[which(palette.breaks.tmp==as.character(mi))]
    color.palette.tmp[which(palette.breaks.tmp==as.character(ma)):length(color.palette.tmp)]=color.palette.tmp[which(palette.breaks.tmp==as.character(ma -cada))]
    
    #################################HEATMAPS#####################
    heats <- lapply(1:2,function(n){
      if(n==1) showleg=T else showleg=F
      m=data_to_plot[[n]]
      return(ComplexHeatmap::pheatmap(m,cluster_rows = F,
                                      cluster_cols = F,
                                      show_colnames = T,
                                      labels_col = labcol,
                                      show_rownames = F,
                                      color = color.palette.tmp,
                                      breaks = palette.breaks.tmp,
                                      use_raster=F,
                                      legend = showleg,
                                      fontsize_number = 9, column_names_side = c("top"),
                                      angle_col = c("0"),
                                      heatmap_legend_param = list(title="log2(CPM+0.01)"),
                                      main = names(data_to_plot)[n]))
    } )
    
    #graphics.off()
    
    # pdf(paste0(plot_base_name,".colorscale_",pname,".pdf"))
    # print(heats[[1]] + heats[[2]])
    # dev.off()
    #
    graphics.off()
    tiff(paste0(plot_base_name,".colorscale_",pname,".tiff"), width = 8, height = 8, units = 'in', res = 300)
    print(heats[[1]] + heats[[2]])
    dev.off()
    
    
  }
  
  
  
}
