## ---------------------------
##
##
## Purpose of script: generate OncoPrint plot for MDS patients
##
## Author: Lucia Lorenzi
##
## Date Created: 2021-11-03
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

raw_data=readxl::read_xlsx("~/2104_gene_panels_FusterCuartero/Analisis_Pame/Resumen_analisis paneles Pame.xlsx")

length(unique(raw_data$Sample))
length(unique(raw_data$Gene))
table(raw_data$`Exonic mutation type`)
colnames(raw_data)
table(raw_data$`Mutation type`)
raw_data$mutation_type=raw_data$`Exonic mutation type`
raw_data$mutation_type[raw_data$`Exonic mutation type`=="NA"] <- "splicing"
table(raw_data$mutation_type)
table(raw_data$Sample,raw_data$Gene,raw_data$mutation_type)
test=table(raw_data$Sample,raw_data$Gene,raw_data$mutation_type)
length(test)

mutations <- list()
for (i in names(table(raw_data$mutation_type))) {
  mutations[[i]]=t(test[,,i])
}
names(mutations) <- c("frameshift_del",
                      "nonframeshit_del",
                      "nonsynonSNV",
                      "splicing",
                      "stopgain")

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

alter_fun_list = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = cbp1[1], col = NA))
  },
  frameshift_del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = cbp1[2], col = NA))
  },
  nonframeshit_del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = cbp1[3], col = NA))
  },
  nonsynonSNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = cbp1[4], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, just = "bottom",gp = gpar(fill = cbp1[5], col = NA))
  },
  stopgain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),just = "top", h*0.33, gp = gpar(fill = cbp1[6], col = NA))
  }
)
col = c("frameshift_del"=cbp1[2],
        "nonframeshit_del"=cbp1[3],
        "nonsynonSNV"=cbp1[4],
        "splicing"=cbp1[5],
        "stopgain"=cbp1[6])
oncoPrint(mutations,alter_fun = alter_fun_list,col = col)

pdf("MDS_mutations_OncoPrint.pdf")
print(oncoPrint(mutations,alter_fun = alter_fun_list,col = col))
dev.off()
