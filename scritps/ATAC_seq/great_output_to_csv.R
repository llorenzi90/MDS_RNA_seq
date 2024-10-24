tsv_files <- list.files("~/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/MDS/ATAC-seq/analyses/GREAT/outputs",
                        full.names = T)


for (tsv in tsv_files) {
  da <- read.delim(tsv,skip=3)
  colnames(da)[1]="Ontology"
  write.csv(da,gsub(".tsv",".csv",tsv),row.names = F)
}
