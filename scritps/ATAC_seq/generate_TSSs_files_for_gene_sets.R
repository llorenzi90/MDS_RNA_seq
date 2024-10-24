library(rtracklayer)
# gtf=readGFF("/home/llorenzi/inflammatory_genes_gencode.v38.annotation.genes.gtf")
gtf <- readGFF("/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/annotation/gencode/gencode.v38.annotation.gtf")
gene_sets_folder <- "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets/"
gene_sets=list.files(gene_sets_folder,pattern = "HALLMARK",full.names = T)

gene_sets
# [1] "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets//HALLMARK_INFLAMMATORY_RESPONSE.txt"    
# [2] "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets//HALLMARK_INTERFERON_ALPHA_RESPONSE.txt"
# [3] "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets//HALLMARK_INTERFERON_GAMMA_RESPONSE.txt"
# [4] "/home/llorenzi/shares/INVESTIGACIO/Cuartero Group/CUARTERO GROUP/references/human/gene_sets//HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt" 
outdir=paste0(gene_sets_folder,"/TSS_files")
setwd(outdir)
for (gs in gene_sets) {
  tgenes <- read.table(gs,skip=2)
  setname <- gsub(".txt","",basename(gs))
  print(table(tgenes$V1%in%gtf$gene_name))
  tgtf <- gtf[gtf$type=="gene"&gtf$gene_name%in%tgenes$V1,c("seqid","start","gene_name")]
  tgtf <- tgtf[match(tgenes$V1,tgtf$gene_name),]
  write.table(tgtf,paste0(setname,".TSSs.txt"),row.names = F,quote = F,sep = "\t")
}

