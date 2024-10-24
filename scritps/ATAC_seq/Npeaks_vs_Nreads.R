dat=read.csv("/Users/llorenzi/Nextcloud/MDS/ATAC-seq/analyses/QC/per_sample_peak_info/per_sample_peak_info_-qval_9/all_samples_sample_peaks_info.-log10qval.9.csv")
tdat=t(dat)
colnames(tdat) <- tdat[1,]
tdat <- tdat[-1,]
tdat <- as.data.frame(tdat)

tdat$N_own_peaks=as.numeric(tdat$N_own_peaks)
tdat$Total_read_pairs_no_dups <- as.numeric(tdat$Total_read_pairs_no_dups)
cor(tdat$Total_read_pairs_no_dups,tdat$N_own_peaks)
pdf("Total_read_pairs_nodups_vs_N_own_peaks.pdf")
plot(tdat$Total_read_pairs_no_dups,tdat$N_own_peaks)
dev.off()

tdat$N_unique_peaks <- as.numeric(tdat$N_unique_peaks)
pdf("Total_read_pairs_nodups_vs_N_unique_peaks.pdf")
plot(tdat$Total_read_pairs_no_dups,tdat$N_unique_peaks)
dev.off()
