while read s; do unzip /scratch/llorenzi/MDS/ATAC-seq/$s/fastqc/*_1_fastqc.zip -d $s.fastqc; reads=$(grep 'Total Sequences' $s.fastqc/*/fastqc_data.txt |grep  -o -w .[0-9]*); echo -e "$s\t$reads" >> total_sequenced_reads.txt; rm -r $s.fastqc; done < ../../MDS_files.txt

while read s; do grep ^SN /scratch/llorenzi/MDS/ATAC-seq/$s/$s.chrM.stats| cut -f 2- > $s.chrM_SN.txt;done <../../MDS_files.txt
while read s; do grep ^SN /scratch/llorenzi/MDS/ATAC-seq/$s/$s.primary_chr.stats| cut -f 2- > $s.primary_chr_SN.txt;done <../../MDS_files.txt
while read s; do grep ^SN /scratch/llorenzi/MDS/ATAC-seq/$s/$s.primary_chr.quality.stats| cut -f 2- > $s.primary_chr.quality_SN.txt;done <../../MDS_files.txt
while read s; do grep ^SN /scratch/llorenzi/MDS/ATAC-seq/$s/$s.primary_chr.quality.concordant.stats| cut -f 2- > $s.primary_chr.quality.concordant_SN.txt;done <../../MDS_files.txt

while read s; do grep ^SN /scratch/llorenzi/MDS/ATAC-seq/$s/$s.primary_chr.quality.concordant.sorted.markedDups.stats| cut -f 2- > $s.primary_chr.quality.concordant.sorted.markedDups_SN.txt;done <../../MDS_files.txt

while read s; do grep ^IS /scratch/llorenzi/MDS/ATAC-seq/$s/$s.primary_chr.quality.concordant.sorted.markedDups.stats| cut -f 2- > $s.primary_chr.quality.concordant.sorted.markedDups_insert_sizes.txt;done <../../MDS_files.txt



