#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -p mem #default is std
#SBATCH --ntasks=10
#SBATCH --time=12:05:30
##SBATCH --mem=16G

#set threads

PPN=${2-10} #this has to match ntasks 

#command line arguments
sdir=$1
sample=$(basename $sdir) 


#indexes and ref files
bowtie2Index="/scratch/llorenzi/references/human/ncbi/bowtie2_index/GRCh38_noalt_as/GRCh38_noalt_as"
primary_chrs_bed="/scratch/llorenzi/references/human/ncbi/primary_chrs_nochrM.bed"

#load modules and set tool paths
bowtie2="/home/llorenzi/bowtie2-2.4.2-linux-x86_64/bowtie2"
picard="/home/llorenzi/picard.jar"
trim_galore="/home/llorenzi/TrimGalore-0.6.6/trim_galore"

cd $sdir

read1=$(ls *_1.fq.gz)
read2=$(ls *_2.fq.gz)
echo "read1: $read1"
echo "read2: $read2"

#1) QC
#echo `date`
#echo -e "################################\n\tFASTQC started\n################################\n"
#module load fastqc/0.11.8
#mkdir fastqc
#fastqc $read1 -o fastqc
#fastqc $read2 -o fastqc
#echo -e "################################\n\tFASTQC done\n################################\n"
echo `date  +'%r'`

#2)trim_galore
#echo `date  +'%r'`
#echo -e "################################\n\tTrimGalore started\n################################\n"

module load conda/current
#conda activate cutadaptenv

#$trim_galore --paired --length 35 --fastqc $read1 $read2


read1_val=$(ls *_val_1.fq.gz)
echo "read1_val $read1_val"

read2_val=$(ls *_val_2.fq.gz)
echo "read2_val $read2_val"

#echo -e "################################\n\tTrimGalore done\n################################\n"
#echo `date  +'%r'`

#3)alignment
echo `date  +'%r'`
echo -e "################################\n\tBowtie2 started\n################################\n"
$bowtie2 --very-sensitive -X 2000 -p $PPN -x $bowtie2Index -1 $read1_val -2 $read2_val -S $sample.bt2.sam
echo -e "################################\n\tBowtie2 done\n################################\n"
echo `date  +'%r'`

sam=$sample.bt2.sam

#4) convert to bam while sorting
echo `date  +'%r'`
echo -e "################################\n\tsamtools sort started\n################################\n"
module load samtools/1.9
samtools sort -@ $PPN -O bam -o $sample.sorted.bam $sam
echo -e "################################\n\tsamtools sort done\n################################\n"
echo `date  +'%r'`


#5) index
echo `date  +'%r'`
echo -e "################################\n\tsamtools index started\n################################\n"
samtools index -@ $PPN $sample.sorted.bam
echo -e "################################\n\tsamtools index done\n################################\n"
echo `date  +'%r'`

#6) mark duplicates
PATH=/home/llorenzi/jdk-16.0.1/bin:$PATH
export PATH

echo `date  +'%r'`
echo -e "################################\n\tpicard Mark duplicates started\n################################\n"
java -jar $picard MarkDuplicates INPUT=$sample.sorted.bam OUTPUT=$sample.sorted.markedDups.bam METRICS_FILE=$sample.dupmatrix TMP_DIR=/tmp/ CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE TAGGING_POLICY=All
echo -e "################################\n\tpicard Mark duplicates done\n################################\n"
echo `date  +'%r'`


#7) MAPQ distribution, raw and proper pairs, filtering and stats
echo `date  +'%r'`
echo -e "################################\n\tsamtools view and stats\n################################\n"

#A) FILTERS: primary chromosomes only (exclude chrM reads), proper pairs and MAPQ > 2
#A.1) exclude chrM reads
samtools view -@ $PPN -b -h -L $primary_chrs_bed $sample.sorted.markedDups.bam > $sample.sorted.markedDups.primary_chr.bam
#A.2) retain only proper pairs and reads with minimum MAPQ=2
samtools view -@ $PPN -b -h -f2 -q2 $sample.sorted.markedDups.primary_chr.bam > $sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam

#B) MAPQ distribution for raw reads and proper pairs with and without chrM
for fi in $sample.sorted.markedDups $sample.sorted.markedDups.primary_chr
do
	samtools view $fi.bam |cut -f5 |sort -n |uniq -c > $fi.MAPQdist
	samtools view -f2 $fi.bam |cut -f5 |sort -n |uniq -c > $fi.proper_pairs.MAPQdist
done

#C) STATS and fragment lenght distributions for all bam files (raw, no chrM and no chrM + qulity filter) 
for fi in $sample.sorted.markedDups $sample.sorted.markedDups.primary_chr $sample.sorted.markedDups.primary_chr.proper_pairs.minq2
do
	samtools stats $fi.bam > $fi.stats
	samtools view $fi.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $fi.fragment_length_count.txt
done

#D) index the final file
samtools index -@ $PPN $sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam


#E) remove intermediate sam and bam files (only interested in the filtered one for downstream analyses)
rm $sam $sample.sorted.bam $sample.sorted.bam.bai $sample.sorted.markedDups.bam $sample.sorted.markedDups.bai $sample.sorted.markedDups.primary_chr.bam

echo -e "################################\n\tsamtools view and stats done\n################################\n"
echo `date  +'%r'`


#8) bigwig
echo `date  +'%r'`
echo -e "################################\n\tbamcoverage started\n################################\n"
conda activate deeptoolsenv

bamCoverage --binSize 10 --normalizeUsing CPM --bam $sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam -o $sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam.CPM.bw
echo -e "################################\n\tbamcoverage done\n################################\n"
echo `date  +'%r'`


#run:
#cat MDS_samples.txt|while read sdir; do bn=$(basename $sdir); sbatch -J $bn.ATAC-seq_fromQC2bigWig ATAC_seq_qc_to_bigWig_human.sh $sdir
