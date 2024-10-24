#!/bin/bash
##SBATCH -e %x.%j.err
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH -p mem #default is std
#SBATCH --ntasks=10
#SBATCH --time=05:05:30
##SBATCH --mem=16G

#set threads
PPN=10 #this has to match ntasks 

#command line arguments
sample=$1

#working dir
wdir="/scratch/llorenzi/MDS/ATAC-seq"

#indexes and ref files
bowtie2Index="/home/llorenzi/references/human/ncbi/bowtie2_index/GRCh38_noalt_as/GRCh38_noalt_as"
primary_chrs_bed="/home/llorenzi/references/human/ncbi/primary_chrs_nochrM.bed"

#load modules and set tool paths
bowtie2="/home/llorenzi/bowtie2-2.4.2-linux-x86_64/bowtie2"
picard="/home/llorenzi/picard.jar"
trim_galore="/home/llorenzi/TrimGalore-0.6.6/trim_galore"


echo $PWD
cd $wdir
#mkdir $samplename
cd $sample

read1=$(ls *_1.fq.gz)
read2=$(ls *_2.fq.gz)
echo "read1: $read1"
echo "read2: $read2"

#1) QC
echo `date`
echo -e "################################\n\tFASTQC started\n################################\n"
module load fastqc/0.11.8
mkdir fastqc
fastqc $read1 -o fastqc
fastqc $read2 -o fastqc
echo -e "################################\n\tFASTQC done\n################################\n"
echo `date  +'%r'`

#2)trim_galore
echo `date  +'%r'`
echo -e "################################\n\tTrimGalore started\n################################\n"

module load conda/current
conda activate cutadaptenv

$trim_galore --paired --length 35 --retain_unpaired --fastqc $read1 $read2

file=${read1%.fq.gz}
echo "file $file"
read1_val=${read1%.fq.gz}_val_1.fq.gz
echo "read1_val $read1_val"
read1_unpaired=${read1%.fq.gz}.unpaired_1.fq.gz
echo "read1_unpaired $read1_unpaired"
read2_val=${read2%.fq.gz}_val_2.fq.gz
echo "read2_val $read2_val"
read2_unpaired=${read2%.fq.gz}.unpaired_2.fq.gz
echo "read2_unpaired $read2_unpaired"

echo -e "################################\n\tTrimGalore done\n################################\n"
echo `date  +'%r'`

#3)alignment
echo `date  +'%r'`
echo -e "################################\n\tBowtie2 started\n################################\n"
$bowtie2 --very-sensitive -p 10 -x $bowtie2Index -X 2000 -N 1 -1 $read1_val -2 $read2_val -S $sample.bt2.sam
#bowtie2 --very-sensitive -x $bowtie2Index -X 2000 -N 1 -U $read1_unpaired -S $samplename.unpaired_1.bt2.sam
#bowtie2 --very-sensitive -x $bowtie2Index -X 2000 -N 1 -U $read2_unpaired -S $samplename.unpaired_2.bt2.sam
echo -e "################################\n\tBowtie2 done\n################################\n"
echo `date  +'%r'`

#4)Filtering steps and corresponding stats
echo `date  +'%r'`
echo -e "### $sample SAMTOOLS filter and stats including MT reads###" 
module load samtools/1.9 

#retain MAPQ 20 or more
samtools view -@ $PPN -b -h -q20 $sample.bt2.sam > $sample.quality.bam
#run samtools stats
samtools stats $sample.quality.bam > $sample.quality.stats
#retain only proper pairs, concordant pairs
samtools view -@ $PPN -b -f2 -h $sample.quality.bam > $sample.quality.concordant.bam
#run samtools stats
samtools stats $sample.quality.concordant.bam > $sample.quality.concordant.stats
#calculate fragment lenght distribution this time retaining chrM reads
samtools view $sample.quality.concordant.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $sample.chrM_incl.fragment_length_count.txt
echo -e "### $sample SAMTOOLS filter and stats including MT reads DONE###" 
echo `date  +'%r'`

#5) filter out mitochondrial reads
echo `date  +'%r'`
echo -e "### $sample SAMTOOLS remove MT reads and stats###" 
#retain only primary chromosomes
samtools view -@ $PPN -b -h -L $primary_chrs_bed $sample.bt2.sam > $sample.primary_chr.bam
samtools stats $sample.primary_chr.bam > $sample.primary_chr.stats
#retain MAPQ 20 or more
samtools view -@ $PPN -b -q20 -h $sample.primary_chr.bam > $sample.primary_chr.quality.bam
samtools stats $sample.primary_chr.quality.bam > $sample.primary_chr.quality.stats
#retain only proper pairs, concordant pairs
samtools view -@ $PPN -b -f2 -h $sample.primary_chr.quality.bam > $sample.primary_chr.quality.concordant.bam
samtools stats $sample.primary_chr.quality.concordant.bam > $sample.primary_chr.quality.concordant.stats
#calculate fragment lenght distribution this time excluding chrM reads
samtools view $sample.primary_chr.quality.concordant.stats | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $sample.chrM_excl.fragment_length_count.txt
echo -e "### $sample SAMTOOLS remove MT reads and stats DONE###" 
echo `date  +'%r'`
#retain chrM in another file (although it was not used later)
#samtools view -@ $PPN -b -h -L /home/llorenzi/references/ncbi/chrM.bed $sample.bt2.sam > $sample.chrM.bam
#samtools stats $sample.chrM.bam > $sample.chrM.stats

#6)samtools sort and index 
echo `date  +'%r'`
echo -e "### $sample SAMTOOLS sort and index primary_chr.quality.concordant.bam ###" 

samtools sort -@ $PPN -o $sample.primary_chr.quality.concordant.sorted.bam $sample.primary_chr.quality.concordant.bam
samtools index $sample.primary_chr.quality.concordant.sorted.bam
echo -e "### $sample SAMTOOLS sort and index primary_chr.quality.concordant.bam DONE###" 
echo `date  +'%r'`
#7) Mark duplicates
echo `date  +'%r'`
echo -e "### $sample picard Mark duplicates ###" 

#picard MarkDuplicates
PATH=/home/llorenzi/jdk-16.0.1/bin:$PATH
export PATH
java -jar $picard MarkDuplicates INPUT=$sample.primary_chr.quality.concordant.sorted.bam OUTPUT=$sample.primary_chr.quality.concordant.sorted.markedDups.bam METRICS_FILE=$sample.primary_chr.quality.concordant.sorted.dupmatrix TMP_DIR=/tmp/ CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=FALSE TAGGING_POLICY=All
#run samtools stats on marked dups
samtools stats $sample.primary_chr.quality.concordant.sorted.markedDups.bam > $sample.primary_chr.quality.concordant.sorted.markedDups.stats
echo -e "### $sample picard Mark duplicates and related stats done###" 
echo `date  +'%r'`

#8) convert bam files to bigWig CPM normalized bam_to_bigWig_CPM.sh 
effGenomeSize=2913022398 #from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# PPN=10
module load conda/current
conda activate deeptoolsenv
echo `date  +'%r'` 
echo -e "### $sample bamcoverage ###" 

bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --effectiveGenomeSize $effGenomeSize --bam $sample.primary_chr.quality.concordant.sorted.bam -o $sample.primary_chr.quality.concordant.sorted.CPM.bw
echo -e "### $sample bamcoverage ended ###" 
echo `date  +'%r'` 
