#!/bin/bash
#SBATCH -o %x.%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=10
#SBATCH --time=05:05:30
PPN=10

path=$1
sample=$(basename $path)
cd $path

#1) sort

echo `date  +'%r'`
echo -e "################################\n\tsamtools sort started\n################################\n"
module load samtools/1.9
samtools sort -@ $PPN -o $sample.hisat2.sorted.bam $sample.hisat2.sam
sortedbam="$sample.hisat2.sorted.bam"
echo -e "################################\n\tsamtools sort done\n################################\n"
echo `date  +'%r'`


#2) index

echo `date  +'%r'`
echo -e "################################\n\tsamtools index started\n################################\n"
samtools index $sortedbam
echo -e "################################\n\tsamtools index done\n################################\n"
echo `date  +'%r'`

#3) stats

echo `date  +'%r'`
echo -e "################################\n\tsamtools stats started\n################################\n"
samtools stats $sortedbam > $sample.hisat2.sorted.stats
echo -e "################################\n\tsamtools stats done\n################################\n"
echo `date  +'%r'`


#4)bigWig
effGenomeSize=2913022398 #from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

echo `date  +'%r'`
echo -e "################################\n\tBigWig generation\n################################\n"
module load conda/current
conda activate deeptoolsenv

bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --filterRNAstrand forward --bam $sortedbam -o $sample.sorted.bam.CPM.fwd.bw

bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --filterRNAstrand reverse --bam $sortedbam -o $sample.sorted.bam.CPM.rev.bw
echo -e "################################\n\tsamtools stats done\n################################\n"
echo `date  +'%r'`

