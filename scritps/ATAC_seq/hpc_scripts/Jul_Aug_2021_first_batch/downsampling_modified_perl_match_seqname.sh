#!/bin/sh
#########################################################################################
# Author: Haibo Liu              
# script_name: downsampling.sh
# Usage: sh downsampling.sh    /path/to/BAM_file  suffix  a_sequence_of_fractions_between_0_and_1 
# Example:  sh downsampling.sh    SRR891270.chr.filtered.sorted.rmdup.bam  .chr.filtered.sorted.rmdup.bam  0.1 0.2 0.3 0.4 0.5  0.6  0.7  0.8  0.9
#                                 
#########################################################################################

USAGE="Usage: sh $0 /path/to/BAM_file a_sequence_of_fractions_between_0_and_1"


## checking number of command-line arguments
if [ $# -lt 3 ]; then
    echo "At least three arguments are needed!"
	echo "$USAGE"
	exit 1
fi


## getting the basename of a bam file
base=`basename $1 $2`


## all arguments from $3 to ${N}, holding subsampling percentage
percentage=("${@:3}")

## randomizing alignments by sorting by read names
samtools sort -l 9 -m 8G  -o  ${base}.name.sorted.sam  -O SAM  -n -@ 8  $1

## getting SAM header
grep  -P '^@(HD|SQ|PG)'  ${base}.name.sorted.sam > ${base}.name.sorted.sam.header


## converting paired read alignment to a single line representing a fragment alignment
grep -v -P '^@(HD|SQ|PG)'  ${base}.name.sorted.sam | awk 'BEGIN{FS=OFS="\t"; ORS="\t"; i =1} \
        {if (i==1) {print $0; a=$1; i=i+1} else if ($1==a) {print $0; i= i+1} else if ($1 !=a) \
        {print "\n"; print $0; i= i+1; a= $1}}'| perl -p -e 's/^\t//; s/\t$//' >  ${base}.name.sorted.oneliner.sam


## subsampling by using Unix built-in shuf

total=`wc -l ${base}.name.sorted.oneliner.sam | awk  '{print $1}' `

while (($#)); do
   subsample_size+=( `echo "${total} * ${1} " |bc | perl -p -e 's/\.\d+//'` )
   shift
done

echo "${subsample_size[@]}"


## number of subsamples
num_subsamples=${#subsample_size[@]}

for (( i=0; i<${num_subsamples}; i++ ));

do
	
	file=${base}.name.sorted.${percentage[${i}]}.sam 
	cp ${base}.name.sorted.sam.header   $file
	shuf -n  ${subsample_size[$i]} ${base}.name.sorted.oneliner.sam  >>  $file
	
	base1=`basename $file .sam`
    	perl -p -e 's/(V300.+?)\s+(V300.+)/$1\n$2/' $file   | samtools view -b -h  -@ 8 - | \
	samtools sort  -l 9 -m 8G  -o  ${base1}.coord.sorted.bam  -O BAM -@ 8 -
    	samtools index  -@ 1  ${base1}.coord.sorted.bam

done

