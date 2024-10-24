#!/bin/bash
die() {
	printf '%s\n' "$1" >&2
    exit 1
    }
	

show_help() {
cat << EOF
		Usage: ${0##*/} [options] [working_dir]...
		Run ATAC-seq basic data processing in CSUC hpc;
		 from QC to alignment to generation of bigWig files for visualization 
		 
		Available options are:
			
			-h			display this help and exit
			-rs runstep		run the analysis starting from this particular step. All steps are run by default (runstep=all).
			-only			if this flag is set (no additional arguments) then only the step indicated by -rs will be run (default FALSE)
			-ppn integer		number of proccesors per node (default 10)
		   
EOF
}
	
runstep=all
PPN=10
only=FALSE

#parse optional arguments	   

while :; do
	case $1 in
		-h|-\?|--help)
			show_help    # Display a usage synopsis.
			exit
			;;
		-rs|--runstep)       # Takes an option argument; ensure it has been specified.
			re='^-'
			if [ "$2" ] &&  ! [[ $2 =~ $re ]]; then
				case $2 in 
				  "all"|"qc"|"trim"|"align"|"sort_index"|"markdups"|"filtstats"|"bigWig")
					  runstep=$2
					  echo "running step $runstep"
				  ;;
				  *)
				    die 'ERROR: non-valid runnning step
					"-rs" requires one of these non-empty option arguments : "all", "qc", "trim", "align", "sort_index", "markdups", "filtstats" or "bigWig" .'
				  ;;
				esac
	            shift 2
	        else
	            die 'ERROR: "-rs" requires one of these non-empty option arguments : "all", "qc", "trim", "align", "sort_index" .'
			fi
			;;
		-ppn)
			if [ "$2" ]; then
				PPN=$2
				re='^[0-9]+$'
				if ! [[ $PPN =~ $re ]] ; then
				   die 'ERROR:  "-ppn" must be followed by an integer' 
				fi
	            shift 2
	        else
	            die 'ERROR: "-ppn" requires an integer as option argument'
			fi
			;;
		-only)
			only=TRUE
			shift
			echo "Only $runstep will be run"
			;;
		-?*)
			printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
			show_help
			shift
	    	;;
	    *)  # Default case: No more options, so break out of the loop.
	    	break
	esac
done	    

echo "PPN: $PPN"
echo "only: $only"
echo "wdir: $1"

if [ $runstep == all ] && [ $only ]
then
	only=FALSE
	echo 'Note that when -rs all then -only must by FALSE, setting -only back to FALSE'
fi
#Now the actual analysis starts 
#NOTE: this is a preliminary version, in the future I would like to make a function out of each step for the sake of organisation

sdir=$1
sample=$(basename $sdir) 


#indexes and ref files.These are human files, they can be changed if working with mouse or another species
bowtie2Index="/scratch/llorenzi/references/human/ncbi/bowtie2_index/GRCh38_noalt_as/GRCh38_noalt_as"
primary_chrs_bed="/scratch/llorenzi/references/human/ncbi/primary_chrs_nochrM.bed"

#set tool paths
bowtie2="/home/llorenzi/bowtie2-2.4.2-linux-x86_64/bowtie2"
picard="/home/llorenzi/picard.jar"
trim_galore="/home/llorenzi/TrimGalore-0.6.6/trim_galore"

cd $sdir

#1) QC
step=qc
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1=$(ls *_1.fq.gz)
	read2=$(ls *_2.fq.gz)
	
	if ! [ -f "$read1" ]; then
		die "input file $read1 not found"
	fi
	read2=$(ls *_2.fq.gz)
	if ! [ -f "$read2" ]; then
		die "input file $read2 not found"
	fi
	echo "read1: $read1"
	echo "read2: $read2"
	
	echo `date`
	echo -e "################################\n\tFASTQC started\n################################\n"
	module load fastqc/0.11.8
	mkdir fastqc
	fastqc $read1 -o fastqc
	fastqc $read2 -o fastqc
	echo -e "################################\n\tFASTQC done\n################################\n"
	echo `date  +'%r'`
	if [ $only ]
	then 
		exit 0
	fi
fi


#2)trim_galore
step=trim
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1=$(ls *_1.fq.gz)
	if ! [ -f "$read1" ]; then
		die "input file $read1 not found"
	fi
	read2=$(ls *_2.fq.gz)
	if ! [ -f "$read2" ]; then
		die "input file $read2 not found"
	fi
	echo "read1: $read1"
	echo "read2: $read2"
	
	echo `date  +'%r'`
	echo -e "################################\n\tTrimGalore started\n################################\n"

	module load conda/current
	conda activate cutadaptenv

	$trim_galore --paired --length 35 --fastqc $read1 $read2
	echo -e "################################\n\tTrimGalore done\n################################\n"
	echo `date  +'%r'`
	if [ $only ]
	then 
		exit 0
	fi
fi
	

#3)alignment
step=align
if [ $runstep == all ] || [ $runstep == $step ]
then
	read1_val=$(ls *val_1.fq.gz)
	if ! [ -f "$read1_val" ]; then
		die "input file $read1_val not found"
	fi
	read2_val=$(ls *val_2.fq.gz)
	if ! [ -f "$read2_val" ]; then
		die "input file $read2_val not found"
	fi
	echo "read1: $read1_val"
	echo "read2: $read2_val"
	
	echo `date  +'%r'`
	echo -e "################################\n\tBowtie2 started\n################################\n"
	$bowtie2 --very-sensitive -X 2000 -p $PPN -x $bowtie2Index -1 $read1_val -2 $read2_val -S $sample.bt2.sam
	echo -e "################################\n\tBowtie2 done\n################################\n"
	echo `date  +'%r'`

	
	if [ $only ]
	then 
		exit 0
	fi
fi
	
#4) convert to bam while sorting and index
step=sort_index
if [ $runstep == all ] || [ $runstep == $step ]
then
	input=$sample.bt2.sam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	echo `date  +'%r'`
	echo -e "################################\n\tsamtools sort started\n################################\n"
	module load samtools/1.9
	
	samtools sort -@ $PPN -O bam -o $sample.sorted.bam $input
	samtools index -@ $PPN $sample.sorted.bam
	
	echo -e "################################\n\tsamtools sort done\n################################\n"
	echo `date  +'%r'`
	
	if [ $only ]
	then 
		exit 0
	fi
fi



#5) mark duplicates
step=markdups
if [ $runstep == all ] || [ $runstep == $step ]
then 
	PATH=/home/llorenzi/jdk-16.0.1/bin:$PATH
	export PATH
	
	input=$sample.sorted.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	echo `date  +'%r'`
	echo -e "################################\n\tpicard Mark duplicates started\n################################\n"
	java -jar $picard MarkDuplicates INPUT=$input OUTPUT=$input.markedDups.bam METRICS_FILE=$sample.dupmatrix TMP_DIR=/tmp/ 	CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=FALSE TAGGING_POLICY=All
	echo -e "################################\n\tpicard Mark duplicates done\n################################\n"
	echo `date  +'%r'`
	
	if [ $only ]
	then 
		exit 0
	fi
fi


#6) MAPQ distribution, raw and proper pairs, filtering and stats
step=filtstats
if [ $runstep == all ] || [ $runstep == $step ]
then 
	input=$sample.sorted.markedDups.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi
	
	module load samtools/1.9
	
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
	
	if [ $only ]
	then 
		exit 0
	fi
fi

#7) bigwig
step=bigWig
if [ $runstep == all ] || [ $runstep == $step ]
then 
	input=$sample.sorted.markedDups.primary_chr.proper_pairs.minq2.bam
	if ! [ -f "$input" ]; then
		die "input file $input not found"
	fi

	module load conda/current
	conda activate deeptoolsenv
	echo `date  +'%r'`
	echo -e "################################\n\tbamcoverage started\n################################\n"
	bamCoverage --binSize 10 --normalizeUsing CPM --bam $input -o $input.CPM.bw
	echo -e "################################\n\tbamcoverage done\n################################\n"
	echo `date  +'%r'`
	
	if [ $only ]
	then 
		exit 0
	fi
fi

		
 
