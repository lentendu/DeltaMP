#!/bin/bash

#$ -S /bin/bash
#$ -o log/Illumina_quality.out
#$ -e log/Illumina_quality.err
#$ -l h_rt=24:00:00
#$ -l h_vmem=6G
#$ -binding linear:1

# Define Variables
while read var val; do unset $var ; if [[ $val == "(["* ]]; then declare -A $var="`echo $val | sed 's/].\"/]=\"/g'`" ; else declare $var="$val" ; fi ; done < config/env.txt
. $BIN/check_previous_step

# load modules
module load DeltaMP/${VERSION[DELTAMP]}

# Raw, primer detected, pair-end and trim read counts per library
while read lib fwd rvs
do
	samp=`awk -v L=$lib '{for(i=2;i<=NF;i++){if($i==L){print $1}}}' config/lib4.list`
	paste <(echo -e "$samp\t$lib") <(awk '$0=="+"{a+=1}END{print a}' libraries/$lib$fwd) <(tac libraries/fastq/log.pandaseq.$lib.txt | grep -m 1 -P "STAT\tREADS" | cut -f 4) <(tac libraries/fastq/log.pandaseq.$lib.txt | grep -m 1 -P "STAT\tOK\t" | cut -f 4) quality_check/$lib.stat
done < config/lib2.list | cat <(paste <(echo -e "Sample\tLibrary\tRaw\tPrimer detected\tPair-end") <(seq $MINQUAL 30 | sed 's/\(.*\)/Trimmed at \1 Phred score/' | tr "\n" "\t" | sed 's/$/\n/')) -> quality_check/$SUBPROJECT.lib_counts.tsv
cd quality_check

# per sample
sed '1d' $SUBPROJECT.lib_counts.tsv | sort -k 1,1 | awk 'BEGIN{FS="\t"}{if(NR==1){S=$1;for(i=3;i<=NF;i++){a[i]=$i};printf "%s",S} else {if($1==S){for(i=3;i<=NF;i++){a[i]+=$i}} else {S=$1;for(i=3;i<=NF;i++){printf "\t%s",a[i];a[i]=$i};printf "\n%s",S}}}END{for(i=3;i<=NF;i++){printf "\t%s",a[i]};printf "\n"}' | mimato | cat <(head -1 $SUBPROJECT.lib_counts.tsv | sed 's/\tLibrary\t/\t/') - > $SUBPROJECT.sample_counts.tsv

# Count check
sed '1d' $SUBPROJECT.sample_counts.tsv | tac | sed '1,4d' | awk -v M=$MIN_DEPTH '$3<M{print $1,$3}' > min_depth.raw_with_primer.test
sed '1d' $SUBPROJECT.sample_counts.tsv | tac | sed '1,4d' | awk -v M=$MIN_DEPTH '$4<M{print $1,$4}' > min_depth.pairend.test
sed '1d' $SUBPROJECT.sample_counts.tsv | tac | sed '1,4d' | awk -v M=$MIN_DEPTH '$5<M{print $1,$5}' > min_depth.trim.test
for i in raw_with_primer pairend trim
do
	if [ -s min_depth.$i.test ]
	then
		echo "The following sample(s) do not have the requested minimum amount of $i reads fixed at $MIN_DEPTH :"
		echo "Sample Reads"
		cat min_depth.$i.test
		echo ""
		echo "Aborting"
		awk 'BEGIN{FS="\t";OFS="\t"};{print $1,$2,$3,$4,$5}' $SUBPROJECT.sample_counts.tsv > $SUBPROJECT.summary.stat.tsv
		dd if=$SUBPROJECT.summary.stat.tsv of=../archives/$SUBPROJECT.outputs/$SUBPROJECT.summary.stat.tsv bs=1M
		. $BIN/list_step_files.sh
		exit 100
	else
		echo "All samples have at least $MIN_DEPTH $i sequences"
	fi
done

# Optimize quality
cut -f 5 $SUBPROJECT.sample_counts.tsv | tac | sed '1,4d' | tac | transpose_tab > ref.tmp
tac $SUBPROJECT.sample_counts.tsv | sed '1,4d' | tac | transpose_tab | sed '1,4d;/^$/d' | tac | cat ref.tmp - | awk -v M=$MIN_DEPTH 'BEGIN{FS="\t"} {if(NR==1){for(i=2;i<=NF;i++){ref[i]=$i};next} ; for(i=2;i<=NF;i++){if($i<M || $i<ref[i]*0.75){next}};print $1;exit}' | sed 's/Trimmed at //;s/ Phred score//' > optimized.quality.txt
rm ref.tmp
QUAL=`cat optimized.quality.txt`
cut -f 1-4,$(($QUAL - $MINQUAL + 5)) $SUBPROJECT.sample_counts.tsv > $SUBPROJECT.summary.stat.tsv
dd if=$SUBPROJECT.summary.stat.tsv of=../archives/$SUBPROJECT.outputs/$SUBPROJECT.read_counts.tsv bs=1M

# Merge raw stat statistics and weblogos
cd ..
gs -q -sDEVICE=pdfwrite -o archives/$SUBPROJECT.outputs/$SUBPROJECT.raw_and_pair-end_reads_statistics.pdf libraries/raw_stat/$SUBPROJECT.raw_reads_with_primer_quality.pdf libraries/raw_stat/$SUBPROJECT.pair-end_reads_quality.pdf libraries/raw_stat/$SUBPROJECT.legend.pdf
gs -q -sDEVICE=pdfwrite -dEPSCrop -o archives/$SUBPROJECT.outputs/$SUBPROJECT.forward.weblogo.pdf libraries/raw_stat/weblogo.*.forward.eps
gs -q -sDEVICE=pdfwrite -dEPSCrop -o archives/$SUBPROJECT.outputs/$SUBPROJECT.reverse.weblogo.pdf libraries/raw_stat/weblogo.*.reverse.eps


# list files and directories
. $BIN/list_step_files.sh

echo END
