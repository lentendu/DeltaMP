#!/bin/bash

# 
# DeltaMP, a flexible, reproducible and resource efficient metabarcoding amplicon pipeline for HPC
# Copyright (C) 2018 Guillaume Lentendu, Christina Weißbecker, Anna Heintz-Buschart, Tesfaye Wubet
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 

## Documentation
declare -a TITLE=("Input data" "Pipeline execution" "Benchmarking" "Demultiplexing" "Quality check" "High quality reads processing" "Outputs" "References" "Raw reads extraction")

echo "Documentation of the project $SUBPROJECT produced by DeltaMP v.${VERSION[DELTAMP]}" 
echo "" 

printf '%'$((${#TITLE[0]}+8))'s\n' |tr " " "#"
echo "### ${TITLE[0]} ###"
printf '%'$((${#TITLE[0]}+8))'s\n' |tr " " "#"

if [ $DEMULTI == "no" ]
then
	DEM="not demultiplexed"
	if [ $TECH == "454" ]
	then
		RAWAVG=`grep "^It took [0-9]* secs to extract [0-9]*.$" log/454_demulti.*.out | sed 's/\.$//' | awk '{sum+=$NF}END{printf "%.0f\n", sum/NR}'`
	else
		RAWAVG
	fi
else
	DEM="demultiplexed"
	if [ $TECH == "454" ]
	then
		RAWAVG=`grep "^It took [0-9]* secs to extract [0-9]*.$" log/454_sff.*.out | sed 's/\.$//' | awk '{sum+=$NF}END{printf "%.0f\n", sum/NR}'`
	else
		RAWAVG=`awk '$1=="Average"{printf "%.0f\n", $2}' quality_check/$SUBPROJECT.summary.stat.tsv`
	fi
fi
NBLIB=$LIB1_SIZE
echo "The analysed $RAW_EXT raw reads originated from $NBLIB $DEM library(ies) representing $SAMP_SIZE samples and containing an average of $RAWAVG reads per library."
echo "" 

printf '%'$((${#TITLE[1]}+8))'s\n' |tr " " "#"
echo "### ${TITLE[1]} ###"
printf '%'$((${#TITLE[1]}+8))'s\n' |tr " " "#"
echo "The input sequences were analysed with DeltaMP version ${VERSION[DELTAMP]}, a metabarcode analysis pipeline for grid engines mainly based on vsearch (version ${VERSION[VSEARCH]}, ${CITATION[VSEARCH]}), MOTHUR (version ${VERSION[MOTHUR]}, ${CITATION[MOTHUR]}) and OBITools (version ${VERSION[OBI]}, ${CITATION[OBI]}) software suites."
CIT=(VSEARCH MOTHUR OBI)
echo ""
echo "The pipeline was executed from the directory $INIT_DIR with the following command:"
echo "deltamp $ARGUMENTS"
echo ""
END_TIME=`date +%s`
TIME_ELAPSED=$(($END_TIME-$START_TIME)) #it should be in seconds anyway
echo "The pipeline took $TIME_ELAPSED seconds to complete (real time, from submission to completion) and used a total of $CPU_HOURS cpu-hours."
echo ""

if [ $REF_SUBPROJECT != "no" ]
then
	printf '%'$((${#TITLE[2]}+8))'s\n' |tr " " "#"
	echo "### ${TITLE[2]} ###"
	printf '%'$((${#TITLE[2]}+8))'s\n' |tr " " "#"
	echo "The initial steps of the pipeline were identical to $REF_SUBPROJECT (until ${LAST_REF/script_/} step)."
	echo "deltamp execution began with ${FIRST_STEP/script_/} step."
	echo ""
	echo "The following steps and parameters differed with the previous SUBPROJECTs:"
	sed 's/\('"$SUBPROJECT"'\)/\1 (current)/' config/tree.summary | column -t -s $'\t'
	echo ""
fi

if [ $DEMULTI == "no" ]
then
	printf '%'$((${#TITLE[3]}+8))'s\n' |tr " " "#"
	echo "### ${TITLE[3]} ###"
	printf '%'$((${#TITLE[3]}+8))'s\n' |tr " " "#"
	if [ $TECH == "454" ]
	then
		echo "The demultiplexing was achieved using QIIME (version ${VERSION[QIIME]}, ${CITATION[QIIME]}) script 'make_per_library_sff.py'."
		CIT+=(QIIME)
		echo ""
		if [ $BDIFFS == "a" ]
		then
			echo "The libraries were demultiplexed allowing $PDIFFS mismatches on the primer sequence and the following number of mismatches on the barcode sequence:" 
			cat <(echo "Library bdiffs") <(cut -d " " -f 1,2 processing/bdiffs.$SUBPROJECT) | column -t 
		else
			echo "The libraries were demultiplexed allowing $PDIFFS mismatches on the primer and $BDIFFS mismatches on the barcode." 
		fi
	else
		echo "The demultiplexing was achieved using cutadapt (version ${VERSION[CUT]}, ${CITATION[CUT]}) allowing until $BDIFFS mismatches on the barcode sequences."
		CIT+=(CUT)
		echo ""
	fi

	echo ""
	echo "More informations on demultiplexing efficiency are available in the file demultiplexing_check.csv ."
	echo ""
	echo "The raw demultiplexed libraries are available in the archive $SUBPROJECT.demultiplex_${RAW_EXT}.tar.gz located at $OUT ."
	echo ""
fi

if [ $CLIPPING != "no" ]
then
	printf '%'$((${#TITLE[8]}+8))'s\n' |tr " " "#"
	echo "### ${TITLE[8]} ###"
	printf '%'$((${#TITLE[8]}+8))'s\n' |tr " " "#"
	FWDSIM=`awk -v F=${#FWD} -v D=$PDIFFS 'BEGIN{printf "%.2g\n", 1-D/F}'`
	RVSSIM=`awk -v R=${#RVS} -v D=$PDIFFS 'BEGIN{printf "%.2g\n", 1-D/R}'`
	echo "Read pairs were extracted from raw libraries if at least one of the two reads hold the expected primer (forward primer for forward library, reverse primer for reverse library) at its 5' end, with a similarity threshold of $FWDSIM and $RVSSIM for the forward and reverse primer, respectively."
	CUTAVG=`grep -m 1 "Total read pairs processed:" libraries/fastq/log.cutadapt.* | awk '{gsub(",","",$NF);sum+=$NF}END{printf "%.0f\n", sum/NR}'`
	echo "An average of $CUTAVG reads was extracted per pair of libraries."
	echo ""
fi

printf '%'$((${#TITLE[4]}+8))'s\n' |tr " " "#"
echo "### ${TITLE[4]} ###"
printf '%'$((${#TITLE[4]}+8))'s\n' |tr " " "#"
if [ $TECH == "454" ]
then
	if [ $DEMULTI == "yes" ] && [ $BDIFFS == "a" ] ; then BDIFFS=0 ; fi
	QUAL=`cut -d "." -f 2 quality_check/trimming.parameters.txt`
	LENGTH=`cut -d "." -f 1 quality_check/trimming.parameters.txt`
	echo "The length and average quality trimming parameters were optimised to $LENGTH nt and $QUAL Phred score, respectively, in order to keep at least $MIN_DEPTH trimmed reads per sample, considering all other provided trimming parameters."
	if [ $DENOISE == "yes" ]
	then
		echo "This optimisation does not take into account the loss of reads due to flow trimming and denoising"
	fi
elif [ $TECH == "Illumina" ]
then
	PAIRAVG=`awk '$1=="Average"{printf "%.0f\n", $4}' quality_check/$SUBPROJECT.summary.stat.tsv`
	echo "An average of $PAIRAVG reads were successfully pair-end assembled per sample using the $PALG algorithm with a threshold of $PTRESH and a minimum overlap of $MIN_OV nucleotides as implemented in pandaseq (version ${VERSION[PANDA]}, ${CITATION[PANDA]})."
	CIT+=(PANDA)
	QUAL=`cat quality_check/optimized.quality.txt`
	echo "The average quality trimming parameter was optimised to $QUAL Phred score in order to keep at least $MIN_DEPTH trimmed reads per sample, considering all other provided trimming parameters."
	LENGTH=$MINLEN
fi
TRIMAVG=`awk 'BEGIN{FS="\t"} {if(NR==1){for(i=1;i<=NF;i++){if($i~"^Trimmed")C=i}};if($1=="Average"){printf "%.0f\n", $C}}' quality_check/$SUBPROJECT.summary.stat.tsv`
if [ $TECH == "454" ] ; then if [ $DENOISE == "yes" ] ; then TRIMDEN="denoised and " ; else TRIMDEN="" ; fi ; else TRIMDEN="" ; fi
echo "Each sample contains in average $TRIMAVG ${TRIMDEN}trimmed reads."
echo ""

printf '%'$((${#TITLE[5]}+8))'s\n' |tr " " "#"
echo "### ${TITLE[5]} ###"
printf '%'$((${#TITLE[5]}+8))'s\n' |tr " " "#"
echo "The reads were then ${TRIMDEN}trimmed with the following parameters:"
printf 'minimum length\t%s\nminimum average Phred score on the trimmed length\t%s\nmaximum number of ambiguities in the sequence\t%s\nmaximum length of homopolymers\t%s\n' "$LENGTH" "$QUAL" "$MAXAMBIG" "$MAXHOMOP" | column -t -s $'\t' 
if [ $TECH == "454" ]
then
	printf 'maximum number of mismatch(es) on the primer sequence\t%s\n' "$PDIFFS" | column -t -s $'\t' 
	if [ $BDIFFS == "a" ]
	then
		echo "maximum number of mismatches allowed on the barcode sequence (in each sequence library separatedly):"
		cat <(echo "Library bdiffs") <(cut -d " " -f 1,2 processing/bdiffs.$SUBPROJECT) | column -t
	else
		printf 'maximum number of mismatches allowed on the barcode sequence\t%s\n' "$BDIFFS"
	fi
	if [ $DENOISE == "yes" ]
	then
		printf "minimum flow length\t%s\nmaximum flow length\t%s\n\n" "$MINFLOW" "$MAXFLOW"
		echo "Reads were trimmed and flows were denoised using FlowClus (version ${VERSION[DEN]}, ${CITATION[DEN]})."
		echo "As FlowClus do not allow for homopolymer extension or reduction during the mismatch detection in the primer sequence, cutadapt (version ${VERSION[CUT]}, ${CITATION[CUT]}) was first use to detect all primer variants in the error range allowed."
		CIT+=(DEN CUT)
	fi
	if [ $ITSX == "no" ]
	then
		echo ""
		echo "The first $LENGTH nucleotides of the reads were kept for further analysis"
	fi
fi
echo ""

if [ $SUBSAMPLE == "yes" ]
then
	SUBSIZE=$(grep "Minimum" quality_check/$SUBPROJECT.summary.stat.tsv | cut -f 3 | cut -d "." -f 1)
	echo "The read count was randomly normalized to $SUBSIZE reads in each samples."
	echo ""
fi
if [ $PRECLUST == "mothur" ]
then
	PRECL=`for i in log/trim.*.out; do grep -B 1 "^pre.cluster removed" $i | grep -o '[0-9]\+' | paste - - | awk '{print $2/$1*100}' ; done | awk '{m+=$1}END{print m/NR}'`
	echo "Dereplicated sequences of each sample were aligned against the reference alignment of the database $DB (version ${VERSION[DB]}, ${CITATION[DB]}) and pre-clustered following the mothur $TARG SOP (${CITATION[SOP]})."
	echo "Approximately 5% of the reads with bad alignment were discarded after the alignment and pre-clustering reduced the number of sequences to analyse by an average of $PRECL %."
	echo ""
	CIT+=(SOP)
elif [ $PRECLUST == "cdhit454" ]
then
	PRECL=`for i in log/trim.*.out; do grep "^ *[0-9].*clusters$" $i ; done | awk '{mean+=$3/$1*100}END{print mean/NR}'`
	echo "The reads were pre-clustered in order to merge reads likely arising from sequencing errors (${CITATION[PRECL]}) and thus reducing the computational load."
	echo "cd-hit-454 was used for that purpose, allowing a maximum of 1 % dissimilarity and with only one base allowd per indel (version ${VERSION[PRECDHIT]}, ${CITATION[PRECDHIT]})."
	echo "This step reduced the number of sequences to analyse by an average of $PRECL %."
	echo ""
	CIT+=(PRECL PRECDHIT)
fi

if [ $TARG == "ITS" ] && [ $ITSX != "no" ]
then
	ITSXMEAN=$(while read num samp; do NAMES=$(tac log/trim.$num.out | sed -n "1,/unique\.seqs/{s/^.*name=\([^)]*\))$/processing\/$samp\/\1/p}") ; echo $(sed -n '$=' $NAMES) $(sed -n '$=' ${NAMES/itsx\.$ITSX\./}) ; done < <(awk '{print NR,$1}' config/lib4.list) | awk '{sum+=$1/$2*100}END{print sum/NR}')
	echo "The $ITSX fragment was detected and extracted from $ITSXMEAN % of chimera free reads using ITSx (version ${VERSION[ITSX]}, ${CITATION[ITSX]}). "
	CIT+=(ITSX)
	echo ""
elif [ $CHIMERA1 == "yes" ]
then
	CHIMMEAN=$(grep "Removed [0-9]* sequences from your name file" log/trim.[0-9]*.out | cut -d " " -f 2 | awk '{c+=$1}END{printf "%.0d\n", c/NR}')
	echo "An average of $CHIMMEAN chimeric reads were detected and removed from each sample using the UCHIME algorithm as implemented in MOTHUR (${CITATION[UCHIME]})."
	CIT+=(UCHIME)
	echo ""
fi

NBREADS=$(awk 'NR>1{sum+=$2}END{print sum}' processing/$SUBPROJECT.unique.count_table)
NBUNIQ=$(sed -n '$=' processing/$SUBPROJECT.names)
echo "Reads from each sample were pooled together, representing a total of $NBREADS high-quality reads, and were further used in all downstream analyses."
echo ""
echo "The reads were dereplicated into $NBUNIQ unique sequences and then sorted in their decreasing abundance order before clustering."
echo ""

while read var val; do unset $var ; if [ $REF_SUBPROJECT == "no" ] ; then declare $var="$val" ; else declare $var="${val//$REF_SUBPROJECT/$SUBPROJECT}" ; fi ; done < config/OTU_env.txt
if [ -s processing/$ACCNOS.accnos ]
then
	NBCHIM=`sed -n '$=' processing/$ACCNOS.accnos`
	NBOTUS=$(($(cut -f 2 processing/$LIST.list | sed '1d')+$NBCHIM))
	NBRM=$(($NBREADS-$(cut -f 2 processing/$NAMES.names | tr "," "\n" | wc -l)))
	echo "The dereplicated reads were cluster into $NBOTUS OTUs using the $CLUST algorithm (version ${VERSION[CLUST]}, ${CITATION[CLUST]})."
	if [ $CLUST == "cd-hit-est" ] && [ $PREV_PATH != "no" ]
	then
		PREV_NB=`sed -n '$='  processing/${PREV_PATH##*/}.*.cdhit.names`
		echo "First, the $PREV_NB representative sequence of the previous subproject ${PREV_PATH##*/}, located in ${PREV_PATH%/*}, were used as OTU seed using cd-hit-est-2d."
		echo "Then, reads which failed to cluster with previous subproject representative sequence were de-novo clustered into OTU with cd-hit-est and add to the previous OTUs, starting the at Otu$(( $PREV_NB + 1 ))."
	fi
	echo "$NBCHIM chimeras were detected among the OTU representative sequences using the UCHIME algorithm as implemented in MOTHUR (${CITATION[UCHIME]})."
	echo "Those $NBCHIM OTUs, representing $NBRM sequences, were removed from the final dataset."
	echo ""
else
	NBOTUS=`cut -f 2 processing/$LIST.list | sed '1d'`
	echo "The dereplicated reads were cluster into $NBOTUS OTUs using the $CLUST algorithm (version ${VERSION[CLUST]}, ${CITATION[CLUST]})."
	if [ $CLUST == "cd-hit-est" ] && [ $PREV_PATH != "no" ]
	then
		PREV_NB=`sed -n '$='  processing/${PREV_PATH##*/}.*.cdhit.names`
		echo "First, the $PREV_NB representative sequence of the previous subproject ${PREV_PATH##*/}, located in ${PREV_PATH%/*}, were used as OTU seed using $CLUST-2d."
		echo "Then, reads which failed to cluster with previous subproject representative sequence were de-novo clustered into OTU with $CLUST and add to the previous OTUs, starting at Otu$(( $PREV_NB + 1 ))."
	fi
	echo "No chimeras could be detected among the OTU representative sequences using the UCHIME algorithm as implemented in MOTHUR (${CITATION[UCHIME]})."
	echo ""
fi
CIT+=(CLUST)
if [ $ASSIGN_ALL == "yes" ]
then
	echo "All dereplicated reads were taxonomically assigned based on the reference sequences of the $DB database (version ${VERSION[DB]}, ${CITATION[DB]}) using the naive bayesian classifier (${CITATION[TAXO]}), as implemented in MOTHUR, at a consensus threshold of 60%."
	echo "The OTUs were assigned to the longest taxonomy path shared by at least 60% of their reads."
else
	echo "The OTU representative sequences (the most abundant sequence in each OTU) were taxonomically assigned based on the reference sequences from the $DB database (version ${VERSION[DB]}, ${CITATION[DB]}) using the naive bayesian classifier (${CITATION[TAXO]}), as implemented in MOTHUR, at a consensus threshold of 60%."
fi
CIT+=(DB TAXO)
if [ $TARG == "ITS" ]
then
	echo "The fungal sequences not assigned til genus were assigned using the dynamic version of $DB database which include singletons. The first assignment against the dynamic $DB database without singletons was replaced by this second assignement only if the taxonomy could be improved (i.e. assigned at deeper taxonomic level)."
	echo "Putative functions were annotated using the FUNGuild fungal database (version ${VERSION[FUN]}, ${CITATION[FUN]})."
fi
CIT+=(FUN)
LABELS=(all abundant all_$TARG_ORG abundant_$TARG_ORG)
NBOTUS=(`for i in "${LABELS[@]}"; do sed '1d' processing/$SUBPROJECT.${i}_OTUs.tsv | sed -n '$=' ; done`)
NBREADS=(`grep "^Total" processing/$SUBPROJECT.read_counts.tsv | cut -f 4- | sed 's/\.[0-9]*//g'`)
echo ""
echo "The final OTU table contains ${NBOTUS[0]} OTUs, representing ${NBREADS[0]} sequences."
echo "From it, ${NBOTUS[1]} OTUs, representing ${NBREADS[1]} sequences, are abundant (present in at least $MIN_SAMP samples and with at least $MIN_DOM sequences)."
echo "${NBOTUS[2]} OTUs, representing ${NBREADS[2]} sequences, were assigned to the target lineage $TARG_ORG."
echo "${NBOTUS[3]} OTUs, representing ${NBREADS[3]} sequences, are abundant and were assigned to the target lineage $TARG_ORG."
echo ""

printf '%'$((${#TITLE[6]}+8))'s\n' |tr " " "#"
echo "### ${TITLE[6]} ###"
printf '%'$((${#TITLE[6]}+8))'s\n' |tr " " "#"
echo "The OTU tables and representative sequences were outputted in the following files:"
for i in "${LABELS[@]}"; do echo -e "${i/_/ } OTUs:\t$SUBPROJECT.${i}_OTUs.tsv\t$SUBPROJECT.${i}_repseq.fasta" ; done | column -t -s $'\t'
echo ""
echo "The full OTU table was additionnaly stored in the json BIOM formated file $SUBPROJECT.json.biom using the BIOM format (version ${VERSION[BIOM]}, ${CITATION[BIOM]})."
CIT+=(BIOM)
if [ $TARG == "ITS" ]
then
	echo "The fungal functional annotations (functional_assignment_level, trophicMode, guild, growthForm and trait) were inserted in the BIOM file as observation metadata."
fi
if [[ $LIB_DIR != "/"* ]]
then
	echo "The BIOM file also includes all sample's metadata associated with their respective accession in the BioProject $LIB_DIR."
fi
echo ""
echo "Per sample read counts for raw, trimmed and the four output variants were written into the file $SUBPROJECT.read_counts.tsv ."
echo ""
echo "The original configuration file was saved under the filename metadata.$SUBPROJECT.tsv ."
echo ""
echo "Software's log files and verbosities from each steps were pooled in their execution order in the file $SUBPROJECT.log ."
echo ""
echo "All additionnal processing files produced from the pooling of each sample high quality reads to the end of the pipeline were archive in $SUBPROJECT.processing.files.tar.gz ."
echo ""

printf '%'$((${#TITLE[7]}+8))'s\n' |tr " " "#"
echo "### ${TITLE[7]} ###"
printf '%'$((${#TITLE[7]}+8))'s\n' |tr " " "#"
for i in "${CIT[@]}"; do echo ${FULLCITATION[$i]} ; done | sort -k 1,1 -u

