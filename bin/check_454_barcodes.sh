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

LIB_ALL=($@)

echo "Library Raw_reads_count Raw_reads_with_forward_primer Raw_reads_with_one_of_the_provided_barcode Provided_barcode_count Found_barcodes_with_forward_primer_in_library Average_read_count_per_provided_barcodes Unused_provided_barcode(s) Unused_barcode(s)_with_forward_primer_in_library" | tr " " "," | tr "_" " "
for i in "${LIB_ALL[@]}"
do
	unset BARCODES BAR_REGEX UNMATCHED UNM_REGEX
	# barcodes with 1 mismatch allowed
	BARCODES=(`awk '$1=="barcode"{print $2}' oligos.$i`)
	if [ ${#BARCODES[0]} -le 6 ]
	then
		BAR_REGEX=(${BARCODES[@]})
	else
		for j in $(seq 0 $((${#BARCODES[@]}-1))) ; do BAR_REGEX[$j]=`for k in $(seq 0 $((${#BARCODES[$j]}-1))); do echo "^${BARCODES[$j]:0:k}[ATCG]\{0,2\}${BARCODES[$j]:k+1}" ; done | tr "\n" "|" | sed 's/|$//;s/|/\\\|/g'` ; done
	fi
	# Keep only sequences beginning with ${#BARCODES[0]} (+/- 2) nucleotides + the five first nucleotides of the forward primer (with 1 mismatch allowed)
	FWD_START=`awk '$1=="forward"{print substr($2,1,5)}' oligos.$i`
	FWD_REGEX=`for j in $(seq 0 $((${#FWD_START}-1))); do echo "^[ATCG]\{$((${#BARCODES[0]}-1)),$((${#BARCODES[0]}+1))\}${FWD_START:0:j}[ATCG]\{0,2\}${FWD_START:j+1}" ; done | tr "\n" "|" | sed 's/|$//;s/|/\\\|/g'`
	grep "$FWD_REGEX" $i.fasta > $i.tmp
	# For each barcode, print barcode and number of sequences (1 mismatch allowed)
	for j in $(seq 0 $((${#BARCODES[@]}-1))) ; do grep -c "^${BAR_REGEX[$j]}" $i.tmp | paste <(echo ${BARCODES[$j]}) - ; done > $i.match.tmp
	# Mean count of barcodes
	BAR_MEAN=`awk '$2>0{t+=$2}END{printf "%.d", t/NR}' $i.match.tmp`
	BAR_TOTAL=`cut -c 1-${#BARCODES[0]} $i.tmp | sort | uniq -c | awk -v M=$BAR_MEAN '$1>=M/10{c+=1}END{print c}'`
	# Unmatched barcodes with at least 10% of mean count of matched barcodes
	UNMATCHED=(`grep -v $(echo ${BAR_REGEX[@]} | sed 's/ /\\\|/g') $i.tmp | cut -c 1-${#BARCODES[0]} | sort | uniq -c | awk -v M=$BAR_MEAN '$1>=M/10{print $2}'`)
	for j in $(seq 0 $((${#UNMATCHED[@]}-1))) ; do UNM_REGEX[$j]=`for k in $(seq 0 $((${#UNMATCHED[$j]}-1))); do echo "^${UNMATCHED[$j]:0:k}[ATCG]\{0,2\}${UNMATCHED[$j]:k+1}" ; done | tr "\n" "|" | sed 's/|$//;s/|/\\\|/g'` ; done
	for j in $(seq 0 $((${#UNMATCHED[@]}-1))) ; do grep -c "^${UNM_REGEX[$j]}" $i.tmp | paste <(echo ${UNMATCHED[$j]}) - ; done > $i.unmatch.tmp
	# Output
	echo "$i $(grep -c ">" $i.fasta) $(sed -n '$=' $i.tmp) $(awk '{t+=$2}END{print t}' $i.match.tmp) ${#BARCODES[@]} $BAR_TOTAL $BAR_MEAN $(awk '$2==0{printf "%s#", $1}' $i.match.tmp | sed 's/#$//') $(tr "\n" "#" < $i.unmatch.tmp | sed 's/\t/(/g;s/#/)#/g;s/#$//')" | tr " " ","
	rm $i.tmp $i.match.tmp $i.unmatch.tmp
done | awk 'BEGIN{FS=","} {if ($0~"#") {
		for (i=1;i<=7;i++) {
			printf "%s,", $i
		}
		x=split($8,a,"#")
		y=split($9,b,"#")
		printf "%s,%s\n", a[1],b[1]
		if(x>y) {
			for (i=2;i<=y;i++) {
				printf ",,,,,,,%s,%s\n", a[i],b[i]
			}
			for (i=y+1;i<=x;i++) {
				printf ",,,,,,,%s,\n", a[i]
			}
		} else {
			for (i=2;i<=x;i++) {
				printf ",,,,,,,%s,%s\n", a[i],b[i]
			}
			for (i=x+1;i<=y;i++) {
				printf ",,,,,,,,%s\n", b[i]
			}			
		}
	} else {
		print $0
	}
}'

