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

mkdir $EXEC/libraries/demultiplex

echo "Library Raw_reads_count Provided_barcode_count Provided_barcode_found_in_library Provided_barcode_total_read_count Provided_barcode(s)_not_found Other_barcode(s)_in_library Other_barcode(s)_total_read_count" | tr " " "\t" | tr "_" " "
while read ORI ORIF LIB_NAME FWD_SUF RVS_SUF
do
	FWD_ORI=${ORI}${ORIF}
	BARCODES=(`awk -v L=$FWD_ORI '$3~L{print $1}' config/barcodes.txt`)
	FWD_LIB=${LIB_NAME}${FWD_SUF/.$RAW_EXT.*/}
	RVS_LIB=${LIB_NAME}${RVS_SUF/.$RAW_EXT.*/}
	cd libraries
	SUMALL=`sed -n '$=' ${FWD_LIB}.fastq | awk '{print $1/4}'`

	# list all potential barcodes (with pair-end exact match control)
	paste <(awk -v L=${#BARCODES[0]} '(NR+2) % 4 == 0{print substr($1,1,L)}' ${FWD_LIB}.fastq) <(awk -v L=${#BARCODES[0]} '(NR+2) % 4 == 0{print substr($1,1,L)}' ${RVS_LIB}.fastq) | awk '$1 == $2{print $1}' | sort | uniq -c > demultiplex/${LIB_NAME}.barcodes.tmp
	
	# matched barcodes
	cd demultiplex
	echo ${BARCODES[@]} | tr " " "\n" | sort | join -2 2 - ${LIB_NAME}.barcodes.tmp > ${LIB_NAME}.match.tmp
	
	# unmatched barcodes
	echo ${BARCODES[@]} | tr " " "\n" | sort | join -v 1 -2 2 - ${LIB_NAME}.barcodes.tmp > ${LIB_NAME}.unmatch.tmp
	UNM=$(tr "\n" "," < ${LIB_NAME}.unmatch.tmp | sed 's/,*$/\n/')
	
	# other valid barcodes ( >= 1% of raw reads)
	echo ${BARCODES[@]} | tr " " "\n" | sort | join -v 2 -2 2 - ${LIB_NAME}.barcodes.tmp | awk -v A=$SUMALL '$2>=(A/100)' > ${LIB_NAME}.other.tmp
	OTH=$(tr "\n" "," < ${LIB_NAME}.other.tmp | sed 's/,*$/\n/')
	if [ -s ${LIB_NAME}.other.tmp ] ; then OTHN=$(awk '{s+=$2}END{print s,s/NR}' ${LIB_NAME}.other.tmp) ; else OTHN=null ; fi
	
	# output
	echo "$LIB_NAME $SUMALL ${#BARCODES[@]} $(sed -n '$=' ${LIB_NAME}.match.tmp) $(awk '{s+=$2}END{print s}' ${LIB_NAME}.match.tmp) ${UNM:-none} ${OTH:-none} $OTHN" | tr " " "\t"
	rm ${LIB_NAME}.barcodes.tmp ${LIB_NAME}.match.tmp ${LIB_NAME}.unmatch.tmp ${LIB_NAME}.other.tmp
	cd $EXEC
done < <(paste <(cut -f 1-2 config/lib1.list) config/lib2.list | sed 's/ftp\.sra[^\t]*\///g')

rmdir $EXEC/libraries/demultiplex
