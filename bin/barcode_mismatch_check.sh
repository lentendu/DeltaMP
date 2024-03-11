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

if [ -z "$1" ]; then 
	echo usage: $0 MOTHUR-oligo-file
	echo "Determine the allowed amount of mismatch(es) to use for demultiplexing (0 to 2)"
	exit
fi

OLIGO=$1 ; shift

# One mismatch dictionary
declare -A ATCG=( [A]={C,G,T} [C]={A,G,T} [G]={A,C,T} [T]={A,C,G} )
while read g
do
	echo $g
	for i in $(seq 0 $((${#g}-1)))
	do
		eval echo ${g:0:$i}${ATCG[${g:$i:1}]}${g:$((i+1)):$((${#g}-$i))}
	done
done < <(awk '$1~"barcode"{print $2}' $OLIGO) | tr " " "\n" > one.mismatch.$OLIGO
ALL=`sed -n '$=' one.mismatch.$OLIGO`
UNIQ=`sort -u one.mismatch.$OLIGO | sed -n '$='`
if [[ $UNIQ -lt $ALL ]]
then
	echo 0
	rm one.mismatch.$OLIGO
	exit
fi

# Two mismatches dictionary
while read g
do
	for i in $(seq 0 $((${#g}-2)))
	do
		for j in $(seq $((i+1)) $((${#g}-1)))
		do
			eval echo ${g:0:$i}${ATCG[${g:$i:1}]}${g:$((i+1)):$(($j-$i-1))}${ATCG[${g:$j:1}]}${g:$((j+1)):$((${#g}-$j))}
		done
	done
done < <(awk '$1~"barcode"{print $2}' $OLIGO) | tr " " "\n" | cat one.mismatch.$OLIGO - > two.mismatch.$OLIGO
ALL=`sed -n '$=' two.mismatch.$OLIGO`
UNIQ=`sort two.mismatch.$OLIGO | uniq | sed -n '$='`
if [[ $UNIQ -lt $ALL ]]
then
	echo 1
else
	echo 2
fi
rm one.mismatch.$OLIGO two.mismatch.$OLIGO
