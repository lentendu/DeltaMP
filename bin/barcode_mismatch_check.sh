#!/bin/bash

# author: Guillaume Lentendu (guillaume.lentendu@ufz.de)

	if [ -z "$1" ]; then 
		echo usage: $0 MOTHUR-oligo-file
		echo "Determine the allowed amount of mismatches (only indels checked) to use for demultiplexing a 454 library (from 0 to 2)"
		exit

	fi

OLIGO=$1 ; shift
BAR=(`awk '$1~"barcode"{print $2}' $OLIGO | sed 's/\(.\)/\1#/g;s/#$//' | tr "\n" " " | sed 's/$/\n/'`)
OUT=0

# One mismatch dictionary
echo -n > one.mismatch.$OLIGO
for g in "${BAR[@]}"
do
	TMP=(`echo $g | sed 's/#/ /g'`)
	for i in $(seq 1 ${#TMP[@]})
	do
		for j in A T C G;
		do
			if [[ $j == ${TMP[i-1]} ]]
			then
				if [[ $i -eq 1 ]]
				then
					echo ${TMP[@]} | sed 's/ //g' >> one.mismatch.$OLIGO
				fi
			else
				CHANGE=(${TMP[@]})
				CHANGE[i-1]=$j
				echo ${CHANGE[@]} | sed 's/ //g' >> one.mismatch.$OLIGO
			fi
		done
	done
done
ALL=`sed -n '$=' one.mismatch.$OLIGO`
UNIQ=`sort one.mismatch.$OLIGO | uniq | sed -n '$='`
if [[ $ALL -le $UNIQ ]]
then
	unset OUT && OUT=1
else
	echo $OUT
	rm one.mismatch.$OLIGO
	exit
fi

# Two mismatches dictionary
echo -n > two.mismatch.$OLIGO
for g in "${BAR[@]}"
do
	TMP=(`echo $g | sed 's/#/ /g'`)
	for i in $(seq 1 ${#TMP[@]})
	do
		for j in A T C G;
		do
			if [[ $j == ${TMP[i-1]} ]]
			then
				if [[ $i -eq 1 ]]
				then
					echo ${TMP[@]} | sed 's/ //g' >> two.mismatch.$OLIGO
				fi
			else
				CHANGE=(${TMP[@]})
				CHANGE[i-1]=$j
				echo ${CHANGE[@]} | sed 's/ //g' >> two.mismatch.$OLIGO
				for h in $(seq 1 ${#TMP[@]})
				do
					if [[ $h -gt $i ]]
					then
						for k in A T C G;
						do
							if [[ $k != ${TMP[h-1]} ]]
							then
								CHANGE[h-1]=$k
								echo ${CHANGE[@]} | sed 's/ //g' >> two.mismatch.$OLIGO
							fi
						done
					fi
				done
			fi
		done
	done
done
ALL=`sed -n '$=' two.mismatch.$OLIGO`
UNIQ=`sort two.mismatch.$OLIGO | uniq | sed -n '$='`
if [[ $ALL -le $UNIQ ]]
then
	unset OUT && OUT=2
fi
echo $OUT
rm two.mismatch.$OLIGO
