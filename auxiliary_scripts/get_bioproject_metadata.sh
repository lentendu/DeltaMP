#!/bin/bash

# author: Guillaume Lentendu (guillaume.lentendu@unine.ch)

if [ -z "$*" ]
then
	echo "usage: ${0##*/} Bioproject_accession"
	exit 0
fi

bioproject=$1

# retrieve bioproject runs metadata
wget -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${bioproject}&result=read_run&fields=all&format=tsv&download=true" | transpose_tab | awk 'NF>1{print}' | transpose_tab > ${bioproject}.libraries.tsv

mkdir metadata.${bioproject} && cd metadata.${bioproject}
# retrieve bioproject samples metadata
awk 'BEGIN{FS=OFS="\t"} {if(NR==1){for(i=1;i<=NF;i++){if($i=="sample_alias"){s=i};if($i=="sample_accession"){a=i}}} else {print $s,$a}}' ../${bioproject}.libraries.tsv | sed 's/ /_/g' | sort -u > ${bioproject}.acc.txt
while read sample accession
do
	wget -nv -O $sample.xml http://www.ebi.ac.uk/ena/browser/api/xml/$accession?download=true
done < ${bioproject}.acc.txt

# format sample's xml metadata to tab delimited
grep "<TAG>" *.xml | sed 's/^.*<TAG>\([^<]*\)<\/TAG>.*$/\1/;s/\//\\\\\\\//g' | sort -u > tags.txt
while read sample accession
do
	TAX=$(sed -n '/TAXON_ID/{s/^.*>\([^<]*\)<.*$/\1/p;}' $sample.xml)
	SCN=$(sed -n '/SCIENTIFIC_NAME/{s/^.*>\([^<]*\)<.*$/\1/p;}' $sample.xml)
	while read tag
	do
		if grep -q ">$tag<" $sample.xml
		then
			sed -n '/>'"$tag"'</{n;s/^.*>\([^<]*\)<.*$/\1/p;}' $sample.xml | sed 's/^\([0-9]*\),\([0-9]*\)$/\1.\2/' | sed -e :a -e '$!N;s/\n/ | /;ta' -e 'P;D'
		else
			echo "NA"
		fi
	done < tags.txt | tr "\n" "\t" | sed 's/\t$//' | paste <(echo -e "$sample\t$accession\t$TAX\t$SCN") -
done < ${bioproject}.acc.txt | sed 's/\t/\"\t\"/g;s/$/\"/;s/\"//' | cat <(paste <(echo -e "#SampleID\tsample_accession\ttaxon_id\tscientific_name") <(tr "\n" "\t" < tags.txt | sed 's/ /_/g;s/\\//g;s/$/\n/')) - > ../${bioproject}.samples.tsv
cd ..
rm -r metadata.${bioproject}
