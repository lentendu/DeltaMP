#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]
then
	echo "usage: ${0##*/} VERSION TARGET_DIRECTORY"
	echo "Download the SSURef NR99 trunc unaligned and aligned versions of the silva database (see https://www.arb-silva.de/) in the current directory, format into the MOTHUR and UTAX formats and create symlinks into the provided target directory."
	exit
fi

if [ -z "$(command -v vsearch)" ]
then
	read -p "vsearch is not availble, so no udb format will be created for this database. Do you whish to continue (y/n)? " -n 1 -r
	echo ""
	if [[ ! $REPLY =~ ^[Yy]$ ]]
	then
	    [[ "$0" = "$BASH_SOURCE" ]] && echo "Aborted" && exit 1 || return 1
	fi
fi

VERSION=$1
TARGET_DIRECTORY=$(readlink -f $2)

wget -O - https://www.arb-silva.de/fileadmin/silva_databases/release_${VERSION}/Exports/SILVA_${VERSION}_SSURef_Nr99_tax_silva_trunc.fasta.gz | dd of=SILVA_${VERSION}_SSURef_Nr99_tax_silva_trunc.fasta.gz bs=1M
wget https://www.arb-silva.de/fileadmin/silva_databases/release_${VERSION}/Exports/SILVA_${VERSION}_SSURef_Nr99_tax_silva_trunc.fasta.gz.md5
wget -O - https://www.arb-silva.de/fileadmin/silva_databases/release_${VERSION}/Exports/SILVA_${VERSION}_SSURef_Nr99_tax_silva_full_align_trunc.fasta.gz | dd of=SILVA_${VERSION}_SSURef_Nr99_tax_silva_fill_align_trunc.fasta.gz bs=1M
wget https://www.arb-silva.de/fileadmin/silva_databases/release_${VERSION}/Exports/SILVA_${VERSION}_SSURef_Nr99_tax_silva_full_align_trunc.fasta.gz.md5

if [ ! -z $(md5sum --quiet -c SILVA_132_SSURef_Nr99_tax_silva_trunc.fasta.gz.md5) ]
then
        echo "MD5 sum mismatch, retry to download."
        exit 1
fi
if [ ! -z $(md5sum --quiet -c SILVA_132_SSURef_Nr99_tax_silva_full_align_trunc.fasta.gz.md5) ]
then
        echo "MD5 sum mismatch, retry to download."
        exit 1
fi
module load pigz
for i in *gz; do dd if=$i bs=1M | unpigz -p 1 | dd of=${i%.*} bs=1M ; done

# vsearch format
sed '/>/s/;/,/g;/>/s/ /;tax=/;/>/s/ /_/g;/>/s/[;]$/;/;/>/!s/U/T/g;/>/G' SILVA_${VERSION}_SSURef_Nr99_tax_silva_trunc.fasta | sed -e :a -e '$!N;/>/!s/\n//;ta' -e 'P;D' > silva_${VERSION}_SSURef_Nr99_UTAX.fasta
if [ ! -z "$(command -v vsearch)" ]
then
	vsearch --makeudb_usearch silva_${VERSION}_SSURef_Nr99_UTAX.fasta --output silva_${VERSION}_SSURef_Nr99_UTAX.udb
	ln -s $PWD/silva_${VERSION}_SSURef_Nr99_UTAX.udb ${TARGET_DIRECTORY}/silva_${VERSION}.udb
fi

# mothur format
awk 'BEGIN{FS=";"}{print $1}' silva_${VERSION}_SSURef_Nr99_UTAX.fasta > silva_${VERSION}_SSURef_Nr99.fasta
sed -n '/>/{s/^>//;s/;tax=/\t/;s/,/;/g;s/[;,]*$/;/;p}' silva_${VERSION}_SSURef_Nr99_UTAX.fasta > silva_${VERSION}_SSURef_Nr99.taxonomy
awk '{print $1}' SILVA_${VERSION}_SSURef_Nr99_tax_silva_full_align_trunc.fasta | sed '/>/!s/U/T/g' > silva_${VERSION}_SSURef_Nr99.align.fasta
cat << EOF > ${TARGET_DIRECTORY}/silva_${VERSION}.txt
VERSION $VERSION
CITATION        Quast et al., 2013
FULLCITATION    Quast C, Pruesse E, Yilmaz P et al. (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research, 41, D590–D596.
EOF
ln -s $PWD/silva_${VERSION}_SSURef_Nr99.fasta ${TARGET_DIRECTORY}/silva_${VERSION}.fasta
ln -s $PWD/silva_${VERSION}_SSURef_Nr99.taxonomy ${TARGET_DIRECTORY}/silva_${VERSION}.taxonomy
ln -s $PWD/silva_${VERSION}_SSURef_Nr99.align.fasta ${TARGET_DIRECTORY}/silva_${VERSION}.align.fasta

echo "The database is now ready to be used in DeltaMP under the name silva_${VERSION} ."
echo "The database containing directory to provide in DeltaMP configuration file is ${TARGET_DIRECTORY} ."
