#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ]
then
	echo "usage: ${0##*/} VERSION TARGET_DIRECTORY"
	echo "Download and unpack the MOTHUR and UTAX version of PR2 database (see https://github.com/vaulot/pr2_database) in the current directory and create symlink into the provided target directory."
	exit
fi

VERSION=$1
TARGET_DIRECTORY=$(readlink -f $2)

wget https://github.com/vaulot/pr2_database/releases/download/$VERSION/pr2_version_${VERSION}_mothur.fasta.gz
wget https://github.com/vaulot/pr2_database/releases/download/$VERSION/pr2_version_${VERSION}_mothur.tax.gz
wget https://github.com/vaulot/pr2_database/releases/download/$VERSION/pr2_version_${VERSION}_UTAX.fasta.gz

for i in pr2_version_${VERSION}_*gz; do gunzip $i ; done

cat << EOF > ${TARGET_DIRECTORY}/pr2_${VERSION}.txt
VERSION ${VERSION}
CITATION        Guillou et al., 2012
FULLCITATION    Guillou L, Bachar D, Audic S et al. (2012) The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41, D597–D604.
EOF
ln -s $PWD/pr2_version_${VERSION}_mothur.fasta ${TARGET_DIRECTORY}/pr2_${VERSION}.fasta
ln -s $PWD/pr2_version_${VERSION}_mothur.tax ${TARGET_DIRECTORY}/pr2_${VERSION}.taxonomy
vsearch --makeudb_usearch pr2_version_${VERSION}_UTAX.fasta --output pr2_version_${VERSION}_UTAX.udb
ln -s $PWD/pr2_version_${VERSION}_UTAX.udb ${TARGET_DIRECTORY}/pr2_${VERSION}.udb

echo "The database is now ready to be used in DeltaMP under the name pr2_${VERSION} ."
echo "The database containing directory to provide in DeltaMP configuration file is ${TARGET_DIRECTORY} ."
