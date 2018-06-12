#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
then
	echo "usage: ${0##*/} RELEASE_DATE MOTHUR_DOI TARGET_DIRECTORY [UTAX_DOI]"
	echo "Download and unpack the MOTHUR (and UTAX) version(s) of UNITE database (https://unite.ut.ee) specified by there respective DOIs (without https://doi.org/ prefix) in the current directory and create symlink into the provided target directory."
	echo "The database version is the release date version (dd.mm.yyyy)."
	exit
fi

if [ -z "$(command -v vsearch)" ]
then
	read -p "vsearch is not availble, so the UTAX version will not be downloaded and no udb format will be created for this database. Do you whish to continue (y/n)? " -n 1 -r
	echo ""
	if [[ ! $REPLY =~ ^[Yy]$ ]]
	then
	    [[ "$0" = "$BASH_SOURCE" ]] && echo "Aborted" && exit 1 || return 1
	fi
fi

RELEASE_DATE=$1
MOTHUR_DOI=$2
TARGET_DIRECTORY=$(readlink -f $3)
UTAX_DOI=$4

# Mothur format
wget -O sh_mothur_release_${RELEASE_DATE}.zip https://data.datacite.org/application/zip-1072808/${MOTHUR_DOI}
unzip -j sh_mothur_release_${RELEASE_DATE}.zip
wget -O sh_mothur_release_s_${RELEASE_DATE}.zip https://data.datacite.org/application/zip/${MOTHUR_DOI}
unzip -j sh_mothur_release_s_${RELEASE_DATE}.zip
cat > ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}_dynamic.txt << EOF
VERSION ${RELEASE_DATE}
CITATION        UNITE Community, 2017
FULLCITATION    UNITE Community (2017): UNITE mothur release. Version ${RELEASE_DATE}. UNITE Community. https://doi.org/${MOTHUR_DOI}
EOF
ln -s $PWD/UNITEv6_sh_dynamic.fasta ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}_dynamic.fasta
ln -s $PWD/UNITEv6_sh_dynamic.tax ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}_dynamic.taxonomy
ln -s $PWD/UNITEv6_sh_dynamic_s.fasta ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}_dynamic_s.fasta
ln -s $PWD/UNITEv6_sh_dynamic_s.tax ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}_dynamic_s.taxonomy

# Vsearch format
if [ ! -z "$(command -v vsearch)" ]
then
	wget -O utax_reference_dataset_${RELEASE_DATE}.zip https://data.datacite.org/application/zip/${UTAX_DOI}
	unzip -j utax_reference_dataset_${RELEASE_DATE}.zip
	cat > ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}.txt << EOF
VERSION ${RELEASE_DATE}
CITATION        UNITE Community, 2017
FULLCITATION    UNITE Community (2017): UNITE USEARCH/UTAX release. Version 01.12.2017. UNITE Community. https://doi.org/${UTAX_DOI}
EOF
	vsearch --makeudb_usearch utax_reference_dataset_${RELEASE_DATE}.fasta --output utax_reference_dataset_${RELEASE_DATE}.udb
	ln -s $PWD/utax_reference_dataset_${RELEASE_DATE}.udb ${TARGET_DIRECTORY}/unite_${RELEASE_DATE}.udb
	VSEARCH_INFO="or unite_${RELEASE_DATE} for vsearch "
fi

echo "The database is now ready to be used in DeltaMP under the name unite_${RELEASE_DATE}_dynamic for mothur ${VSEARCH_INFO}based taxonomic identification."
echo "The database containing directory to provide in DeltaMP configuration file is ${TARGET_DIRECTORY} ."
