#!/bin/bash

# 
# DeltaMP, a flexible, reproducible and resource efficient metabarcoding amplicon pipeline for HPC
# Copyright (C) 2018 Guillaume Lentendu, Tesfaye Wubet
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

if [ -z "$1" ] || [ -z "$2" ]
then
	echo "usage: ${0##*/} VERSION TARGET_DIRECTORY"
	echo "Download and unpack the MOTHUR and UTAX version of PR2 database (see https://github.com/vaulot/pr2_database) in the current directory and create symlink into the provided target directory."
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

VERSION=$1
TARGET_DIRECTORY=$(readlink -f $2)

wget https://github.com/vaulot/pr2_database/releases/download/$VERSION/pr2_version_${VERSION}_mothur.fasta.gz
wget https://github.com/vaulot/pr2_database/releases/download/$VERSION/pr2_version_${VERSION}_mothur.tax.gz

for i in pr2_version_${VERSION}_*gz; do gunzip $i ; done

cat << EOF > ${TARGET_DIRECTORY}/pr2_${VERSION}.txt
VERSION ${VERSION}
CITATION        Guillou et al., 2012
FULLCITATION    Guillou L, Bachar D, Audic S et al. (2012) The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41, D597–D604.
EOF
ln -s $PWD/pr2_version_${VERSION}_mothur.fasta ${TARGET_DIRECTORY}/pr2_${VERSION}.fasta
ln -s $PWD/pr2_version_${VERSION}_mothur.tax ${TARGET_DIRECTORY}/pr2_${VERSION}.taxonomy

if [ ! -z "$(command -v vsearch)" ]
then
	wget https://github.com/vaulot/pr2_database/releases/download/$VERSION/pr2_version_${VERSION}_UTAX.fasta.gz
	gunzip pr2_version_${VERSION}_UTAX.fasta.gz
	vsearch --makeudb_usearch pr2_version_${VERSION}_UTAX.fasta --output pr2_version_${VERSION}_UTAX.udb
	ln -s $PWD/pr2_version_${VERSION}_UTAX.udb ${TARGET_DIRECTORY}/pr2_${VERSION}.udb
fi

echo "The database is now ready to be used in DeltaMP under the name pr2_${VERSION} ."
echo "The database containing directory to provide in DeltaMP configuration file is ${TARGET_DIRECTORY} ."
