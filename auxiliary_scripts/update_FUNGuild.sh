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

if [ -z "$1" ]
then
	echo "usage: ${0##*/} TARGET_DIRECTORY"
	echo "Download and format the FUNGuild database (see https://github.com/UMNFuN/FUNGuild) in the current directory and create symlink into the provided target directory."
	exit
fi

TARGET_DIRECTORY=$(readlink -f $1)

wget http://stbates.org/funguild_db.php

grep -v "^<" funguild_db.php | sed 's/} , {/} \n {/g' | sed 's/^[^,]*, //;s/  *:  */:/g;s/  *,  */\t/g;s/\"//g;s/}$//;s/} ].*$//' | awk 'BEGIN{FS="\t"}{for(i=1;i<=NF;i++){print NR"\t"$i}}' | sed 's/:/\t/' | sort -t $'\t' -k 1,1n -k 2,2 | awk 'BEGIN{FS="\t"} {if(NR==1){prev=$1;printf "%s",$3} else {if($1==prev){printf "\t%s",$3} else {prev=$1;printf "\n%s",$3}}}END{printf "\n"}' | awk 'BEGIN{FS="\t";OFS="\t"} {print $6,$7,$9,$4,$3,$8,$2,$5,$1}' | cat <(echo -e "taxon\ttaxonomicLevel\ttrophicMode\tguild\tgrowthForm\ttrait\tconfidenceRanking\tnotes\tcitationSource") - > funguild_db.tsv
ln -s $PWD/funguild_db.tsv ${TARGET_DIRECTORY}/FUNGuild.tsv
