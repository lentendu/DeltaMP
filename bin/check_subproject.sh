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

if [ ! -d $PROUT ]
then
	echo "project directory $PROUT does not exist."
	exit 1
else
	MD5CONFIG=$(md5sum $METAD)
	SUB_LIST=($(for i in $PROUT/${PROJECT}_${TECH}_${TARG}_*/config/conf*; do echo "${MD5CONFIG% *} $i" | md5sum -c | sed -n '/: OK$/s/\/[^\/]*\/[^\/]*: OK$//p' ; done))
	if [ -z $SUB_LIST ]
	then
		echo "project directory $PROUT does not contain any subproject ${PROJECT}_${TECH}_${TARG}_xxx with a configuration file identical to $METAD."
		exit 1
	elif [ ${#SUB_LIST[@]} -eq 1 ]
	then
		echo $SUB_LIST
	else
		LAST_SUB=$(for i in ${SUB_LIST[@]} ; do tree=$i/config/tree.summary; if [ -f $tree ] ; then stat -c "%s %n" $tree ; fi ; done | sort -k 1,1nr | cut -f 2 -d " ")
		echo -e "# The same configuration file was found in ${#SUB_LIST[@]} subproject with the following checkpointing relationship:\n" >&2
		column -t -s $'\t' $LAST_SUB >&2
		echo -e "\nPlease select the subproject you want to use:\n$(echo ${SUB_LIST[@]##*/} | awk '{printf "1) %s",$1;for(i=2;i<=NF;i++){printf "\\n%s) %s",i,$i}}')" >&2
		n=""
		while true
		do
		    read -p 'Select option: ' n
		    if [ "$n" -eq "$n" ] && [ "$n" -gt 0 ] && [ "$n" -le "${#SUB_LIST[@]}" ]
			then
				echo ${SUB_LIST[$((n-1))]##*/} ; break
			else
				echo "Please give a number between 1 and ${#SUB_LIST[@]} or type 'no' to cancel" >&2
			fi
		done
	fi
fi
