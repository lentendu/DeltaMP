#!/bin/bash

MY_NAME=`grep "N $JOB_NAME " $EXEC/config/steps.final | awk '{print $NF}'`
find $EXEC/* -type d | sort > $EXEC/config/${MY_NAME##*/}.dir
find $EXEC ! -type d | cat - <(echo $EXEC/config/${MY_NAME##*/}.files) | sort > $EXEC/config/${MY_NAME##*/}.files
