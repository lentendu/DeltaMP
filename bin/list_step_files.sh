#!/bin/bash

MY_TMP_NAME=`grep "N $JOB_NAME " $EXEC/config/steps.final | awk '{print $NF}'`
MY_TMP_NAME2=${MY_TMP_NAME##*/}
find $EXEC/* -type d | sort > $EXEC/config/${MY_TMP_NAME2%%)}.dir
find $EXEC ! -type d | cat - <(echo $EXEC/config/${MY_TMP_NAME2%%)}.files) | sort > $EXEC/config/${MY_TMP_NAME2%%)}.files
