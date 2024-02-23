#!/bin/bash

IC_USED=$1
FOLD_NO=$2
sed -i "s/^IC_USED=.*/IC_USED=$IC_USED/" PBS/03b_fit_model.pbs
sed -i "s/^IC_USED=.*/IC_USED=$IC_USED/" PBS/03c_merge_fits.pbs
if [ -z "$FOLD_NO"]
then
  sed -i "s/^FOLD_NO=.*/FOLD_NO=\$NULL/" PBS/03b_fit_model.pbs
  sed -i "s/^FOLD_NO=.*/FOLD_NO=\$NULL/" PBS/03c_merge_fits.pbs
else
  sed -i "s/^FOLD_NO=.*/FOLD_NO=$FOLD_NO/" PBS/03b_fit_model.pbs
  sed -i "s/^FOLD_NO=.*/FOLD_NO=$FOLD_NO/" PBS/03c_merge_fits.pbs
fi
FullJobID1=`qsub PBS/03b_fit_model.pbs`
jobID1=$(echo $FullJobID1 | sed 's/[.a-z]*//g')
qsub -W depend=afterok:$jobID1 PBS/03c_merge_fits.pbs
