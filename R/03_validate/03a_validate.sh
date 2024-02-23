#!/bin/bash

IC_USED=$1
FOLD_NO=$2
sed -i "s/^IC_USED=.*/IC_USED=$IC_USED/" PBS/03a_data_idxs.pbs
if [ -z "$FOLD_NO"]
then
  sed -i "s/^FOLD_NO=.*/FOLD_NO=\$NULL/" PBS/03a_data_idxs.pbs
else
  sed -i "s/^FOLD_NO=.*/FOLD_NO=$FOLD_NO/" PBS/03a_data_idxs.pbs
fi
qsub PBS/03a_data_idxs.pbs
