#!/usr/bin/env bash

if [ -f path.sh ]; then . ./path.sh; fi
logdir="$DATA_ROOT/log/createBinauralFeatureDev"
mkdir -p $logdir


azimuths="`seq 0 5 90` `seq 270 5 359`"

for az in $azimuths; do
    logfile="create_dev_features_az$az.log"
    qsub -l mem=8G,rmem=4G -j y -o $logdir/$logfile \
    local/run_matlab.sh "f1_createBinauralFeatureDev($az)"
done


