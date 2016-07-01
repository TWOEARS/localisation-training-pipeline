#!/usr/bin/env bash

if [ -f path.sh ]; then . ./path.sh; fi
logdir="$DATA_ROOT/log/createBinauralFeatureTrain"
mkdir -p $logdir

azRes=5

# Full range
presetList="MCT-DIFFUSE"
azimuths=`seq 0 $azRes 359`

# Front hemifield only
# presetList="MCT-DIFFUSE-FRONT CLEAN-FRONT"
# azimuths="`seq 0 $azRes 90` `seq 270 $azRes 359`"

for preset in $presetList; do
  for az in $azimuths; do
    logfile="create_features_${preset}_${azRes}deg_az$az.log"
    qsub -l mem=8G,rmem=4G -j y -o $logdir/$logfile \
    local/run_matlab.sh "f1_createBinauralFeatureTrain($az, '$preset', $azRes)"
  done
done


