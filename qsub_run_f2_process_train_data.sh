#!/usr/bin/env bash

if [ -f path.sh ]; then . ./path.sh; fi
logdir="$DATA_ROOT/log/processBinauralFeatureTrain"
mkdir -p $logdir

#presetList="MCT-DIFFUSE MCT-DIFFUSE-FRONT CLEAN CLEAN-FRONT"
presetList="MCT-DIFFUSE"

# Full 360-deg range
#presetList="MCT-DIFFUSE CLEAN"

# Front 180-deg range
#presetList="MCT-DIFFUSE-FRONT CLEAN-FRONT"

nChannels=32
azRes=5

for preset in $presetList; do
  for ch in `seq 1 $nChannels`; do

    logfile="process_features_${preset}_${azRes}deg_${nChannels}channels_channel$ch.log"
      qsub -l mem=8G,rmem=4G -j y -o $logdir/$logfile \
        local/run_matlab.sh "f2_processBinauralFeatureTrain($ch, '$preset', $azRes)"

  done
done


