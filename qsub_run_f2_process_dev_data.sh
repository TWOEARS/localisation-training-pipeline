#!/usr/bin/env bash

if [ -f path.sh ]; then . ./path.sh; fi
logdir="$DATA_ROOT/log/processBinauralFeatureDev"
mkdir -p $logdir

#presetList="MCT-DIFFUSE MCT-DIFFUSE-FRONT CLEAN CLEAN-FRONT"
presetList="MCT-DIFFUSE"

azRes=5
nChannels=32

for preset in $presetList; do

  for ch in `seq 1 $nChannels`; do

    logfile="process_dev_features_${preset}_${nChannels}channels_channel$ch.log"
      qsub -l mem=8G,rmem=4G -j y -o $logdir/$logfile \
        local/run_matlab.sh "f2_processBinauralFeatureDev($ch, '$preset', $azRes)"

  done
done


