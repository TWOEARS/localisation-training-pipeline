#!/bin/bash
# 

if [ -f path.sh ]; then . ./path.sh; fi
logdir="$DATA_ROOT/log/trainDNNs"
mkdir -p $logdir

#presetList="MCT-DIFFUSE MCT-DIFFUSE-FRONT CLEAN CLEAN-FRONT"
presetList="MCT-DIFFUSE"

#featureList="ild-cc itd-ild"
featureList="ild-cc"

nChannels=32
azRes=5

for preset in $presetList; do
  for feature in $featureList; do
    for ch in `seq $nChannels`; do
    #for ch in 18; do

      logfile="train_DNNs_${preset}_${feature}_${azRes}deg_${nChannels}channels_channel${ch}.log"
      qsub -l mem=8G,rmem=8G,h_rt=6:00:00 -j y -o $logdir/$logfile \
      local/run_matlab.sh "f3_trainDNNs($ch, '$preset', '$feature', $azRes)"

    done
  done
done

