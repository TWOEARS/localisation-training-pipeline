#!/bin/bash
# 

if [ -f path.sh ]; then . ./path.sh; fi
logdir="$DATA_ROOT/log/trainGMMs"
mkdir -p $logdir

#presetList="MCT-DIFFUSE MCT-DIFFUSE-FRONT CLEAN CLEAN-FRONT"
presetList="MCT-DIFFUSE"

#featureList="itd-ild cc-ild"
featureList="itd-ild"

nChannels=32
azRes=5

for feature in $featureList; do
  for preset in $presetList; do
    
      logfile="train_GMMs_${preset}_${feature}_${azRes}deg_${nChannels}channels.log"
      qsub -l mem=10G,rmem=8G,h_rt=8:00:00 -j y -o $logdir/$logfile \
      local/run_matlab.sh "f3_trainGMMs('$preset', '$feature', $azRes)"

  done
done

