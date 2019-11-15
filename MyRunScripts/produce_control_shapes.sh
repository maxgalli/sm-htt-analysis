#!/bin/bash
renice -n 19 -u `whoami`

ERA=$1
CHANNELS=${@:2}

BINNING=shapes/binning.yaml

source MyUtils/clean.sh
source MyUtils/setup_cvmfs_sft.sh
source MyUtils/setup_python.sh
source MyUtils/setup_samples.sh $ERA

# Produce shapes
python3 MyRunScripts/produce_control_shapes_${ERA}.py \
    --directory $ARTUS_OUTPUTS \
    --datasets $KAPPA_DATABASE \
    --binning $BINNING \
    --fake-factor-friend-directory $ARTUS_FRIENDS_FAKE_FACTOR \
    --mt-friend-directory $ARTUS_FRIENDS_MT \
    --et-friend-directory $ARTUS_FRIENDS_ET \
    --tt-friend-directory $ARTUS_FRIENDS_TT \
    --em-friend-directory $ARTUS_FRIENDS_EM \
    --num-threads 30 \
    --channels ${CHANNELS}
