#!/bin/bash
renice -n 19 -u `whoami`

ERA=$1
CHANNELS=${@:2}

BINNING=shapes/binning.yaml

source scripts/utils/clean.sh
source scripts/utils/setup_cvmfs_root.sh
source scripts/utils/setup_python.sh
source scripts/utils/setup_samples.sh $ERA

# Produce shapes
python3 main/produce_control_shapes.py \
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
