#!/bin/bash

ERA=$1
STXS_FIT=$2

source utils/setup_cmssw.sh

# Set stack size to unlimited, otherwise the T2W tool throws a segfault if
# combining all eras
ulimit -s unlimited

if [ $STXS_FIT == "inclusive" ]
then
    combine -M MaxLikelihoodFit -m 125 -d ${ERA}_workspace.root \
        --robustFit 1 -n $ERA \
        --minimizerAlgoForMinos=Minuit2,Migrad
    python combine/check_mlfit.py mlfit${ERA}.root
fi

if [ $STXS_FIT == "stxs_stage0" ] || [ $STXS_FIT == "stxs_stage1" ]
then
    combineTool.py -M MultiDimFit -m 125 -d ${ERA}_workspace.root \
        --algo singles -t -1 --expectSignal 1 \
        --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP \
        --robustFit 1 -n $ERA \
        --minimizerAlgoForMinos=Minuit2,Migrad
fi
