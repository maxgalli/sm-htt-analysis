#!/bin/bash

CHANNEL=$1

source utils/setup_cvmfs_sft.sh

export KERAS_BACKEND=theano
export OMP_NUM_THREADS=24
export THEANO_FLAGS=gcc.cxxflags=-march=corei7

mkdir -p ml/${CHANNEL}

python htt-ml/training/keras_training.py ml/${CHANNEL}_training_config.yaml 0
python htt-ml/training/keras_training.py ml/${CHANNEL}_training_config.yaml 1
