#!/bin/bash

ERA=$1
CHANNEL=$2
[[ -z $3 ]] || method=$3
source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh

[[ -z $method ]] && outdir=ml/out/${ERA}_${CHANNEL} || outdir=ml/out/${ERA}_${CHANNEL}_${method}
python htt-ml/application/export_keras_to_json.py  ${outdir}/dataset_config.yaml ml/templates/${ERA}_${CHANNEL}_application.yaml
