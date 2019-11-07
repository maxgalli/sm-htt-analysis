#!/bin/bash

# Error handling to ensure that script is executed from top-level directory of
# this repository
for DIRECTORY in shapes datacards combine plotting utils
do
    if [ ! -d "$DIRECTORY" ]; then
        echo "[FATAL] Directory $DIRECTORY not found, you are not in the top-level directory of the analysis repository?"
        exit 1
    fi
done

# Clean-up workspace
./utils/clean.sh

# Parse arguments
ERA=$1               # options: 2016, 2017
CHANNELS=${@:2}      # options: em, et, mt, tt

# Create shapes of systematics
./shapes/produce_control_shapes.sh $ERA $CHANNELS

# Apply blinding strategy
#./shapes/apply_blinding.sh $ERA

# Convert shapes to synced format
#./shapes/convert_to_synced_shapes.sh $ERA


# added by myself
exit
