#!/bin/bash

set -eu

dirs=`find $SCRATCH/blood_ABF/swarm_pipe/ -name "Re_*_theta_*" -type d`

for d in $dirs; do
    if [ ! -d $d/extracts ]; then
	./run.sh $d
    fi
done
