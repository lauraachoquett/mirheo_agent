#!/bin/sh

set -eu

Ht=0.30

for seed in `seq 4200 4210`; do
    ./run.daint.sh \
	 --no-visc \
	 --seed=$seed \
	 --Ht=$Ht \
	 --control \
	 --case_name=bifurcation_control

    ./run.daint.sh \
	 --no-visc \
	 --seed=$seed \
	 --Ht=$Ht \
	 --case_name=bifurcation_no_control
done
