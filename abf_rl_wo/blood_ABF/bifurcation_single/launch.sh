#!/bin/sh

set -eu

Ht=0.30
mean_vel=5_mm_per_s

for seed in `seq 4200 4210`; do
    ./run.daint.sh \
	 --no-visc \
	 --seed=$seed \
	 --Ht=$Ht \
	 --mean_vel=$mean_vel \
	 --control \
	 --case_name=bifurcation_control

    ./run.daint.sh \
	 --no-visc \
	 --seed=$seed \
	 --Ht=$Ht \
	 --mean_vel=$mean_vel \
	 --case_name=bifurcation_no_control
done
