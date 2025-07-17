#!/bin/bash

set -eu

seed=123112
Ht=0.20
v=5_mm_per_s

for theta in `seq 0 5 85`; do
    echo "theta = $theta"

    ./run.daint.sh \
	--no-visc \
	--NX=4 \
	--NY=2 \
	--NZ=2 \
	--ABF-ic=rnd,300,$seed \
	--seed=$seed \
	--R=25_um \
	--L=200_um \
	--Ht=$Ht \
	--mean_vel=$v \
	--theta-deg=$theta

done
