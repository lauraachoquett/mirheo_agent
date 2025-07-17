#!/bin/sh

set -eu

Ht=0.30
mean_vel=5_mm_per_s
seed=4567937

NX=8
NY=2
NZ=2

for theta in 0 50; do
    ./run.daint.sh \
	--NX=$NX --NY=$NY --NZ=$NZ \
	 --no-visc \
	 --seed=$seed \
	 --Ht=$Ht \
	 --mean_vel=$mean_vel \
	 --theta-deg=$theta \
	 --ABF-ic="rnd,100,$seed" \
	 --num-restarts=2
done
