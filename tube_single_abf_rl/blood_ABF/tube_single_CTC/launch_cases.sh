#!/bin/bash

set -eu

mean_vel=3.5 # mm/s
Re=0.5
auf=1

Ht=0.15

for alpha in 0.0 1.6; do
    for seed in `seq 42345 42349`; do
	echo "alpha=$alpha seed=$seed"
	./run.daint.sh \
	    --Re=$Re \
	    --Ht=$Ht \
	    --auf=$auf \
            --control_param=$alpha \
            --mean_vel=$mean_vel \
	    --seed=$seed \
	    --case_name="tube_single_CTC"
    done
done
