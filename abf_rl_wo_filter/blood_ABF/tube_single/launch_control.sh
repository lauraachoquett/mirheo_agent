#!/bin/bash

# launch multiple seeded simulations for several control parameters.

set -eu

mean_vel=3.5_mm_per_s
Ht=0.20
Re=0.5
auf=1

for alpha in 0.0 0.01 0.025 0.05 0.1 0.2 0.4 0.8 1.6 3.2 6.4 12.8 25.6; do
    for seed in `seq 42345 42349`; do
	echo "alpha=$alpha seed=$seed"
	./run.daint.sh \
	    --Re=$Re \
	    --Ht=$Ht \
	    --auf=$auf \
            --control_param=$alpha \
            --mean_vel=$mean_vel \
	    --seed=$seed \
	    --case_name="tube_single_control_auf_${auf}_V_${mean_vel}_Ht_${Ht}"
    done
done
