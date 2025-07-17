#!/bin/bash

# launch multiple seeded simulations for several ABF pitches for a given Ht.

set -eu

ABF_L=10
Ht=0.20
Re=0.5

for ABF_P in 2 2.5 3 3.5 4 4.5 5; do
    for seed in `seq 42345 42349`; do
	echo "P=$ABF_P seed=$seed"
	./run.daint.sh \
	    --Re=$Re \
	    --Ht=$Ht \
	    --ABF_L=$ABF_L \
	    --ABF_P=$ABF_P \
	    --seed=$seed \
	    --case_name="tube_single_pitch_L_${ABF_L}_Ht_${Ht}"
    done
done
