#!/bin/bash

# launch multiple seeded simulations for several hematocrit, for a given ABF.

set -eu

abf_L=10
abf_P=4
Re=0.5

: ${control_param=0.0}

for Ht in `seq 0.0 0.05 0.45`; do
    for seed in `seq 42345 42349`; do
	echo "Ht=$Ht seed=$seed Re=$Re"
	./run.daint.sh \
	    --Re=$Re \
	    --Ht=$Ht \
	    --abf_L=$abf_L \
	    --abf_P=$abf_P \
	    --seed=$seed \
	    --control_param=$control_param \
	    --case_name="tube_single_Ht_L_${abf_L}_P_${abf_P}_cp_${control_param}" \
	    --num-restarts=2
    done
done
