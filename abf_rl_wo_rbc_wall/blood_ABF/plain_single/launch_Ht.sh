#!/bin/bash

# launch multiple seeded simulations for several hematocrit, for a given ABF.

set -eu

ABF_L=10
ABF_P=4
Re=0.5

# freq=100_Hz
# magn_m=1e-12_N_m_per_T

freq=1000_Hz
magn_m=1e-11_N_m_per_T


for Ht in `seq 0.0 0.05 0.45`; do
    for seed in `seq 12345 12354`; do
	echo "Ht=$Ht seed=$seed Re=$Re"
	./run.daint.sh \
	    --Re=$Re \
	    --Ht=$Ht \
	    --freq=$freq \
	    --magn_m=$magn_m \
	    --ABF_L=$ABF_L \
	    --ABF_P=$ABF_P \
	    --seed=$seed \
	    --case_name="plain_single_Ht_freq_${freq}_m_${magn_m}"
    done
done
