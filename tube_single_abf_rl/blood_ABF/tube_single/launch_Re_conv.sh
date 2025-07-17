#!/bin/bash

# launch simulations for several Reynolds numbers for a given ABF and Ht.

set -eu

Ht=0.15
seed=12345

for Re in 32 16 8 4 2 1 0.5 0.25 0.125; do
    echo "Ht=$Ht seed=$seed Re=$Re"
    ./run.daint.sh \
	--Re=$Re \
	--Ht=$Ht \
	--seed=$seed \
	--case_name="tube_single_Re"
done
