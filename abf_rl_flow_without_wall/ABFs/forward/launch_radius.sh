#!/usr/bin/sh

set -eu

. mir.load

Re=0.1
Ma=0.02
freq="5_Hz"
ABF_L=27

ABF_length=15_um
ABF_pitch=5_um
ABF_thickness=1_um

Lx=150_um
Ly=150_um
Lz=150_um

Nx=3
Ny=3
Nz=3

for R in `seq 1.5 0.5 4.0`; do
    ABF_radius="${R}_um"
    echo $ABF_radius
    ./run.daint.sh \
	 --Lx=$Lx --Ly=$Ly --Lz=$Lz \
	 --Nx=$Nx --Ny=$Ny --Nz=$Nz \
	 --freq=$freq --Re=$Re --Ma=$Ma \
	 --ABF-L=$ABF_L \
	 --ABF-length=$ABF_length \
	 --ABF-radius=$ABF_radius \
	 --ABF-pitch=$ABF_pitch \
	 --ABF-thickness=$ABF_thickness \
	 --case=ABF_forward_radius
done
