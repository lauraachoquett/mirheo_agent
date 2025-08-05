#!/usr/bin/sh

set -eu

. mir.load

Re=0.1
Ma=0.02
freq="5_Hz"
ABF_L=20

ABF_length=15_um
ABF_radius=2.5_um
ABF_pitch=5_um
ABF_thickness=1_um

Lx=50_um
Ly=50_um
Lz=50_um

Nx=1
Ny=1
Nz=1

for f in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
    freq="${f}_Hz"
    echo $freq
    ./run.daint.sh \
	 --Lx=$Lx --Ly=$Ly --Lz=$Lz \
	 --Nx=$Nx --Ny=$Ny --Nz=$Nz \
	 --freq=$freq --Re=$Re --Ma=$Ma \
	 --ABF-L=$ABF_L \
	 --ABF-length=$ABF_length \
	 --ABF-radius=$ABF_radius \
	 --ABF-pitch=$ABF_pitch \
	 --ABF-thickness=$ABF_thickness \
	 --case=ABF_forward_wobbling
done
