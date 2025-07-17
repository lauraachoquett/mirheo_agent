#!/bin/bash

set -eu

Re=0.1
Ma=0.02
R=5
bead_radius="5_um"
freq="5_Hz"
magn_m="1e-11_N_m_per_T"

L="100_um"

Nx=1
Ny=1
Nz=1

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

usage()
{
    cat <<EOF
usage: $0
       [-h | --help] Print this help message and exit
       [--Re=<Re>] Set the simulation Reynolds number (default: $Re)
       [--Ma=<Ma>] Set the simulation Mach number (default: $Ma)
       [--R=<R>] Set the simulation radius of the bead (default: $R)
       [--bead-radius=<bead_radius>] Set the radius of the beads (default: $bead_radius)
       [--freq=<f>] Set the frequency of rotation of the magnetic field (default: $freq)
       [--L=<L>] Size of the periodic domain (default: $L)
       [--Nx=<Nx>] Number of ranks in the x direction (default: $Nx)
       [--Ny=<Ny>] Number of ranks in the y direction (default: $Ny)
       [--Nz=<Nz>] Number of ranks in the z direction (default: $Nz)
EOF
}

# parse optional arguments
while test $# -ne 0; do
    case "$1" in
	--Re=*)            Re="${1#*=}"            ; shift ;;
	--Ma=*)            Ma="${1#*=}"            ; shift ;;
	--R=*)             R="${1#*=}"             ; shift ;;
	--bead_radius=*)   bead_radius="${1#*=}"   ; shift ;;
	--freq=*)          freq="${1#*=}"          ; shift ;;
	--L=*)             L="${1#*=}"             ; shift ;;
	--Nx=*)            Nx="${1#*=}"            ; shift ;;
	--Ny=*)            Ny="${1#*=}"            ; shift ;;
	--Nz=*)            Nz="${1#*=}"            ; shift ;;
	-h|--help)
	    usage; exit 0 ;;
	-*|--*)
	    echo "Error: unsupported option $1"
	    usage; exit 1 ;;
	*)
	    echo "Error: no positional arguments required."
	    usage; exit 1 ;;
    esac
done


src_dir=$SCRIPTPATH
run_dir=freq_${freq}_R_${bead_radius}

freeze_tools_dir=$src_dir/../../../tools/freeze/

params="parameters.pkl"
beads_coords="beads_frozen.txt"

mkdir -p $run_dir

cp $src_dir/parameters.py $run_dir
cp $src_dir/main.py $run_dir
cd $run_dir

mir.run -n 1 ./parameters.py \
	--Re $Re \
	--Ma $Ma \
	--L $L \
	--R-sim $R \
	--R $bead_radius \
	--frequency $freq \
	--magn-m $magn_m \
	--out-params $params

mir.run -n 1 $freeze_tools_dir/ellipsoid.py \
	--num-density 10 \
	--semi-axes $R $R $R \
	--niter 1000 \
	--out $beads_coords

NRANKS=`python -c "print(2 * $Nx * $Ny * $Nz)"`

mir.run -n $NRANKS ./main.py \
	$params \
	--ranks $Nx $Ny $Nz \
	--beads-coords $beads_coords
