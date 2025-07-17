#!/bin/bash

set -eu

Re=0.1
Ma=0.02
ABF_L=20
freq="5_Hz"

ABF_length=15_um
ABF_radius=2.5_um
ABF_pitch=5_um
ABF_thickness=1_um

Lx="100_um"
Ly="100_um"
Lz="100_um"

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
       [--ABF-L=<ABF_L>] Set the simulation ABF length (default: $ABF_L)
       [--ABF-length=<ABF_length>] Set the ABF length (default: $ABF_length)
       [--ABF-radius=<ABF_radius>] Set the ABF radius (default: $ABF_radius)
       [--ABF-pitch=<ABF_pitch>] Set the ABF pitch (default: $ABF_pitch)
       [--ABF-thickness=<ABF_thickness>] Set the ABF thickness (default: $ABF_thickness)
       [--freq=<f>] Set the frequency of rotation of the magnetic field (default: $freq)
       [--Lx=<Lx>] Size of the periodic domain along x (default: $Lx)
       [--Ly=<Ly>] Size of the periodic domain along y (default: $Ly)
       [--Lz=<Lz>] Size of the periodic domain along z (default: $Lz)
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
	--ABF-L=*)         ABF_L="${1#*=}"         ; shift ;;
	--ABF-length=*)    ABF_length="${1#*=}"    ; shift ;;
	--ABF-radius=*)    ABF_radius="${1#*=}"    ; shift ;;
	--ABF-pitch=*)     ABF_pitch="${1#*=}"     ; shift ;;
	--ABF-thickness=*) ABF_thickness="${1#*=}" ; shift ;;
	--freq=*)          freq="${1#*=}"          ; shift ;;
	--Lx=*)            Lx="${1#*=}"            ; shift ;;
	--Ly=*)            Ly="${1#*=}"            ; shift ;;
	--Lz=*)            Lz="${1#*=}"            ; shift ;;
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
run_dir=freq_${freq}_L_${ABF_length}_R_${ABF_radius}_P_${ABF_pitch}_h_${ABF_thickness}

ABF_tools_dir=$src_dir/../tools
echo "ABF tools dir : $ABF_tools_dir"
echo "src_dir : $src_dir"

params="parameters.pkl"
ABF_mesh="ABF_P_2.ply"
ABF_coords="ABF_frozen.txt"

mkdir -p $run_dir

cp -r $src_dir/geometry $run_dir
cp $src_dir/parameters.py $run_dir
cp $src_dir/main.py $run_dir
cd $run_dir


mir.run -n 1 python3 ./parameters.py \
	--Re $Re \
	--Ma $Ma \
	--L $Lx $Ly $Lz \
	--ABF-L-sim $ABF_L \
	--ABF-length $ABF_length \
	--ABF-pitch $ABF_pitch \
	--ABF-radius $ABF_radius \
	--ABF-thickness $ABF_thickness \
	--frequency $freq \
	--out-params $params \
	--out-mesh $ABF_mesh

echo "Re=$Re, Ma=$Ma, Lx=$Lx, Ly=$Ly, Lz=$Lz"
echo "ABF_length=$ABF_length, ABF_pitch=$ABF_pitch, ABF_radius=$ABF_radius, ABF_thickness=$ABF_thickness"
echo "params=$params, mesh=$ABF_mesh"

mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	$ABF_mesh \
	--out-mesh $ABF_mesh \
	--out-coords $ABF_coords

echo "Frozen ok"

NRANKS=`python -c "print(2 * $Nx * $Ny * $Nz)"`

mir.run -n $NRANKS ./main.py \
	$params \
	--ranks $Nx $Ny $Nz \
	--ABF-coords $ABF_coords
