#!/bin/bash

set -eu


# default arguments

Re=0.5
RA=4
Ht=0.20
ABF_L=10
ABF_P=4
seed=12345
dry_run=false
ABF_radius=2_um
freq=100_Hz
magn_m=1e-12_N_m_per_T

Nx=1
Ny=1
Nz=1

Lx=50_um
Ly=20_um
Lz=20_um

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

usage()
{
    cat <<EOF
usage: $0
       [-h | --help] Print this error message
       [--dry-run] Print the parameters and exit
       [--Lx=<Lx>] Domain length along x (default: $Lx)
       [--Ly=<Ly>] Domain length along y (default: $Ly)
       [--Lz=<Lz>] Domain length along z (default: $Lz)
       [--Nx=<Nx>] Number of ranks along x (default: $Nx)
       [--Ny=<Ny>] Number of ranks along y (default: $Ny)
       [--Nz=<Nz>] Number of ranks along z (default: $Nz)
       [--freq=<freq>] Rotation frequency of the magnetic field (default: $freq)
       [--magn_m=<magn_m>] Magnetic moment of the ABF (default: $magn_m)
       [--Re=<Re>] Set the simulation Reynolds number (default: $Re)
       [--RA=<RA>] Set the simulation reduced radius of the RBC (default: $RA)
       [--Ht=<Ht>] Set the hematocrit (default: $Ht)
       [--seed=<seed>] Set the random seed to place the RBCs initially (default: $seed)
       [--ABF_L=<ABF_L>] The ABF tail length (arbitrary units, see in the data mesh folder.) (default: $ABF_L)
       [--ABF_P=<ABF_P>] The ABF pitch (arbitrary units, see in the data mesh folder.) (default: $ABF_P)
EOF
}

# parse optional arguments
while test $# -ne 0; do
    case "$1" in
	--dry-run)         dry_run=true;            shift ;;
	--Lx=*)            Lx="${1#*=}";            shift ;;
	--Ly=*)            Ly="${1#*=}";            shift ;;
	--Lz=*)            Lz="${1#*=}";            shift ;;
	--Nx=*)            Nx="${1#*=}";            shift ;;
	--Ny=*)            Ny="${1#*=}";            shift ;;
	--Nz=*)            Nz="${1#*=}";            shift ;;
	--freq=*)          freq="${1#*=}";          shift ;;
	--magn_m=*)        magn_m="${1#*=}";        shift ;;
	--Re=*)            Re="${1#*=}";            shift ;;
	--RA=*)            RA="${1#*=}";            shift ;;
	--Ht=*)            Ht="${1#*=}";            shift ;;
	--seed=*)          seed="${1#*=}";          shift ;;
	--ABF_L=*)         ABF_L="${1#*=}";         shift ;;
	--ABF_P=*)         ABF_P="${1#*=}";         shift ;;
	-h|--help)                                               usage; exit 0 ;;
	-*|--*) echo "Error: unsupported option $1";             usage; exit 1 ;;
	*)      echo "Error: no positional arguments required."; usage; exit 1 ;;
    esac
done

srcdir=$SCRIPTPATH
rundir=plain_single_Ht_${Ht}_freq_${freq}_ABF_radius_${ABF_radius}


ABF_data_dir=$srcdir/../data
ABF_tools_dir=$srcdir/../../ABFs/tools


origin_ABF_mesh=$ABF_data_dir/helix_head_P_${ABF_P}_L_${ABF_L}.ply

ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply
parameters="parameters.pkl"

mkdir -p $rundir
cp $srcdir/parameters.py $rundir
cp $srcdir/generate_ic.py $rundir
cp $srcdir/main.py $rundir
cd $rundir

preprocess()
{
    mir.run -n 1 ./parameters.py \
	    $origin_ABF_mesh \
	    --Re $Re \
	    --RA $RA \
	    --ABF-radius $ABF_radius \
	    --L $Lx $Ly $Lz \
	    --freq $freq \
	    --magn-m $magn_m \
	    --out-prms $parameters \
	    --out-ABF-mesh $ABF_mesh

    mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	    $ABF_mesh \
	    --out-mesh $ABF_mesh \
	    --out-coords $ABF_coords
}

num_nodes=`python -c "print($Nx * $Ny * $Nz)"`
num_ranks=`python -c "print($num_nodes * 2)"`

rbc_ic()
{
    mir.run -n $num_ranks ./generate_ic.py \
	    $parameters \
	    --ABF-coords $ABF_coords \
            --Ht $Ht \
	    --seed $seed \
	    --ranks $Nx $Ny $Nz
}

run()
{
    mir.run -n $num_ranks ./main.py \
	    $parameters \
	    --ABF-coords $ABF_coords \
	    --ranks $Nx $Ny $Nz
}

preprocess

if [ $dry_run = true ]; then
    echo "dry-run enabled, exiting."
    exit 0
fi

rbc_ic
run
