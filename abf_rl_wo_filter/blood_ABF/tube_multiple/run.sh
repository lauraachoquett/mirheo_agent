#!/bin/bash

set -eu


# default arguments

Re=0.5
RA=5
Ht=0.20
ABF_IC="rnd,20,424242"
mean_vel=0_mm_per_s
ABF_L=10
ABF_P=4
seed=12345
dry_run=false
no_visc=false
L=100_um
R=10_um
ABF_radius=1.63445_um
freq=1000_Hz
magn_m=1e-11_N_m_per_T
theta_deg=0

NX=1
NY=1
NZ=1


usage()
{
    cat <<EOF
usage: $0
       [-h | --help] Print this error message
       [--dry-run] Print the parameters and exit
       [--no-visc] Disable membrane viscosity
       [--Re=<Re>] Set the simulation Reynolds number (default: $Re)
       [--RA=<RA>] Set the simulation reduced radius of the RBC (default: $RA)
       [--Ht=<Ht>] Set the hematocrit (default: $Ht)
       [--mean_vel=<mean_vel>] Set the mean blood velocity (default: $mean_vel)
       [--freq=<f>] Rotation frequency of the magnetic field (default: $freq)
       [--magn_m=<magn_m>] Magnetic moment of the ABF (default: $magn_m)
       [--theta-deg=<theta_deg>] Angle between magnetic rotation axis and the cylinder axis (default: $theta_deg)
       [--L=<L>] Length of the pipe (default: $L)
       [--R=<R>] Radius of the pipe (default: $R)
       [--ABF-ic=<ABF_IC>] Initial positions of ABFs (default: $ABF_IC)
       [--seed=<seed>] Set the random seed to place the RBCs initially (default: $seed)
       [--ABF_L=<ABF_L>] The ABF tail length (arbitrary units, see in the data mesh folder.) (default: $ABF_L)
       [--ABF_P=<ABF_P>] The ABF pitch (arbitrary units, see in the data mesh folder.) (default: $ABF_P)
       [--NX=<NX>] number of ranks in the X direction (default: $NX)
       [--NY=<NY>] number of ranks in the Y direction (default: $NY)
       [--NZ=<NZ>] number of ranks in the Z direction (default: $NZ)
EOF
}

# parse optional arguments
while test $# -ne 0; do
    case "$1" in
	--NX=*)            NX="${1#*=}";            shift ;;
	--NY=*)            NY="${1#*=}";            shift ;;
	--NZ=*)            NZ="${1#*=}";            shift ;;
	--dry-run)         dry_run=true;            shift ;;
	--no-visc)         no_visc=true;            shift ;;
	--Re=*)            Re="${1#*=}";            shift ;;
	--RA=*)            RA="${1#*=}";            shift ;;
	--Ht=*)            Ht="${1#*=}";            shift ;;
	--mean_vel=*)      mean_vel="${1#*=}";      shift ;;
	--freq=*)          freq="${1#*=}";          shift ;;
	--magn_m=*)        magn_m="${1#*=}";        shift ;;
	--theta-deg=*)     theta_deg="${1#*=}";     shift ;;
	--L=*)             L="${1#*=}";             shift ;;
	--R=*)             R="${1#*=}";             shift ;;
	--ABF-ic=*)        ABF_IC="${1#*=}";        shift ;;
	--seed=*)          seed="${1#*=}";          shift ;;
	--ABF_L=*)         ABF_L="${1#*=}";         shift ;;
	--ABF_P=*)         ABF_P="${1#*=}";         shift ;;
	-h|--help)                                               usage; exit 0 ;;
	-*|--*) echo "Error: unsupported option $1";             usage; exit 1 ;;
	*)      echo "Error: no positional arguments required."; usage; exit 1 ;;
    esac
done

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
srcdir=$SCRIPTPATH

ABF_data_dir=$srcdir/../data
ABF_tools_dir=$srcdir/../../ABFs/tools

origin_ABF_mesh=$ABF_data_dir/helix_head_P_${ABF_P}_L_${ABF_L}.ply

rundir=Re_${Re}_Ht_${Ht}_V_${mean_vel}_theta_${theta_deg}_IC_${ABF_IC}_seed_${seed}

extra_args=""

if [ $no_visc = true ]; then
    rundir=${rundir}_novisc
    extra_args="${extra_args} --no-visc"
fi


mkdir -p $rundir
cp $srcdir/parameters.py $rundir
cp $srcdir/generate_ic.py $rundir
cp $srcdir/main.py $rundir
cp $ABF_tools_dir/generate_frozen.py $rundir
cp $srcdir/initial_pos.py $rundir
cd $rundir

ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply
parameters=parameters.pkl

nnodes=`python -c "print($NX * $NY * $NZ)"`
nranks=`python -c "print(2 * $nnodes)"`

preprocess()
{
    mir.run -n 1 ./parameters.py \
	    $origin_ABF_mesh \
	    --rbc_mesh_res 4 \
	    --Re $Re \
	    --RA $RA \
	    --ABF-radius $ABF_radius \
	    --L $L --R $R \
	    --mean-vel $mean_vel \
	    --freq $freq \
	    --magn-m $magn_m \
	    --out-prms $parameters \
	    --out-ABF-mesh $ABF_mesh

    mir.run -n 1 ./generate_frozen.py \
	    $ABF_mesh \
	    --out-mesh $ABF_mesh \
	    --out-coords $ABF_coords
}

ic()
{
    mir.run -n $nranks ./generate_ic.py \
	    $parameters \
	    --ranks $NX $NY $NZ \
	    --ABF-coords $ABF_coords \
            --Ht $Ht \
	    --ABF-ic $ABF_IC \
	    --seed $seed
}

run()
{
    mir.run -n $nranks ./main.py \
	    $parameters \
	    --theta-deg $theta_deg \
	    --ranks $NX $NY $NZ \
	    --ABF-coords $ABF_coords \
	    $extra_args
}

preprocess

if [ $dry_run = true ]; then
    echo "dry-run enabled, exiting."
    exit 0
fi

ic
run
