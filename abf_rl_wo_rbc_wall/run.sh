#!/bin/bash

set -eu


# default arguments

Re=0.5
RA=6
Ht=0.20
auf=0 # (relative) initial radial position of the abf
mean_vel=0_mm_per_s
control_param=0
abf_L=10
abf_P=4
seed=67890
dry_run=false
L=100_um
R=10_um
abf_radius=1.63445_um
freq=1000_Hz
magn_m=1e-11_N_m_per_T
description='helix - Different length scale - Wall - RBC - without flow'
usage()
{
    cat <<EOF
usage: $0
       [-h | --help] Print this error message
       [--dry-run] Print the parameters and exit
       [--Re=<Re>] Set the simulation Reynolds number (default: $Re)
       [--RA=<RA>] Set the simulation reduced radius of the rbc (default: $RA)
       [--Ht=<Ht>] Set the hematocrit (default: $Ht)
       [--mean_vel=<mean_vel>] Set the mean blood velocity (default: $mean_vel)
       [--freq=<f>] Rotation frequency of the magnetic field (default: $freq)
       [--magn_m=<magn_m>] Magnetic moment of the abf (default: $magn_m)
       [--L=<L>] Length of the pipe (default: $L)
       [--R=<R>] Radius of the pipe (default: $R)
       [--auf=<auf>] abf uncenterd factor: initai radial position of the abf, in [0,1] (default: $auf)
       [--control_param=<control_param>] Set the control parameter used to keep the abf near the center line (default: $control_param)
       [--seed=<seed>] Set the random seed to place the rbcs initially (default: $seed)
       [--abf_L=<abf_L>] The abf tail length (arbitrary units, see in the data mesh folder.) (default: $abf_L)
       [--abf_P=<abf_P>] The abf pitch (arbitrary units, see in the data mesh folder.) (default: $abf_P)
       [--prefix=<prefix>] location of output directory (default: $prefix)
EOF
}

# parse optional arguments
while test $# -ne 0; do
    case "$1" in
	--dry-run)         dry_run=true;            shift ;;
	--Re=*)            Re="${1#*=}";            shift ;;
	--RA=*)            RA="${1#*=}";            shift ;;
	--Ht=*)            Ht="${1#*=}";            shift ;;
	--mean_vel=*)      mean_vel="${1#*=}";      shift ;;
	--freq=*)          freq="${1#*=}";          shift ;;
	--magn_m=*)        magn_m="${1#*=}";        shift ;;
	--L=*)             L="${1#*=}";             shift ;;
	--R=*)             R="${1#*=}";             shift ;;
	--auf=*)           auf="${1#*=}";           shift ;;
	--control_param=*) control_param="${1#*=}"; shift ;;
	--seed=*)          seed="${1#*=}";          shift ;;
	--abf_L=*)         abf_L="${1#*=}";         shift ;;
	--abf_P=*)         abf_P="${1#*=}";         shift ;;
        --prefix=*)        prefix="${1#*=}";        shift ;;
	-h|--help)                                               usage; exit 0 ;;
	-*|--*) echo "Error: unsupported option $1";             usage; exit 1 ;;
	*)      echo "Error: no positional arguments required."; usage; exit 1 ;;
    esac
done

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
srcdir=$SCRIPTPATH
prefix=$srcdir


abf_data_dir=$srcdir/data
abf_tools_dir=$srcdir/ABFs/tools
policy_dir=$srcdir/../policy_file

origin_abf_mesh=$abf_data_dir/helix_head_P_${abf_P}_L_${abf_L}.ply

# origin_abf_mesh=$abf_data_dir/helix_head_P_2.5.ply

rundir=$prefix/Re_${Re}_Ht_${Ht}_L_${abf_L}_P_${abf_P}_V_${mean_vel}_alpha_${control_param}_auf_${auf}_seed_${seed}

mkdir -p $rundir
cp $srcdir/parameters.py $rundir
cp $srcdir/generate_ic.py $rundir
cp $srcdir/main.py $rundir
cp $srcdir/run.sh $rundir
cp $abf_tools_dir/generate_frozen.py $rundir
cp $srcdir/load_NN.py $rundir
cp $srcdir/TD3.py $rundir
cp $srcdir/utils.py $rundir

cd $rundir

abf_coords=abf_coords.txt
abf_mesh=abf_mesh.ply
parameters=parameters.pkl


preprocess()
{
    CUDA_VISIBLE_DEVICES=1 mir.run -n 2 python3 ./parameters.py \
	    $origin_abf_mesh \
	    --Re $Re \
	    --RA $RA \
	    --Ht $Ht \
	    --abf-radius $abf_radius \
	    --L $L --R $R \
	    --mean-vel $mean_vel \
	    --freq $freq \
	    --magn-m $magn_m \
	    --out-prms $parameters \
	    --out-abf-mesh $abf_mesh

    CUDA_VISIBLE_DEVICES=1 mir.run -n 2 ./generate_frozen.py \
	    $abf_mesh \
	    --out-mesh $abf_mesh \
	    --out-coords $abf_coords
}

rbc_ic()
{
    CUDA_VISIBLE_DEVICES=1 mir.run -n 2 python3 ./generate_ic.py \
	    $parameters \
	    --abf-coords $abf_coords \
            --Ht $Ht \
	    --abf-uncenter-factor $auf \
	    --seed $seed
}

run()
{
    CUDA_VISIBLE_DEVICES=1 mir.run -n 2 python3 ./main.py \
	    --parameters $parameters \
	    --abf-coords $abf_coords \
        --policy $policy_dir/models_07_10_10_33/agent \
	    --no-visc \
		--with_rbc
}

preprocess

if [ $dry_run = true ]; then
    echo "dry-run enabled, exiting."
    exit 0
fi

rbc_ic
run