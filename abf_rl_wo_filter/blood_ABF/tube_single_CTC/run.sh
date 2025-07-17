#!/bin/bash

set -eu


# default arguments

Re=0.5
Ht=0.20
auf=0 # (relative) initial radial position of the ABF
mean_vel=0 #mm/s
control_param=0
ABF_L=10
ABF_P=4
seed=12345
dry_run=false

usage()
{
    cat <<EOF
usage: $0
       [-h | --help] Print this error message
       [--dry-run] Print the parameters and exit
       [--Re=<Re>] Set the simulation Reynolds number (default: $Re)
       [--Ht=<Ht>] Set the hematocrit (default: $Ht)
       [--mean_vel=<mean_vel>] Set the mean blood velocity (in mm/s) (default: $mean_vel)
       [--auf=<auf>] ABF uncenterd factor: initai radial position of the ABF, in [0,1] (default: $auf)
       [--control_param=<control_param>] Set the control parameter used to keep the ABF near the center line (default: $control_param)
       [--seed=<seed>] Set the random seed to place the RBCs initially (default: $seed)
       [--ABF_L=<ABF_L>] The ABF tail length (arbitrary units, see in the data mesh folder.) (default: $ABF_L)
       [--ABF_P=<ABF_P>] The ABF pitch (arbitrary units, see in the data mesh folder.) (default: $ABF_P)
EOF
}

# parse optional arguments
while test $# -ne 0; do
    case "$1" in
	-h|--help)
	    usage
	    exit 0
	    ;;
	--dry-run)
	    dry_run=true
	    shift
	    ;;
	--Re=*)
	    Re="${1#*=}"
	    shift
	    ;;
	--Ht=*)
	    Ht="${1#*=}"
	    shift
	    ;;
	--mean_vel=*)
	    mean_vel="${1#*=}"
	    shift
	    ;;
	--auf=*)
	    auf="${1#*=}"
	    shift
	    ;;
	--control_param=*)
	    control_param="${1#*=}"
	    shift
	    ;;
	--seed=*)
	    seed="${1#*=}"
	    shift
	    ;;
	--ABF_L=*)
	    ABF_L="${1#*=}"
	    shift
	    ;;
	--ABF_P=*)
	    ABF_P="${1#*=}"
	    shift
	    ;;
	-*|--*)
	    echo "Error: unsupported option $1"
	    usage
	    exit 1
	    ;;
	*)
	    echo "Error: no positional arguments required."
	    usage
	    exit 1
	    ;;
    esac
done

rbc_mesh_dir=../../../data/rbc_meshes
ABF_data_dir=../data
ABF_tools_dir=../../ABFs/tools

L=100 #um
R=10 #um

ABF_radius=3.2689 # um

origin_ABF_mesh=$ABF_data_dir/helix_head_P_${ABF_P}_L_${ABF_L}.ply

ABF_freq=1000 # Hz (!)

ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply
CTC_coords=CTC_coords.txt
CTC_mesh=CTC_mesh.ply
parameters="parameters.pkl"

preprocess()
{
    mir.run -n 1 ./parameters.py \
	    $rbc_mesh_dir/ini_face_l3.off \
	    $rbc_mesh_dir/ref_face_l3.off \
	    $origin_ABF_mesh \
	    --Re $Re \
	    --ABF-radius $ABF_radius \
	    --L $L --R $R \
	    --freq $ABF_freq \
	    --out-prms $parameters \
	    --out-ABF-mesh $ABF_mesh \
    	    --out-CTC-mesh $CTC_mesh

    mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	    $ABF_mesh \
	    --out-mesh $ABF_mesh \
	    --out-coords $ABF_coords

    mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	    $CTC_mesh \
	    --out-mesh $CTC_mesh \
	    --out-coords $CTC_coords
}

rbc_ic()
{
    mir.run -n 2 ./generate_ic.py \
	    $parameters \
	    --ABF-coords $ABF_coords \
	    --CTC-coords $CTC_coords \
            --Ht $Ht \
	    --ABF-uncenter-factor $auf \
	    --seed $seed
}

run()
{
    mir.run -n 2 ./main.py \
	    $parameters \
	    --ABF-coords $ABF_coords \
    	    --CTC-coords $CTC_coords \
	    --control-param $control_param
}

preprocess

if [ $dry_run = true ]; then
    echo "dry-run enabled, exiting."
    exit 0
fi

rbc_ic
run
