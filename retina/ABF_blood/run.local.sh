#!/bin/bash

set -eu

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
srcdir=$SCRIPTPATH

data_dir=$srcdir/../data

ABF_data_dir=$srcdir/../../microswimmers/blood_ABF/data
ABF_tools_dir=$srcdir/../../microswimmers/ABFs/tools

parameters=parameters.pkl
ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply

Ht=0.25

ABF_L=10
ABF_P=4
ABF_radius=2_um
freq=1_kHz
magn_m=1e-11_N_m_per_T

entry_velocity=0.1_mm_per_s
entry_radius=10_um

Nx=2
Ny=1
Nz=1

origin_ABF_mesh=$ABF_data_dir/helix_head_P_${ABF_P}_L_${ABF_L}.ply

rundir=Ht_${Ht}_freq_${freq}_ABF_radius_${ABF_radius}

mkdir -p $rundir
cp $srcdir/parameters.py $rundir
cp $srcdir/generate_ic.py $rundir
cp $srcdir/main.py $rundir
cp $srcdir/load_NN.py $rundir
cp -r $srcdir/../utils $rundir
cd $rundir

num_gpus=`python -c "print($Nx * $Ny * $Nz)"`
num_ranks=`python -c "print(2 * $num_gpus)"`

params()
{
    mir.run -n 1 ./parameters.py \
	    $origin_ABF_mesh \
	    --sdf $data_dir/retina.sdf \
	    --forces $data_dir/forces.sdf \
	    --Ht $Ht \
	    --entry-velocity $entry_velocity \
	    --entry-radius $entry_radius \
	    --ABF-radius $ABF_radius \
	    --freq $freq \
	    --magn-m $magn_m \
	    --out-params $parameters \
	    --out-ABF-mesh $ABF_mesh

    mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	    $ABF_mesh \
	    --out-mesh $ABF_mesh \
	    --out-coords $ABF_coords
}

ic()
{
    mir.run -n $num_ranks ./generate_ic.py \
	    --params $parameters \
	    --ABF-coords $ABF_coords \
	    --ranks $Nx $Ny $Nz \
	    --Ht $Ht \
	    --scale-ini 0.3
}

run()
{
    mir.run -n $num_ranks ./main.py \
	    --params $parameters \
	    --ABF-coords $ABF_coords \
	    --ranks $Nx $Ny $Nz \
            --policy $data_dir/policy.json \
	    --no-visc
}

params
ic
run
