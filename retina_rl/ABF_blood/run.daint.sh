#!/bin/bash

set -eu

: ${policy=policy02}

script_dir=`pwd`

data_dir=${script_dir}/../data

ABF_data_dir=$script_dir/../../microswimmers/blood_ABF/data
ABF_tools_dir=$script_dir/../../microswimmers/ABFs/tools

params="parameters.pkl"
ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply

Ht=0.25
entry_radius=10_um
entry_velocity=1_mm_per_s

ABF_L=10
ABF_P=4
ABF_radius=2_um
freq=1_kHz
magn_m=1e-11_N_m_per_T

Nx=8
Ny=8
Nz=1
gpus_per_node=1

num_gpus=`python -c "print($Nx * $Ny * $Nz)"`
num_nodes=`python -c "print(($num_gpus + $gpus_per_node - 1) // $gpus_per_node)"`
num_ranks=`python -c "print(2 * $num_gpus)"`

origin_ABF_mesh=$ABF_data_dir/helix_head_P_${ABF_P}_L_${ABF_L}.ply

run_dir=$SCRATCH/retina/ABF_blood_${Nx}_${Ny}_${Nz}_Ht_${Ht}_control_${policy}

mkdir -p $run_dir
cd $run_dir

cp -r $script_dir/../utils $run_dir/
cp $script_dir/parameters.py $run_dir/
cp $script_dir/generate_ic.py $run_dir/
cp $script_dir/main.py $run_dir/
cp $script_dir/load_NN.py $run_dir/
cp $data_dir/*.sdf $run_dir/
cp $data_dir/${policy}.json $run_dir/


sbatch_pre_file=batch_pre.sh
checkpoint=checkpoint_0

cat > $sbatch_pre_file <<EOF
#!/bin/bash -l
#
#SBATCH --job-name=retina_0
#SBATCH --time=24:00:00
#SBATCH --nodes=$num_nodes
#SBATCH --partition=normal
#SBATCH --constraint=gpu

. mir.load

mir.run -n 1 ./parameters.py \
        $origin_ABF_mesh \
	--sdf retina.sdf \
	--forces forces.sdf \
	--Ht $Ht \
	--entry-radius $entry_radius \
	--entry-velocity $entry_velocity \
        --ABF-radius $ABF_radius \
        --freq $freq \
        --magn-m $magn_m \
	--out-params $params \
        --out-ABF-mesh $ABF_mesh

mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
         $ABF_mesh \
         --out-mesh $ABF_mesh \
         --out-coords $ABF_coords

mir.run -n $num_ranks ./generate_ic.py \
	--params $params \
        --ABF-coords $ABF_coords \
	--ranks $Nx $Ny $Nz \
	--Ht $Ht \
	--scale-ini 0.3

mir.run -n $num_ranks ./main.py \
	--params $params \
        --ABF-coords $ABF_coords \
	--ranks $Nx $Ny $Nz \
	--no-visc \
        --policy ${policy}.json \
        --checkpoint-dir $checkpoint
EOF

jobid=$(sbatch $sbatch_pre_file | awk '{print $NF}')
echo "Launched job $jobid"

for run_id in $(seq 1); do
    sbatch_run_file=batch_run_${run_id}.sh
    restart=$checkpoint
    checkpoint="checkpoint_${run_id}"

    cat > $sbatch_run_file <<EOF
#!/bin/bash -l
#SBATCH --job-name=retina_${run_id}
#SBATCH --time=24:00:00
#SBATCH --nodes=$num_nodes
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --dependency=afterany:$jobid

. mir.load

mir.run -n $num_ranks ./main.py \
	--params $params \
        --ABF-coords $ABF_coords \
	--ranks $Nx $Ny $Nz \
	--no-visc \
        --policy ${policy}.json \
	--restart-dir $restart \
	--checkpoint-dir $checkpoint

EOF
    jobid=$(sbatch $sbatch_run_file | awk '{print $NF}')
    echo "Launched job $jobid"
done
