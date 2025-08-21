#!/bin/bash
set -eu

# === Paramètres par défaut ===
Ht=0.25
entry_radius=10_um
entry_velocity=1_mm_per_s

control_param=0
abf_L=10
abf_P=4
dry_run=false

abf_radius=2_um
freq=1000_Hz
magn_m=1e-11_N_m_per_T
description='Retina'

Nx=2
Ny=2
Nz=1

((num_gpus = Nx * Ny * Nz))
((num_ranks = 2 * num_gpus))

# === Parsing des arguments ===
usage() {
    cat <<EOF
Usage: $0 [options]
Options:
  --dry-run
  --Ht=...
  --freq=...
  --magn_m=...
  --control_param=...
  --abf_L=...
  --abf_P=...
  --prefix=...
EOF
}

while test $# -ne 0; do
    case "$1" in
        --dry-run)         dry_run=true;            shift ;;
        --Ht=*)            Ht="${1#*=}";            shift ;;
        --freq=*)          freq="${1#*=}";          shift ;;
        --magn_m=*)        magn_m="${1#*=}";        shift ;;
        --control_param=*) control_param="${1#*=}"; shift ;;
        --abf_L=*)         abf_L="${1#*=}";         shift ;;
        --abf_P=*)         abf_P="${1#*=}";         shift ;;
        --prefix=*)        prefix="${1#*=}";        shift ;;
        -h|--help)         usage; exit 0 ;;
        -*|--*)            echo "Error: unsupported option $1"; usage; exit 1 ;;
        *)                 echo "Error: no positional arguments allowed."; usage; exit 1 ;;
    esac
done

# === Répertoires ===
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
srcdir=$SCRIPTPATH
prefix=$SCRATCH/koumoutsakos_lab/Lab/lchoquet

abf_data_dir=$srcdir/data
abf_tools_dir=$srcdir/ABFs/tools
policy_dir=$srcdir/../policy_file
retina_dir=$srcdir/forces
dir_path=$srcdir/data/retina2D_path_time_2025-05-28_20-40-05.npy

origin_abf_mesh=$abf_data_dir/helix_head_P_${abf_P}_L_${abf_L}.ply

rundir=$prefix/Ht_${Ht}_L_${abf_L}_P_${abf_P}_jobs_${num_gpus}_run_seas_test
mkdir -p "$rundir"

# === Fichiers à copier ===
cp $srcdir/{parameters.py,generate_ic.py,main.py,run_job.sh,load_NN.py,TD3.py,utils.py} $rundir
cp $abf_tools_dir/generate_frozen.py $rundir
# cp -r $abf_data_dir/generate_ic $rundir

parameters=parameters.pkl
abf_coords=abf_coords.txt
abf_mesh=abf_mesh.ply
checkpoint=checkpoint_0

# === Génération du script sbatch ===
batch="$rundir/sbatch.sh"

. mir.load
mamba activate /n/home12/lchoquet/myenvMS


cat > "$batch" <<EOS
#!/bin/bash
#SBATCH --partition=seas_gpu
#SBATCH --nodes=$num_gpus            # num_ranks nœuds (un par sous-domaine)
#SBATCH --ntasks-per-node=2           # 2 ranks MPI par nœud
#SBATCH --gres=gpu:1                  # 1 GPU alloué par nœud
#SBATCH -t 0-01:30
#SBATCH --job-name=gpus
#SBATCH --mem=10G
#SBATCH --output=$rundir/output.log
#SBATCH --error=$rundir/error.log




echo "=== ENV ACTIVATED ==="

cd $rundir

echo "=== PREPROCESSING ==="
srun --mpi=pmi2 -n 1 /n/home12/lchoquet/myenvMS/bin/python3 ./parameters.py \
        $origin_abf_mesh \
        --sdf $retina_dir/retina.sdf \
        --forces $retina_dir/forces.sdf \
        --Ht $Ht \
        --entry-radius $entry_radius \
        --entry-velocity $entry_velocity \
        --abf-radius $abf_radius \
        --freq $freq \
        --magn-m $magn_m \
        --out-params $parameters \
        --out-abf-mesh $abf_mesh \
        --dir_path $dir_path


echo "=== PARAMETERS CREATED ==="
srun --mpi=pmi2 -n 1 /n/home12/lchoquet/myenvMS/bin/python3 ./generate_frozen.py \
    $abf_mesh \
    --out-mesh $abf_mesh \
    --out-coords $abf_coords

echo "=== FROZEN GENERATED ==="

EOS

if [ "$dry_run" = true ]; then
    echo "[DRY RUN] Prétraitement seulement, script sbatch généré ici : $batch"
    exit 0
fi

cat >> "$batch" <<EOS

echo "=== GENERATING RBCs ICs ==="
srun --mpi=pmi2 -n $num_ranks /n/home12/lchoquet/myenvMS/bin/python3 ./generate_ic.py \
    --params "$parameters" \
    --abf-coords "$abf_coords" \
    --ranks $Nx $Ny $Nz \
    --Ht $Ht \
    --scale-ini 0.3


echo "=== RUNNING SIMULATION ==="
srun --mpi=pmi2 -n $num_ranks /n/home12/lchoquet/myenvMS/bin/python3 ./main.py \
    --params $parameters \
    --abf-coords $abf_coords \
    --policy $policy_dir/models_07_10_10_33/agent \
    --no-visc \
    --ranks $Nx $Ny $Nz \
    --checkpoint-dir $checkpoint 


EOS

# === Lancer le job ===
echo "Soumission du job SLURM depuis : $batch"
exec sbatch "$batch"

# for run_id in $(seq 1); do
#     sbatch_run_file=batch_run_${run_id}.sh
#     restart=$checkpoint
#     checkpoint="checkpoint_${run_id}"

#     cat > $sbatch_run_file <<EOF
# #/bin/bash -l
# #SBATCH --job-name=retina_${run_id}
# #SBATCH --time=24:00:00
# #SBATCH --nodes=$num_nodes
# #SBATCH --partition=normal
# #SBATCH --constraint=gpu
# #SBATCH --dependency=afterany:$jobid

# . mir.load

# mir.run -n $num_ranks ./main.py \
# 	--params $params \
#     --ABF-coords $ABF_coords \
# 	--ranks $Nx $Ny $Nz \
# 	--no-visc \
#     --policy ${policy}.json \
# 	--restart-dir $restart \
# 	--checkpoint-dir $checkpoint

# EOF
#     jobid=$(sbatch $sbatch_run_file | awk '{print $NF}')
#     echo "Launched job $jobid"
# done


