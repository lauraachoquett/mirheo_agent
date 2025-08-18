#!/bin/bash
set -eu

# === Paramètres par défaut ===
Re=0.5
RA=6
Ht=0.20
auf=0 # (relative) initial radial position of the abf
mean_vel=1_mm_per_s

control_param=0
abf_L=10
abf_P=4
seed=444441
dry_run=false


Lx=100_um
Ly=30_um
Lz=30_um
abf_radius=1.63445_um
freq=1000_Hz
magn_m=1.6e-11_N_m_per_T
description='Helix - Test on FASRC with 2 GPUs - With Policy and Action - WORKING ! - Transition between ranks'

Nx=2
Ny=1
Nz=1

((num_gpus = Nx * Ny * Nz))
((num_ranks = 2 * num_gpus))

# === Parsing des arguments ===
usage() {
    cat <<EOF
Usage: $0 [options]
Options:
  --dry-run
  --Re=...
  --RA=...
  --Ht=...
  --mean_vel=...
  --freq=...
  --magn_m=...
  --auf=...
  --control_param=...
  --seed=...
  --abf_L=...
  --abf_P=...
  --prefix=...
EOF
}

while test $# -ne 0; do
    case "$1" in
        --dry-run)         dry_run=true;            shift ;;
        --Re=*)            Re="${1#*=}";            shift ;;
        --RA=*)            RA="${1#*=}";            shift ;;
        --Ht=*)            Ht="${1#*=}";            shift ;;
        --mean_vel=*)      mean_vel="${1#*=}";      shift ;;
        --freq=*)          freq="${1#*=}";          shift ;;
        --magn_m=*)        magn_m="${1#*=}";        shift ;;
        --Lx=*)            Lx="${1#*=}";            shift ;;
        --Ly=*)            Ly="${1#*=}";            shift ;;
        --Lz=*)            Lz="${1#*=}";            shift ;;
        --auf=*)           auf="${1#*=}";           shift ;;
        --control_param=*) control_param="${1#*=}"; shift ;;
        --seed=*)          seed="${1#*=}";          shift ;;
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
prefix=$srcdir

abf_data_dir=$srcdir/data
abf_tools_dir=$srcdir/ABFs/tools
policy_dir=$srcdir/../policy_file
origin_abf_mesh=$abf_data_dir/helix_head_P_${abf_P}_L_${abf_L}.ply

rundir=$prefix/Re_${Re}_Ht_${Ht}_L_${abf_L}_P_${abf_P}_V_${mean_vel}_alpha_${control_param}_auf_${auf}_seed_${seed}_jobs_${num_gpus}_test
mkdir -p "$rundir"

# === Fichiers à copier ===
cp $srcdir/{parameters.py,generate_ic.py,main.py,run_job.sh,load_NN.py,TD3.py,utils.py} $rundir
cp $abf_tools_dir/generate_frozen.py $rundir

parameters=parameters.pkl
abf_coords=abf_coords.txt
abf_mesh=abf_mesh.ply

# === Génération du script sbatch ===
batch="$rundir/sbatch.sh"

cat > "$batch" <<EOS
#!/bin/bash
#SBATCH --partition=seas_gpu
#SBATCH --nodes=$num_gpus            # num_ranks nœuds (un par sous-domaine)
#SBATCH --ntasks-per-node=2           # 2 ranks MPI par nœud
#SBATCH --gres=gpu:1                  # 1 GPU alloué par nœud
#SBATCH -t 0-02:00
#SBATCH --job-name=gpus
#SBATCH --mem=10G
#SBATCH --output=$rundir/output.log
#SBATCH --error=$rundir/error.log


. mir.load
module load cmake/3.27.5-fasrc01
module load cuda/12.2.0-fasrc01
module load Mambaforge/23.3.1-fasrc01
module load python/3.10.12-fasrc01
module load gmp/6.2.1-fasrc01
module load mpfr/4.2.0-fasrc01
module load mpc/1.3.1-fasrc01
module load gcc/9.5.0-fasrc01
module load openmpi/4.1.4-fasrc01
module load zlib/1.2.11-fasrc01
module load szip/2.1.1-fasrc01
module load hdf5/1.10.7-fasrc01

source $(mamba info --base)/etc/profile.d/conda.sh
mamba activate /n/home12/lchoquet/myenvMS

export CUDA_LAUNCH_BLOCKING=1

echo "=== ENV ACTIVATED ==="

cd $rundir

echo "=== PREPROCESSING ==="
srun --mpi=pmi2 -n 1 /n/home12/lchoquet/myenvMS/bin/python3 ./parameters.py \
	    $origin_abf_mesh \
	    --Re $Re \
	    --RA $RA \
	    --abf-radius $abf_radius \
	    --L $Lx $Ly $Lz \
	    --mean-vel $mean_vel \
	    --freq $freq \
	    --magn-m $magn_m \
	    --out-prms $parameters \
	    --out-abf-mesh $abf_mesh

echo "=== PARAMETERS CREATED ==="
srun --mpi=pmi2 -n 1 /n/home12/lchoquet/myenvMS/bin/python3 ./generate_frozen.py \
    $abf_mesh \
    --out-mesh $abf_mesh \
    --out-coords $abf_coords

EOS

if [ "$dry_run" = true ]; then
    echo "[DRY RUN] Prétraitement seulement, script sbatch généré ici : $batch"
    exit 0
fi

cat >> "$batch" <<EOS

echo "=== GENERATING RBCs ICs ==="

srun --mpi=pmi2 -n 2 /n/home12/lchoquet/myenvMS/bin/python3 ./generate_ic.py \
    $parameters \
    --abf-coords $abf_coords \
    --Ht $Ht \
    --seed $seed \

echo "=== RUNNING SIMULATION ==="
srun --mpi=pmi2 -n $num_ranks /n/home12/lchoquet/myenvMS/bin/python3 ./main.py \
    --parameters $parameters \
    --abf-coords $abf_coords \
    --policy $policy_dir/models_07_10_10_33/agent \
    --no-visc \
    --ranks $Nx $Ny $Nz

EOS

# === Lancer le job ===
echo "Soumission du job SLURM depuis : $batch"
exec sbatch "$batch"