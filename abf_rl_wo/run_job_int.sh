#!/bin/bash
set -eu

# === Paramètres par défaut ===
Re=0.5
RA=6
Ht=0.20
auf=0
mean_vel=0_mm_per_s
control_param=0
abf_L=10
abf_P=4
seed=123
dry_run=false
L=100_um
R=10_um
abf_radius=1.63445_um
freq=1000_Hz
magn_m=1e-11_N_m_per_T
with_rbc=False
description='line - Different length scale'


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
  --L=...
  --R=...
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
        --L=*)             L="${1#*=}";             shift ;;
        --R=*)             R="${1#*=}";             shift ;;
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

rundir=$prefix/Re_${Re}_Ht_${Ht}_L_${abf_L}_P_${abf_P}_V_${mean_vel}_alpha_${control_param}_auf_${auf}_seed_${seed}_jobs
mkdir -p "$rundir"

# === Fichiers à copier ===
cp $srcdir/{parameters.py,generate_ic.py,main.py,run_job_int.sh,load_NN.py,TD3.py,utils.py} $rundir
cp $abf_tools_dir/generate_frozen.py $rundir

parameters=parameters.pkl
abf_coords=abf_coords.txt
abf_mesh=abf_mesh.ply

# === Chargement de l’environnement ===
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

cd $rundir

echo "=== PREPROCESSING ==="
srun --mpi=pmi2 -n 2 /n/home12/lchoquet/myenvMS/bin/python3 ./parameters.py \
    $origin_abf_mesh \
    --Re $Re --RA $RA --Ht $Ht \
    --abf-radius $abf_radius --L $L --R $R \
    --mean-vel $mean_vel --freq $freq \
    --magn-m $magn_m \
    --out-prms $parameters --out-abf-mesh $abf_mesh

echo "=== PARAMETERS CREATED ==="
srun --mpi=pmi2 -n 2 /n/home12/lchoquet/myenvMS/bin/python3 ./generate_frozen.py \
    $abf_mesh \
    --out-mesh $abf_mesh \
    --out-coords $abf_coords

if [ "$dry_run" = true ]; then
    echo "[DRY RUN] Prétraitement seulement, fin du script."
    exit 0
fi

echo "auf = $auf"

echo "=== GENERATING RBCs ICs ==="

srun --mpi=pmi2 -n 2 /n/home12/lchoquet/myenvMS/bin/python3 ./generate_ic.py \
    $parameters \
    --abf-coords $abf_coords \
    --Ht $Ht \
    --abf-uncenter-factor $auf \
    --seed $seed

echo "=== RUNNING SIMULATION ==="
srun --mpi=pmi2 -n 2 /n/home12/lchoquet/myenvMS/bin/python3 ./main.py \
    --parameters $parameters \
    --abf-coords $abf_coords \
    --policy $policy_dir/models_07_10_10_33/agent \
    --no-visc