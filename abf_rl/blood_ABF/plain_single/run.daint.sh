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
case_name=plain_blood_ABF

Nx=1
Ny=1
Nz=1

Lx=75_um
Ly=30_um
Lz=30_um

num_restarts=0
restart=false

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
       [--case_name=<case_name>] Subdirectory name (default: $case_name)
       [--num-restarts=<num_restarts>] The number of additional jobs of 24 hours to launch (due to time limit in slurm) (default: $num_restarts).
       [--restart] Restart from a previous simulation (default: $restart).
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
	--case_name=*)     case_name="${1#*=}";     shift ;;
	--num-restarts=*)  num_restarts="${1#*=}";  shift ;;
	--restart)         restart=true;            shift ;;
	-h|--help)                                               usage; exit 0 ;;
	-*|--*) echo "Error: unsupported option $1";             usage; exit 1 ;;
	*)      echo "Error: no positional arguments required."; usage; exit 1 ;;
    esac
done

if [ $dry_run = true ]; then
    cat <<EOF
Re    = $Re
Ht    = $Ht
ABF_L = $ABF_L
ABF_P = $ABF_P
seed  = $seed
case_name = $case_name
The simulation is not launched. Exiting.
EOF
    exit 0
fi

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
srcdir=$SCRIPTPATH

rbc_mesh_dir=$srcdir/../../../data/rbc_meshes
ABF_data_dir=$srcdir/../data
ABF_tools_dir=$srcdir/../../ABFs/tools

origin_ABF_mesh=$ABF_data_dir/helix_head_P_${ABF_P}_L_${ABF_L}.ply

ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply
parameters="parameters.pkl"

num_nodes=`python -c "print($Nx * $Ny * $Nz)"`
num_ranks=`python -c "print($num_nodes * 2)"`

rundir=$SCRATCH/blood_ABF/$case_name/Re_${Re}_Ht_${Ht}_L_${ABF_L}_P_${ABF_P}_freq_${freq}_ABF_radius_${ABF_radius}_seed_${seed}

mkdir -p $rundir
cp $srcdir/parameters.py $rundir
cp $srcdir/generate_ic.py $rundir
cp $srcdir/main.py $rundir

cd $rundir
sbatch=$rundir/sbatch.sh

b0=$rundir/sbatch0.sh
checkpoint_dir="checkpoint0"

if [ $restart = false ]; then
    sim_index=0
    cat > $b0 <<EOF
#!/bin/bash -l
#SBATCH --job-name=Ht_${Ht}_${freq}
#SBATCH --time=24:00:00
#SBATCH --nodes=$num_nodes
#SBATCH --partition=normal
#SBATCH --constraint=gpu

. mir.load

mir.run -n 1 ./parameters.py \
	$origin_ABF_mesh \
	--RBC-res 3 \
	--ABF-radius $ABF_radius \
	--L $Lx $Ly $Lz \
	--freq $freq \
	--Re $Re \
	--RA $RA \
	--magn-m $magn_m \
	--out-prms $parameters \
	--out-ABF-mesh $ABF_mesh

mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	$ABF_mesh \
	--out-mesh $ABF_mesh \
	--out-coords $ABF_coords

mir.run -n $num_ranks ./generate_ic.py \
	$parameters \
	--ABF-coords $ABF_coords \
        --Ht $Ht \
	--seed $seed \
	--ranks $Nx $Ny $Nz

mir.run -n $num_ranks ./main.py \
	$parameters \
	--ABF-coords $ABF_coords \
	--ranks $Nx $Ny $Nz \
	--checkpoint-dir $checkpoint_dir

EOF
else
    sim_index=`find $rundir -type d -name "checkpoint*" | wc -l`
    prev_index=`python -c "print($sim_index-1)"`
    restart_dir="checkpoint${prev_index}"
    checkpoint_dir="checkpoint${sim_index}"
    b0="sbatch${sim_index}.sh"

    cat > $b0 <<EOF
#!/bin/bash -l
#SBATCH --job-name=Ht_${Ht}_${freq}
#SBATCH --time=24:00:00
#SBATCH --nodes=$num_nodes
#SBATCH --partition=normal
#SBATCH --constraint=gpu

. mir.load

mir.run -n $num_ranks ./main.py \
	$parameters \
	--ABF-coords $ABF_coords \
	--ranks $Nx $Ny $Nz \
	--restart-dir $restart_dir \
	--checkpoint-dir $checkpoint_dir

EOF
fi


jobid=`sbatch $b0 | awk '{print $NF}'`
echo "Launched job $jobid"

start_index=`python -c "print($sim_index + 1)"`
end_index=`python -c "print($sim_index + $num_restarts)"`

for i in `seq $start_index $end_index`; do
    bnext="sbatch${i}.sh"

    restart_dir=$checkpoint_dir
    checkpoint_dir="checkpoint${i}"

    cat > $bnext <<EOF
#!/bin/bash -l
#SBATCH --job-name=Ht_${Ht}_${freq}
#SBATCH --time=24:00:00
#SBATCH --nodes=$num_nodes
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --dependency=afterany:$jobid

. mir.load

mir.run -n $num_ranks ./main.py \
	$parameters \
	--ABF-coords $ABF_coords \
	--ranks $Nx $Ny $Nz \
	--restart-dir $restart_dir \
	--checkpoint-dir $checkpoint_dir

EOF

    jobid=`sbatch $bnext | awk '{print $NF}'`
    echo "Launched job $jobid"
done
