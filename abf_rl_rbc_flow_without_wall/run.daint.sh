#!/bin/bash

set -eu

# default arguments

Re=0.5
Ht=0.20
RA=5
auf=0 # (relative) initial radial position of the abf
mean_vel=0_mm_per_s
control_param=0
abf_L=10
abf_P=4
seed=12345
case_name=tube_single
dry_run=false
L=100_um
R=10_um
abf_radius=1.63445_um
freq=1000_Hz
magn_m=1e-11_N_m_per_T

num_restarts=0
restart=false

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
       [--auf=<auf>] abf uncenterd factor: initial radial position of the abf, in [0,1] (default: $auf)
       [--control_param=<control_param>] Set the control parameter used to keep the abf near the center line (default: $control_param)
       [--seed=<seed>] Set the random seed to place the rbcs initially (default: $seed)
       [--abf_L=<abf_L>] The abf tail length (arbitrary units, see in the data mesh folder.) (default: $abf_L)
       [--abf_P=<abf_P>] The abf pitch (arbitrary units, see in the data mesh folder.) (default: $abf_P)
       [--case_name=<case_name>] The sub-folder to which the simulation will be saved (default: $case_name)
       [--num-restarts=<num_restarts>] The number of additional jobs of 24 hours to launch (due to time limit in slurm) (default: $num_restarts).
       [--restart] Restart from a previous simulation (default: $restart).
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
	--case_name=*)     case_name="${1#*=}";     shift ;;
	--num-restarts=*)  num_restarts="${1#*=}";  shift ;;
	--restart)         restart=true;            shift ;;
	-h|--help)                                          usage; exit 0 ;;
	-*|--*) echo "Error: unsupported option $1";        usage; exit 1 ;;
	*) echo "Error: no positional arguments required."; usage; exit 1 ;;
    esac
done

if [ $dry_run = true ]; then
    cat <<EOF
Re    = $Re
Ht    = $Ht
abf_L = $abf_L
abf_P = $abf_P
seed  = $seed
case_name = $case_name
The simulation is not launched. Exiting.
EOF
    exit 0
fi

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
srcdir=$SCRIPTPATH

abf_data_dir=$srcdir/../data
abf_tools_dir=$srcdir/../../ABFs/tools

origin_abf_mesh=$abf_data_dir/helix_head_P_${abf_P}_L_${abf_L}.ply

abf_coords=abf_coords.txt
abf_mesh=abf_mesh.ply
parameters=parameters.pkl

rundir=$SCRATCH/blood_abf/$case_name/Re_${Re}_Ht_${Ht}_L_${abf_L}_P_${abf_P}_V_${mean_vel}_alpha_${control_param}_auf_${auf}_seed_${seed}

mkdir -p $rundir
cp $srcdir/parameters.py $rundir
cp $abf_tools_dir/generate_frozen.py $rundir
cp $srcdir/generate_ic.py $rundir
cp $srcdir/main.py $rundir

cd $rundir


b0="sbatch0.sh"
checkpoint="checkpoint0"

if [ $restart = false ]; then
    sim_index=0
    cat > $b0 <<EOF
#!/bin/bash -l
#SBATCH --job-name=Ht_${Ht}
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu

. mir.load

mir.run -n 1 ./parameters.py \
	$origin_abf_mesh \
	--abf-radius $abf_radius \
	--L $L --R $R \
	--freq $freq \
	--magn-m $magn_m \
	--Re $Re \
	--mean-vel $mean_vel \
	--out-prms $parameters \
	--out-abf-mesh $abf_mesh

mir.run -n 1 ./generate_frozen.py \
	$abf_mesh \
	--out-mesh $abf_mesh \
	--out-coords $abf_coords

mir.run -n 2 ./generate_ic.py \
	$parameters \
	--abf-coords $abf_coords \
        --Ht $Ht \
	--abf-uncenter-factor $auf \
	--seed $seed

mir.run -n 2 ./main.py \
	$parameters \
	--abf-coords $abf_coords \
	--control-param $control_param \
	--checkpoint-dir $checkpoint

EOF
else
    sim_index=`find $rundir -type d -name "checkpoint*" | wc -l`
    prev_index=`python -c "print($sim_index-1)"`
    restart="checkpoint${prev_index}"
    checkpoint="checkpoint${sim_index}"
    b0="sbatch${sim_index}.sh"

    cat > $b0 <<EOF
#!/bin/bash -l
#SBATCH --job-name=Ht_${Ht}
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu

. mir.load

mir.run -n 2 ./main.py \
	$parameters \
	--abf-coords $abf_coords \
	--control-param $control_param \
	--restart-dir $restart \
	--checkpoint-dir $checkpoint

EOF
fi


jobid=`sbatch $b0 | awk '{print $NF}'`
echo "Launched job $jobid"

start_index=`python -c "print($sim_index + 1)"`
end_index=`python -c "print($sim_index + $num_restarts)"`

for i in `seq $start_index $end_index`; do
    bnext="sbatch${i}.sh"

    restart=$checkpoint
    checkpoint="checkpoint${i}"

    cat > $bnext <<EOF
#!/bin/bash -l
#!/bin/bash -l
#SBATCH --job-name=Ht_${Ht}_${i}
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --dependency=afterany:$jobid

. mir.load

mir.run -n 2 ./main.py \
	$parameters \
	--abf-coords $abf_coords \
	--control-param $control_param \
	--restart-dir $restart \
	--checkpoint-dir $checkpoint
EOF

    jobid=`sbatch $bnext | awk '{print $NF}'`
    echo "Launched job $jobid"
done
