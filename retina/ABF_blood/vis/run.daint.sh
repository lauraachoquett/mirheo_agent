#!/bin/bash

set -eu

scene_dir=$(pwd)
sim_dir=$1; shift
out_dir=$1; shift

sbatch<<EOF
#!/bin/bash -l
#SBATCH --job-name="ParaView"
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint=gpu
#SBATCH --account=s929
#========================================
# load modules
module load daint-gpu
module load ParaView

srun -n 1 -N 1 --cpu_bind=sockets pvbatch scene_bg.py $sim_dir $scene_dir $out_dir

EOF
