#!/bin/bash

set -eu

basedir=$1; shift
dst=$basedir/extracts

sbatch<<EOF
#!/bin/bash -l
#SBATCH --job-name="render_swarm"
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint=gpu
#SBATCH --account=s929
#========================================
# load modules
module load daint-gpu
module load ParaView

srun -n 1 -N 1 --cpu_bind=sockets pvbatch scene.py $basedir $dst

EOF
