#!/bin/bash

set -eu

Re=0.1
Ma=0.02
ABF_L=27
freq="5_Hz"

ABF_length=15_um
ABF_radius=5_um
ABF_pitch=5_um
ABF_thickness=2_um

Lx="150_um"
Ly="150_um"
Lz="150_um"

case="ABF_forward"

Nx=3
Ny=3
Nz=3

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

usage()
{
    cat <<EOF
usage: $0
       [-h | --help] Print this help message and exit
       [--Re=<Re>] Set the simulation Reynolds number (default: $Re)
       [--Ma=<Ma>] Set the simulation Mach number (default: $Ma)
       [--ABF-L=<ABF_L>] Set the simulation ABF length (default: $ABF_L)
       [--ABF-length=<ABF_length>] Set the ABF length (default: $ABF_length)
       [--ABF-radius=<ABF_radius>] Set the ABF radius (default: $ABF_radius)
       [--ABF-pitch=<ABF_pitch>] Set the ABF pitch (default: $ABF_pitch)
       [--ABF-thickness=<ABF_thickness>] Set the ABF thickness (default: $ABF_thickness)
       [--freq=<f>] Set the frequency of rotation of the magnetic field (default: $freq)
       [--Lx=<Lx>] Size of the periodic domain along x (default: $Lx)
       [--Ly=<Ly>] Size of the periodic domain along y (default: $Ly)
       [--Lz=<Lz>] Size of the periodic domain along z (default: $Lz)
       [--Nx=<Nx>] Number of ranks in the x direction (default: $Nx)
       [--Ny=<Ny>] Number of ranks in the y direction (default: $Ny)
       [--Nz=<Nz>] Number of ranks in the z direction (default: $Nz)
       [--case=<case>] case name (default: $case)
EOF
}

# parse optional arguments
while test $# -ne 0; do
    case "$1" in
	--Re=*)            Re="${1#*=}"            ; shift ;;
	--Ma=*)            Ma="${1#*=}"            ; shift ;;
	--ABF-L=*)         ABF_L="${1#*=}"         ; shift ;;
	--ABF-length=*)    ABF_length="${1#*=}"    ; shift ;;
	--ABF-radius=*)    ABF_radius="${1#*=}"    ; shift ;;
	--ABF-pitch=*)     ABF_pitch="${1#*=}"     ; shift ;;
	--ABF-thickness=*) ABF_thickness="${1#*=}" ; shift ;;
	--freq=*)          freq="${1#*=}"          ; shift ;;
	--Lx=*)            Lx="${1#*=}"            ; shift ;;
	--Ly=*)            Ly="${1#*=}"            ; shift ;;
	--Lz=*)            Lz="${1#*=}"            ; shift ;;
	--Nx=*)            Nx="${1#*=}"            ; shift ;;
	--Ny=*)            Ny="${1#*=}"            ; shift ;;
	--Nz=*)            Nz="${1#*=}"            ; shift ;;
	--case=*)          case="${1#*=}"          ; shift ;;
	-h|--help)
	    usage; exit 0 ;;
	-*|--*)
	    echo "Error: unsupported option $1"
	    usage; exit 1 ;;
	*)
	    echo "Error: no positional arguments required."
	    usage; exit 1 ;;
    esac
done


src_dir=$SCRIPTPATH
run_dir=$SCRATCH/$case/freq_${freq}_L_${ABF_length}_R_${ABF_radius}_P_${ABF_pitch}_h_${ABF_thickness}/

ABF_tools_dir=$src_dir/../tools/

params="parameters.pkl"
ABF_mesh="ABF.ply"
ABF_coords="ABF_frozen.txt"

mkdir -p $run_dir

cp -r $src_dir/geometry $run_dir
cp $src_dir/parameters.py $run_dir
cp $src_dir/main.py $run_dir
cd $run_dir


NNODES=`python -c "print($Nx * $Ny * $Nz)"`
NRANKS=`python -c "print(2 * $NNODES)"`

batch_file="sbatch.txt"

cat > $batch_file <<EOF
#!/bin/bash -l
#SBATCH --job-name=L_${ABF_length}_R_${ABF_radius}_P_${ABF_pitch}_h_${ABF_thickness}
#SBATCH --time=24:00:00
#SBATCH --nodes=$NNODES
#SBATCH --partition=normal
#SBATCH --constraint=gpu

. mir.load

mir.run -n 1 ./parameters.py \
	--Re $Re \
	--Ma $Ma \
	--L $Lx $Ly $Lz \
	--ABF-L-sim $ABF_L \
	--ABF-length $ABF_length \
	--ABF-pitch $ABF_pitch \
	--ABF-radius $ABF_radius \
	--ABF-thickness $ABF_thickness \
	--frequency $freq \
	--out-params $params \
	--out-mesh $ABF_mesh

mir.run -n 1 $ABF_tools_dir/generate_frozen.py \
	$ABF_mesh \
	--out-mesh $ABF_mesh \
	--out-coords $ABF_coords

mir.run -n $NRANKS ./main.py \
	$params \
	--ranks $Nx $Ny $Nz \
	--ABF-coords $ABF_coords
EOF

sbatch $batch_file
