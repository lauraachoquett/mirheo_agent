#!/bin/bash

set -eu

data=../data
tools=../tools

origin_ABF_mesh=$data/helix_P_5.0.ply

ABF_coords=ABF_coords.txt
ABF_mesh=ABF_mesh.ply

preprocess()
{
    mir.run -n 1 $tools/generate_frozen.py \
	    $origin_ABF_mesh \
	    --out-mesh $ABF_mesh \
	    --out-coords $ABF_coords
}

run()
{
    mir.run -n 2 ./main.py \
	    $ABF_mesh \
	    $ABF_coords \
	    $@
}

preprocess

unrestricted=-0.0001
run --U $unrestricted 0 0 --W 0.1 0 0 $@

run --U 1 0 0 --W 0 0 0 $@
run --U 0 1 0 --W 0 0 0 $@
#run --U 0 0 1 --W 0 0 0 $@

run --U 0 0 0 --W 0.1 0 0 $@
run --U 0 0 0 --W 0 0.1 0 $@
#run --U 0 0 0 --W 0 0 0.1 $@
