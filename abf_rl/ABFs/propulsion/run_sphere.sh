#!/bin/bash

set -eu

data=../data
tools=../tools

origin_mesh=$data/sphere.ply

coords=coords.txt
mesh=mesh.ply

preprocess()
{
    mir.run -n 1 $tools/generate_frozen.py \
	    $origin_mesh \
	    --out-mesh $mesh \
	    --out-coords $coords \
	    --length 10
}

run()
{
    mir.run -n 2 ./main.py \
	    $mesh \
	    $coords \
	    $@
}

preprocess

run --U 1 0 0 --W 0 0 0 $@
run --U 0 0 0 --W 0.2 0 0 $@
