#!/bin/bash

set -eu

dirs=`find $SCRATCH/blood_ABF/swarm_pipe/ -name "Re_*_theta_*" -type d`

here=`pwd`

for d in $dirs; do
    if [ ! -f $d/movie.mp4 ]; then
	if [ -d $d/extracts ]; then
	    (
		cd $d
		~/ffmpeg-4.4-i686-static/ffmpeg -r 60 -hide_banner -loglevel warning -stats -pattern_type glob -i 'extracts/*.png' -c:v libx264 -pix_fmt yuv420p movie.mp4
		cd $here
	    )
	fi
    fi
    #
done
