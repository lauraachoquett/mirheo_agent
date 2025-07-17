#!/bin/bash

set -eu

module load daint-gpu
module load FFmpeg

ffmpeg -r 15 -hide_banner -loglevel warning -stats -i extracts/RenderView1_%06d.png -c:v mpeg4 movie.mp4
