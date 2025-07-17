#!/bin/bash

set -eu

~/ffmpeg-4.4-i686-static/ffmpeg -r 60 -hide_banner -loglevel warning -stats -pattern_type glob -i 'extracts/*.png' -c:v libx264 -pix_fmt yuv420p movie.mp4
