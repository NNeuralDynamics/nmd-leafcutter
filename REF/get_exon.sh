#!/bin/bash

gtf="$1"
output="$2"

Rscript /scratch/thakur.san/nmd-leafcutter/leafcutter/scripts/gtf_to_exons.R "$gtf" "${output}.txt.gz"