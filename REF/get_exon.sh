#!/bin/bash

gtf="$1"
output="$2"

Rscript ./gtf_to_exons.R "$gtf" "${output}.txt.gz"