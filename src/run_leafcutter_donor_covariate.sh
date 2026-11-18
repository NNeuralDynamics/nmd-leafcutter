#!/bin/bash
# Run LeafCutter (donor covariate version) for every cell type in config/, one after another.
#
# Paths are relative to the project folder (nmd-leafcutter/), found from where this
# script lives, so you can run it from any folder:
#   bash src/run_all_celltypes.sh
#   bash /scratch/thakur.san/nmd-leafcutter/src/run_all_celltypes.sh
#
# Threads for the DS step: uses $THREADS if set, else SLURM_CPUS_PER_TASK, else 1.
#   THREADS=8 bash src/run_all_celltypes.sh

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base="$(cd "${script_dir}/.." && pwd)"      # this script sits in nmd-leafcutter/src/
cd "$base"

config_dir="config"
output_base="leafcutter_donor_results_v1"
exon_file="REF/gencode.v45.primary_assembly.txt.gz"
leafcutter_script="src/run_leafcutter_donor_covariate.py"   # finds ./leafcutter or ../leafcutter on its own
threads="${THREADS:-${SLURM_CPUS_PER_TASK:-1}}"

shopt -s nullglob
bam_lists=("${config_dir}"/*.txt)
if [ "${#bam_lists[@]}" -eq 0 ]; then
  echo "No BAM lists found in ${base}/${config_dir}" >&2
  exit 1
fi

for bam_list in "${bam_lists[@]}"; do
  ct=$(basename "$bam_list" .txt)
  mkdir -p "${output_base}/${ct}"

  echo "========================================="
  echo "Running LeafCutter for cell type: $ct"
  echo "  Project:  $base"
  echo "  BAM list: $bam_list"
  echo "  Threads:  $threads"
  echo "  Started:  $(date)"
  echo "========================================="

  python "$leafcutter_script" \
    --bam_file_list "$bam_list" \
    --output_dir "${output_base}/${ct}" \
    --cell_type "$ct" \
    --exon_file "$exon_file" \
    --threads "$threads"

  echo "Finished $ct at $(date)"
done

