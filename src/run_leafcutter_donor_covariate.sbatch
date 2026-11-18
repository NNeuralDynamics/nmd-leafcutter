#!/bin/bash
#SBATCH --job-name=lc_donor_cov
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/lc_donor_cov_%A_%a.out
#SBATCH --error=logs/lc_donor_cov_%A_%a.err
#
# Run from the project folder (nmd-leafcutter/). All paths below are relative to it,
# and SLURM resolves the --output/--error paths above relative to it as well.
#
#   cd /scratch/thakur.san/nmd-leafcutter
#   mkdir -p logs
#   N=$(ls config/*.txt | wc -l)
#   sbatch --array=0-$((N-1)) src/run_leafcutter_donor_covariate.sbatch

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

config_dir="config"
output_base="leafcutter_donor_results_v1"
exon_file="REF/gencode.v45.primary_assembly.txt.gz"
leafcutter_script="src/run_leafcutter_donor_covariate.py"   # finds ./leafcutter or ../leafcutter on its own

# Environment (needs both Python and R with leafcutter installed)
source /shared/EL9/explorer/anaconda3/2024.06/etc/profile.d/conda.sh
conda activate splicing

mapfile -t bam_lists < <(ls "${config_dir}"/*.txt | sort)
if [ "${SLURM_ARRAY_TASK_ID}" -ge "${#bam_lists[@]}" ]; then
  echo "Array task ${SLURM_ARRAY_TASK_ID} has no BAM list (only ${#bam_lists[@]} in ${config_dir})" >&2
  exit 1
fi
bam_list="${bam_lists[$SLURM_ARRAY_TASK_ID]}"
ct=$(basename "$bam_list" .txt)

mkdir -p "${output_base}/${ct}"

echo "========================================="
echo "Job ${SLURM_ARRAY_JOB_ID}, task ${SLURM_ARRAY_TASK_ID} on $(hostname)"
echo "Running LeafCutter for cell type: $ct"
echo "  Project:  $(pwd)"
echo "  BAM list: $bam_list"
echo "  Threads:  ${SLURM_CPUS_PER_TASK}"
echo "  Started:  $(date)"
echo "========================================="

python "$leafcutter_script" \
  --bam_file_list "$bam_list" \
  --output_dir "${output_base}/${ct}" \
  --cell_type "$ct" \
  --exon_file "$exon_file"

echo "Finished $ct at $(date)"