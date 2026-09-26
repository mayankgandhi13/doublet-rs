#!/bin/bash
# Runs validation/compare_doubletfinder.R for every benchmark dataset and
# seed as a Slurm array (16 datasets x 3 seeds = 48 tasks).
#
#   sbatch validation/hpc/compare.sbatch    # from the repo root
#
# DoubletFinder builds a dense cells x cells distance matrix, so the largest
# datasets need tens of GB; memory is requested for the worst case.
#SBATCH --job-name=doubletrs-compare
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-47%12
#SBATCH --output=/scratch/%u/doublet-rs/logs/%x-%A_%a.out

set -euo pipefail
module load miniconda3/25.9.1
eval "$(conda shell.bash hook)"
conda activate doubletrs
cd "$SLURM_SUBMIT_DIR"

datasets=(
  cline-ch HEK-HMEC-MULTI hm-12k hm-6k HMEC-orig-MULTI HMEC-rep-MULTI
  J293t-dm mkidney-ch nuc-MULTI pbmc-1A-dm pbmc-1B-dm pbmc-1C-dm
  pbmc-2ctrl-dm pbmc-2stim-dm pbmc-ch pdx-MULTI
)
seeds=(1 2 3)

dataset=${datasets[$((SLURM_ARRAY_TASK_ID / ${#seeds[@]}))]}
seed=${seeds[$((SLURM_ARRAY_TASK_ID % ${#seeds[@]}))]}
echo "dataset=$dataset seed=$seed node=$(hostname)"

Rscript validation/compare_doubletfinder.R "validation/data/$dataset.rds" 0.09 10 "$seed"
