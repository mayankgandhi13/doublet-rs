#!/bin/bash
# One-time setup on Northeastern's Explorer cluster: creates the `doubletrs`
# conda environment (R, Seurat, Rust), installs DoubletFinder and the
# doubletrs R package into it, then runs the doubletrs tests.
#
#   mkdir -p /scratch/$USER/doublet-rs/logs
#   sbatch validation/hpc/setup.sbatch      # from the repo root
#SBATCH --job-name=doubletrs-setup
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/scratch/%u/doublet-rs/logs/%x-%j.out

set -euo pipefail
module load miniconda3/25.9.1
eval "$(conda shell.bash hook)"
cd "$SLURM_SUBMIT_DIR"

if conda env list | grep -qE '^doubletrs\s'; then
  conda env update -n doubletrs -f validation/hpc/environment.yml
else
  conda env create -f validation/hpc/environment.yml
fi
conda activate doubletrs
export NOT_CRAN=true

# Install into the environment's own library, not ~/R.
R_LIB=$(Rscript -e 'cat(file.path(R.home(), "library"))')
Rscript -e "
if (!requireNamespace('DoubletFinder', quietly = TRUE)) {
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade = 'never', lib = '$R_LIB')
}
"
R CMD INSTALL --library="$R_LIB" bindings/r
Rscript -e 'print(vapply(c("Seurat", "DoubletFinder", "PRROC", "doubletrs"), requireNamespace, logical(1), quietly = TRUE))'
(cd bindings/r/tests && Rscript testthat.R)
