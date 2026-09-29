#!/usr/bin/env Rscript
# Checks that doubletrs' DoubletFinder-compatible functions give exactly the
# same results as DoubletFinder itself, with the same seed.
#
# Usage (from the repo root):
#   Rscript validation/compare_interface.R <dataset.rds> [n_pcs] [seed]
#
# Compares doubletFinder() pANN and calls, every pANN vector from
# paramSweep(), summarizeSweep() BCreal, and the pK chosen by find.pK(),
# and reports run times for both packages.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
n_pcs <- if (length(args) >= 2) as.integer(args[2]) else 10L
seed <- if (length(args) >= 3) as.integer(args[3]) else 1L
PCs <- seq_len(n_pcs)
name <- sub("\\.rds$", "", basename(path))

data <- readRDS(path)
counts <- as(as(as(data[[1]], "dMatrix"), "generalMatrix"), "CsparseMatrix")
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
seu <- CreateSeuratObject(counts)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = max(50, n_pcs), verbose = FALSE)
cat(sprintf("%s: %d cells, PCs 1:%d, seed %d\n\n", name, ncol(seu), n_pcs, seed))

quiet <- function(expr) suppressMessages(suppressWarnings(capture.output(res <- expr)))
timed <- function(expr) {
  t <- system.time(quiet(value <- expr))[["elapsed"]]
  list(value = value, seconds = t)
}

# ---- doubletFinder() ------------------------------------------------------------

pN <- 0.25
pK <- 0.09
nExp <- round(0.075 * ncol(seu))
col_pann <- paste("pANN", pN, pK, nExp, sep = "_")
col_call <- paste("DF.classifications", pN, pK, nExp, sep = "_")

set.seed(seed)
df <- timed(DoubletFinder::doubletFinder(seu, PCs = PCs, pN = pN, pK = pK, nExp = nExp))
set.seed(seed)
rs <- timed(doubletrs::doubletFinder(seu, PCs = PCs, pN = pN, pK = pK, nExp = nExp))

pann_df <- df$value@meta.data[[col_pann]]
pann_rs <- rs$value@meta.data[[col_pann]]
finder <- data.frame(
  check = "doubletFinder()",
  identical_values = sprintf("%d / %d", sum(pann_df == pann_rs), length(pann_df)),
  max_abs_diff = max(abs(pann_df - pann_rs)),
  calls_identical = identical(df$value@meta.data[[col_call]], rs$value@meta.data[[col_call]]),
  doubletfinder_s = df$seconds,
  doubletrs_s = rs$seconds
)

# ---- paramSweep() / summarizeSweep() / find.pK() -----------------------------------

set.seed(seed)
sw_df <- timed(DoubletFinder::paramSweep(seu, PCs = PCs, sct = FALSE))
set.seed(seed)
sw_rs <- timed(doubletrs::paramSweep(seu, PCs = PCs, sct = FALSE))
stopifnot(identical(names(sw_df$value), names(sw_rs$value)))

all_df <- unlist(lapply(sw_df$value, `[[`, "pANN"))
all_rs <- unlist(lapply(sw_rs$value, `[[`, "pANN"))
sweep <- data.frame(
  check = sprintf("paramSweep() (%d pN-pK pairs)", length(sw_df$value)),
  identical_values = sprintf("%d / %d", sum(all_df == all_rs), length(all_df)),
  max_abs_diff = max(abs(all_df - all_rs)),
  calls_identical = NA,
  doubletfinder_s = sw_df$seconds,
  doubletrs_s = sw_rs$seconds
)

st_df <- DoubletFinder::summarizeSweep(sw_df$value, GT = FALSE)
st_rs <- doubletrs::summarizeSweep(sw_rs$value, GT = FALSE)
pdf(NULL)
quiet(bc_df <- DoubletFinder::find.pK(st_df))
invisible(dev.off())
bc_rs <- doubletrs::find.pK(st_rs, plot = FALSE)
best <- function(bc) as.character(bc$pK[which.max(bc$BCmetric)])
summary_row <- data.frame(
  check = "summarizeSweep() BCreal",
  identical_values = sprintf("%d / %d", sum(st_df$BCreal == st_rs$BCreal), nrow(st_df)),
  max_abs_diff = max(abs(st_df$BCreal - st_rs$BCreal)),
  calls_identical = NA,
  doubletfinder_s = NA,
  doubletrs_s = NA
)

result <- rbind(finder, sweep, summary_row)
print(result, row.names = FALSE, digits = 4)
cat(sprintf("\nfind.pK() chooses pK = %s (DoubletFinder) and pK = %s (doubletrs)\n", best(bc_df), best(bc_rs)))

dir.create("validation/results", showWarnings = FALSE)
out <- file.path("validation/results", sprintf("%s_seed%d_interface.csv", name, seed))
write.csv(cbind(dataset = name, seed = seed, result,
                chosen_pk_df = best(bc_df), chosen_pk_rs = best(bc_rs)),
          out, row.names = FALSE)
cat(sprintf("Saved %s\n", out))
