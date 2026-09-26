#!/usr/bin/env Rscript
# Compare doublet-rs with DoubletFinder on a dataset with known doublets.
#
# Usage:
#   Rscript validation/compare_doubletfinder.R <dataset.rds> [pK] [n_pcs] [seed]
#
# Datasets come from the Xi & Li (2021) doublet-detection benchmark,
# https://zenodo.org/records/4062232. Each .rds is list(counts, labels),
# with labels "singlet" / "doublet".
#
# Part A (exactness): runs DoubletFinder, then repeats its steps with the
# same seed to recover the merged PCs it used, and checks that
# doubletrs::compute_pann() on those PCs reproduces DoubletFinder's pANN.
#
# Part B (accuracy): runs the doubletrs pipeline (Rust doublet simulation,
# same Seurat preprocessing, Rust kNN) and compares AUPRC / AUROC against
# the known doublets for DoubletFinder, doubletrs exact, and doubletrs HNSW.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(DoubletFinder)
  library(doubletrs)
  library(PRROC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: compare_doubletfinder.R <dataset.rds> [pK] [n_pcs] [seed]")
path <- args[1]
pK <- if (length(args) >= 2) as.numeric(args[2]) else 0.09
n_pcs <- if (length(args) >= 3) as.integer(args[3]) else 10L
seed <- if (length(args) >= 4) as.integer(args[4]) else 1L
pN <- 0.25
PCs <- seq_len(n_pcs)
name <- sub("\\.rds$", "", basename(path))

# ---- Load data -------------------------------------------------------------

data <- readRDS(path)
counts <- as(as(as(data[[1]], "dMatrix"), "generalMatrix"), "CsparseMatrix")
truth <- as.character(data[[2]]) == "doublet"
if (is.null(rownames(counts))) rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
stopifnot(length(truth) == ncol(counts))
cat(sprintf("%s: %d genes x %d cells, %d true doublets\n",
            name, nrow(counts), ncol(counts), sum(truth)))

# Standard Seurat preprocessing; DoubletFinder reuses these parameters.
seu <- CreateSeuratObject(counts)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = max(50, n_pcs), verbose = FALSE)
stopifnot(identical(colnames(seu), colnames(counts)))

raw <- LayerData(seu, assay = "RNA", layer = "counts")
n_real <- ncol(raw)
n_doublets <- round(n_real / (1 - pN) - n_real)
k <- round((n_real + n_doublets) * pK)
n_exp <- sum(truth) # threshold at the true count, so calls are comparable

# Merged real + artificial PCs, exactly as DoubletFinder computes them
# (non-SCT branch of doubletFinder(), same parameters taken from `seu`).
merged_pcs <- function(data_wdoublets) {
  cmd <- seu@commands
  m <- CreateSeuratObject(counts = data_wdoublets)
  m <- NormalizeData(m,
    normalization.method = cmd$NormalizeData.RNA@params$normalization.method,
    scale.factor = cmd$NormalizeData.RNA@params$scale.factor,
    margin = cmd$NormalizeData.RNA@params$margin, verbose = FALSE)
  m <- FindVariableFeatures(m,
    selection.method = cmd$FindVariableFeatures.RNA$selection.method,
    loess.span = cmd$FindVariableFeatures.RNA$loess.span,
    clip.max = cmd$FindVariableFeatures.RNA$clip.max,
    mean.function = cmd$FindVariableFeatures.RNA$mean.function,
    dispersion.function = cmd$FindVariableFeatures.RNA$dispersion.function,
    num.bin = cmd$FindVariableFeatures.RNA$num.bin,
    binning.method = cmd$FindVariableFeatures.RNA$binning.method,
    nfeatures = cmd$FindVariableFeatures.RNA$nfeatures,
    mean.cutoff = cmd$FindVariableFeatures.RNA$mean.cutoff,
    dispersion.cutoff = cmd$FindVariableFeatures.RNA$dispersion.cutoff, verbose = FALSE)
  m <- ScaleData(m,
    features = cmd$ScaleData.RNA$features,
    model.use = cmd$ScaleData.RNA$model.use,
    do.scale = cmd$ScaleData.RNA$do.scale,
    do.center = cmd$ScaleData.RNA$do.center,
    scale.max = cmd$ScaleData.RNA$scale.max,
    block.size = cmd$ScaleData.RNA$block.size,
    min.cells.to.block = cmd$ScaleData.RNA$min.cells.to.block, verbose = FALSE)
  m <- RunPCA(m,
    features = cmd$ScaleData.RNA$features,
    npcs = length(PCs),
    rev.pca = cmd$RunPCA.RNA$rev.pca,
    weight.by.var = cmd$RunPCA.RNA$weight.by.var, verbose = FALSE)
  m@reductions$pca@cell.embeddings[, PCs]
}

top_calls <- function(pann) {
  calls <- rep(FALSE, length(pann))
  calls[order(pann, decreasing = TRUE)[seq_len(n_exp)]] <- TRUE
  calls
}

# ---- Part A: exactness against DoubletFinder --------------------------------

cat("\n== Part A: doubletrs on DoubletFinder's own PCs ==\n")
set.seed(seed)
t_df <- system.time(
  seu_df <- suppressMessages(doubletFinder(seu, PCs = PCs, pN = pN, pK = pK, nExp = n_exp))
)[["elapsed"]]
pann_df <- seu_df@meta.data[[paste("pANN", pN, pK, n_exp, sep = "_")]]
rm(seu_df); invisible(gc())

# Same RNG draws as doubletFinder(): two sample() calls right after set.seed().
set.seed(seed)
cells1 <- sample(colnames(raw), n_doublets, replace = TRUE)
cells2 <- sample(colnames(raw), n_doublets, replace = TRUE)
df_doublets <- (raw[, cells1] + raw[, cells2]) / 2
colnames(df_doublets) <- paste0("X", seq_len(n_doublets))
pcs_df <- merged_pcs(cbind(raw, df_doublets))
rm(df_doublets)

pann_rs_on_df <- unname(compute_pann(pcs_df, n_real, k, method = "exact"))
diff <- abs(pann_rs_on_df - pann_df)
part_a <- data.frame(
  dataset = name,
  seed = seed,
  cells = n_real,
  k = k,
  cells_identical_pann = sum(diff == 0),
  max_abs_diff = max(diff),
  calls_identical = identical(top_calls(pann_rs_on_df), top_calls(pann_df))
)
print(part_a, row.names = FALSE)

# ---- Part B: accuracy of the doubletrs pipeline -----------------------------

cat("\n== Part B: accuracy against known doublets ==\n")
t_sim <- system.time(sim <- simulate_doublets(raw, n_doublets, seed = seed))[["elapsed"]]
t_pre <- system.time(pcs_rs <- merged_pcs(cbind(raw, sim)))[["elapsed"]]
t_exact <- system.time(pann_exact <- compute_pann(pcs_rs, n_real, k, "exact"))[["elapsed"]]
t_hnsw <- system.time(pann_hnsw <- compute_pann(pcs_rs, n_real, k, "hnsw"))[["elapsed"]]

score <- function(method, pann, seconds) {
  pos <- pann[truth]
  neg <- pann[!truth]
  calls <- top_calls(pann)
  data.frame(
    dataset = name,
    seed = seed,
    method = method,
    auprc = pr.curve(scores.class0 = pos, scores.class1 = neg)$auc.integral,
    auroc = roc.curve(scores.class0 = pos, scores.class1 = neg)$auc,
    precision_at_n = sum(calls & truth) / n_exp,
    total_seconds = seconds
  )
}
part_b <- rbind(
  score("DoubletFinder", pann_df, t_df),
  score("doubletrs (exact)", unname(pann_exact), t_sim + t_pre + t_exact),
  score("doubletrs (hnsw)", unname(pann_hnsw), t_sim + t_pre + t_hnsw)
)
print(part_b, row.names = FALSE, digits = 3)
cat(sprintf("\nHNSW vs exact pANN correlation: %.4f\n", cor(pann_exact, pann_hnsw)))
cat(sprintf("doubletrs time split: simulate %.2fs, Seurat preprocessing %.2fs, exact kNN %.2fs, HNSW kNN %.2fs\n",
            t_sim, t_pre, t_exact, t_hnsw))

# ---- Save ------------------------------------------------------------------

out_dir <- "validation/results" # run from the repo root
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
run_id <- sprintf("%s_seed%d", name, seed)
write.csv(part_a, file.path(out_dir, paste0(run_id, "_exactness.csv")), row.names = FALSE)
write.csv(part_b, file.path(out_dir, paste0(run_id, "_accuracy.csv")), row.names = FALSE)
saveRDS(
  list(truth = truth, pann_df = pann_df, pann_rs_on_df = pann_rs_on_df,
       pann_exact = pann_exact, pann_hnsw = pann_hnsw,
       params = list(pN = pN, pK = pK, n_pcs = n_pcs, seed = seed, k = k)),
  file.path(out_dir, paste0(run_id, "_scores.rds"))
)
cat(sprintf("\nSaved results to %s/\n", out_dir))
