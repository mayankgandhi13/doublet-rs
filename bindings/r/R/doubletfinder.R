# DoubletFinder-compatible interface.
#
# These functions keep DoubletFinder's names, arguments and outputs, so a
# workflow switches by replacing library(DoubletFinder) with
# library(doubletrs). Parent cells are drawn with the same R random number
# calls, so with the same seed the results are identical. The logic follows
# DoubletFinder by Chris McGinnis (CC0), with the distance-matrix and
# neighbour steps done in Rust.

#' Detect doublets in a Seurat object
#'
#' Drop-in replacement for `DoubletFinder::doubletFinder()`. Simulates
#' artificial doublets, preprocesses the merged real and artificial cells
#' with the parameters already used on `seu`, computes pANN (the proportion
#' of artificial nearest neighbours) for every cell, and calls the `nExp`
#' highest-scoring cells doublets. With the same seed, results are
#' identical to DoubletFinder's, but the dense cells x cells distance matrix
#' is never built.
#'
#' @param seu A Seurat object on which `NormalizeData()`,
#'   `FindVariableFeatures()`, `ScaleData()` and `RunPCA()` have been run
#'   (or `SCTransform()` and `RunPCA()` when `sct = TRUE`).
#' @param PCs Principal components to use, e.g. `1:10`.
#' @param pN Artificial doublets as a proportion of the merged real and
#'   artificial data.
#' @param pK Neighbourhood size as a proportion of the merged data. Choose
#'   it with [paramSweep()], [summarizeSweep()] and [find.pK()].
#' @param nExp Number of cells to call doublets, e.g. the expected doublet
#'   rate times the number of cells, adjusted with [modelHomotypic()].
#' @param reuse.pANN Name of an existing pANN metadata column. If given, no
#'   pANN is computed and only new calls for `nExp` are made.
#' @param sct Whether `seu` was preprocessed with `SCTransform()`.
#' @param annotations Not supported yet; must be `NULL`.
#' @inheritParams simulate_doublets
#'
#' @return `seu` with metadata columns `pANN_<pN>_<pK>_<nExp>` and
#'   `DF.classifications_<pN>_<pK>_<nExp>` (`"Singlet"` or `"Doublet"`),
#'   named as DoubletFinder names them.
#' @export
#' @examplesIf requireNamespace("Seurat", quietly = TRUE)
#' seu <- SeuratObject::pbmc_small
#' set.seed(1)
#' seu <- doubletFinder(seu, PCs = 1:5, pN = 0.25, pK = 0.1, nExp = 5)
#' table(seu$DF.classifications_0.25_0.1_5)
doubletFinder <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = NULL,
                          sct = FALSE, annotations = NULL,
                          threads = doubletrs_threads()) {
  check_seurat()
  if (!is.null(annotations)) {
    stop("`annotations` is not supported by doubletrs yet.", call. = FALSE)
  }
  nExp <- as_count(nExp, "nExp")
  calls_col <- paste("DF.classifications", pN, pK, nExp, sep = "_")

  if (!is.null(reuse.pANN)) {
    seu@meta.data[, calls_col] <- top_calls(seu@meta.data[, reuse.pANN], nExp)
    return(seu)
  }

  real.cells <- rownames(seu@meta.data)
  data <- seurat_counts(seu)[, real.cells]
  n_real <- length(real.cells)
  n_doublets <- round(n_real / (1 - pN) - n_real)
  pcs <- merged_pcs(seu, with_doublets(data, n_doublets, threads), PCs, sct)
  k <- round(nrow(pcs) * pK)
  pann <- unname(compute_pann(pcs, n_real, k, threads = threads))

  seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pann
  seu@meta.data[, calls_col] <- top_calls(pann, nExp)
  seu
}

#' pN-pK parameter sweep
#'
#' Drop-in replacement for `DoubletFinder::paramSweep()`. For each pN from
#' 0.05 to 0.3, simulates doublets once and computes pANN at every pK from
#' 0.0005 to 0.3. Each cell's neighbours are found once, up to the largest
#' pK, instead of sorting all distances. Datasets with more than 10,000
#' cells are down-sampled to 10,000, as in DoubletFinder.
#'
#' @inheritParams doubletFinder
#' @param num.cores Accepted for compatibility with DoubletFinder and
#'   ignored; use `threads`.
#'
#' @return A named list with one data frame (column `pANN`, one row per
#'   cell) per pN-pK pair, named `pN_<pN>_pK_<pK>`, for [summarizeSweep()].
#' @export
#' @examplesIf requireNamespace("Seurat", quietly = TRUE)
#' seu <- SeuratObject::pbmc_small
#' set.seed(1)
#' sweep <- paramSweep(seu, PCs = 1:5)
#' stats <- summarizeSweep(sweep)
#' bcmvn <- find.pK(stats)
paramSweep <- function(seu, PCs = 1:10, sct = FALSE, num.cores = 1,
                       threads = doubletrs_threads()) {
  check_seurat()
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)

  # Drop pK values giving fewer than one neighbour.
  n_cells <- nrow(seu@meta.data)
  min.cells <- round(n_cells / (1 - 0.05) - n_cells)
  pK <- pK[round(pK * min.cells) >= 1]

  counts <- seurat_counts(seu)
  if (n_cells > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:n_cells, 10000, replace = FALSE)]
    data <- counts[, real.cells]
  } else {
    real.cells <- rownames(seu@meta.data)
    data <- counts
  }
  n_real <- ncol(data)

  sweep <- list()
  for (n in seq_along(pN)) {
    n_doublets <- round(n_real / (1 - pN[n]) - n_real)
    pcs <- merged_pcs(seu, with_doublets(data, n_doublets, threads), PCs, sct)
    ks <- round(nrow(pcs) * pK)
    pann <- rs_pann_sweep(as_pcs(pcs), n_real, as.integer(ks), as_threads(threads))
    for (j in seq_along(pK)) {
      sweep[[paste("pN", pN[n], "pK", pK[j], sep = "_")]] <-
        data.frame(pANN = pann[, j], row.names = real.cells)
    }
  }
  sweep
}

#' Summarise a parameter sweep
#'
#' Drop-in replacement for `DoubletFinder::summarizeSweep()`: the bimodality
#' coefficient of each pANN distribution (from a Gaussian kernel density
#' estimate), and optionally the AUC against known doublets.
#'
#' @param sweep.list Output of [paramSweep()].
#' @param GT Whether ground-truth doublet calls are available.
#' @param GT.calls If `GT = TRUE`, a vector of `"Singlet"` / `"Doublet"`
#'   for each cell, in the order of the sweep.
#'
#' @return A data frame with columns `pN`, `pK` (factors), `BCreal`, and
#'   `AUC` when `GT = TRUE`, for [find.pK()].
#' @export
summarizeSweep <- function(sweep.list, GT = FALSE, GT.calls = NULL) {
  params <- do.call(rbind, strsplit(sub("^pN_", "", names(sweep.list)), "_pK_"))
  pN <- as.numeric(unique(params[, 1]))
  pK <- as.numeric(unique(params[, 2]))

  stats <- data.frame(
    pN = factor(rep(pN, each = length(pK))),
    pK = factor(rep(pK, length(pN)))
  )
  if (GT) {
    stats$AUC <- 0
  }
  stats$BCreal <- 0

  for (i in seq_along(sweep.list)) {
    pann <- sweep.list[[i]]$pANN
    gkde <- stats::approxfun(KernSmooth::bkde(pann, kernel = "normal"))
    x <- seq(from = min(pann), to = max(pann), length.out = length(pann))
    stats$BCreal[i] <- bimodality_coefficient(gkde(x))

    if (GT) {
      # As DoubletFinder: logistic regression on a random half of the
      # cells, AUC on the other half.
      meta <- data.frame(SinDub = factor(GT.calls, levels = c("Doublet", "Singlet")), pANN = pann)
      train <- sample(seq_len(nrow(meta)), round(nrow(meta) / 2), replace = FALSE)
      test <- seq_len(nrow(meta))[-train]
      model <- stats::glm(SinDub ~ pANN, family = stats::binomial(link = "logit"),
                          data = meta, subset = train)
      prob <- stats::predict(model, newdata = meta[test, ], type = "response")
      stats$AUC[i] <- auc(prob, meta$SinDub[test] == "Singlet")
    }
  }
  stats
}

#' Choose pK from a parameter sweep
#'
#' Drop-in replacement for `DoubletFinder::find.pK()`: the mean-variance
#' normalised bimodality coefficient (BCmvn) for each pK across all pN. The
#' pK with the highest `BCmetric` is DoubletFinder's recommended choice.
#'
#' @param sweep.stats Output of [summarizeSweep()].
#' @param plot Whether to plot BCmvn against pK, as DoubletFinder does.
#'
#' @return A data frame with columns `ParamID`, `pK`, `MeanBC`, `VarBC` and
#'   `BCmetric` (plus `MeanAUC` when the sweep was summarised with ground
#'   truth).
#' @export
find.pK <- function(sweep.stats, plot = TRUE) {
  has_auc <- "AUC" %in% colnames(sweep.stats)
  pks <- unique(sweep.stats$pK)
  bc <- lapply(pks, function(p) sweep.stats$BCreal[sweep.stats$pK == p])
  bc.mvn <- data.frame(ParamID = seq_along(pks), pK = pks)
  if (has_auc) {
    bc.mvn$MeanAUC <- vapply(pks, function(p) mean(sweep.stats$AUC[sweep.stats$pK == p]), 0)
  }
  bc.mvn$MeanBC <- vapply(bc, mean, 0)
  bc.mvn$VarBC <- vapply(bc, function(b) stats::sd(b)^2, 0)
  bc.mvn$BCmetric <- bc.mvn$MeanBC / bc.mvn$VarBC

  if (plot) {
    old <- graphics::par(mar = c(4, 4, 1, 1))
    on.exit(graphics::par(old))
    graphics::plot(bc.mvn$ParamID, bc.mvn$BCmetric, pch = 16, col = "#41b6c4", cex = 0.75,
                   xlab = "pK index", ylab = "BCmvn")
    graphics::lines(bc.mvn$ParamID, bc.mvn$BCmetric, col = "#41b6c4")
  }
  bc.mvn
}

#' Estimate the proportion of homotypic doublets
#'
#' Drop-in replacement for `DoubletFinder::modelHomotypic()`: the expected
#' fraction of doublets formed by two cells of the same type, given cell
#' type annotations. Homotypic doublets are hard to detect, so the expected
#' doublet count is often scaled by `1 - modelHomotypic(annotations)`.
#'
#' @param annotations Cell type or cluster label for each cell.
#' @return A number between 0 and 1: the sum of squared type frequencies.
#' @export
#' @examples
#' modelHomotypic(c("T", "T", "B", "B"))
modelHomotypic <- function(annotations) {
  freq <- table(annotations) / length(annotations)
  sum(freq^2)
}

# ---- Internal helpers --------------------------------------------------------

check_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Packages 'Seurat' and 'SeuratObject' are needed for this function.", call. = FALSE)
  }
}

# LayerData() reads both Seurat v5 and older objects with SeuratObject >= 5.
# (DoubletFinder's GetAssayData(slot = ) branch for old objects is defunct.)
seurat_counts <- function(seu) {
  SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")
}

# Real cells followed by `n_doublets` artificial doublets, drawn with the
# same R random number calls as DoubletFinder.
with_doublets <- function(data, n_doublets, threads) {
  data <- as_dgc(data)
  cell_a <- sample.int(ncol(data), n_doublets, replace = TRUE)
  cell_b <- sample.int(ncol(data), n_doublets, replace = TRUE)
  doublets <- build_doublets(data, cell_a, cell_b, threads)
  colnames(doublets) <- paste0("X", seq_len(n_doublets))
  cbind(data, doublets)
}

# Preprocesses merged real + artificial counts with the parameters used on
# `seu`, as DoubletFinder does, and returns cells x PCs.
merged_pcs <- function(seu, data_wdoublets, PCs, sct) {
  cmd <- seu@commands
  m <- SeuratObject::CreateSeuratObject(counts = data_wdoublets)
  if (sct) {
    m <- Seurat::SCTransform(m, verbose = FALSE)
    m <- Seurat::RunPCA(m, npcs = length(PCs), verbose = FALSE)
  } else {
    m <- Seurat::NormalizeData(m,
      normalization.method = cmd$NormalizeData.RNA@params$normalization.method,
      scale.factor = cmd$NormalizeData.RNA@params$scale.factor,
      margin = cmd$NormalizeData.RNA@params$margin, verbose = FALSE
    )
    m <- Seurat::FindVariableFeatures(m,
      selection.method = cmd$FindVariableFeatures.RNA$selection.method,
      loess.span = cmd$FindVariableFeatures.RNA$loess.span,
      clip.max = cmd$FindVariableFeatures.RNA$clip.max,
      mean.function = cmd$FindVariableFeatures.RNA$mean.function,
      dispersion.function = cmd$FindVariableFeatures.RNA$dispersion.function,
      num.bin = cmd$FindVariableFeatures.RNA$num.bin,
      binning.method = cmd$FindVariableFeatures.RNA$binning.method,
      nfeatures = cmd$FindVariableFeatures.RNA$nfeatures,
      mean.cutoff = cmd$FindVariableFeatures.RNA$mean.cutoff,
      dispersion.cutoff = cmd$FindVariableFeatures.RNA$dispersion.cutoff, verbose = FALSE
    )
    m <- Seurat::ScaleData(m,
      features = cmd$ScaleData.RNA$features,
      model.use = cmd$ScaleData.RNA$model.use,
      do.scale = cmd$ScaleData.RNA$do.scale,
      do.center = cmd$ScaleData.RNA$do.center,
      scale.max = cmd$ScaleData.RNA$scale.max,
      block.size = cmd$ScaleData.RNA$block.size,
      min.cells.to.block = cmd$ScaleData.RNA$min.cells.to.block, verbose = FALSE
    )
    m <- Seurat::RunPCA(m,
      features = cmd$ScaleData.RNA$features,
      npcs = length(PCs),
      rev.pca = cmd$RunPCA.RNA$rev.pca,
      weight.by.var = cmd$RunPCA.RNA$weight.by.var, verbose = FALSE
    )
  }
  SeuratObject::Embeddings(m, "pca")[, PCs, drop = FALSE]
}

top_calls <- function(pann, nExp) {
  if (nExp > length(pann)) {
    stop(sprintf("`nExp` (%d) is larger than the number of cells (%d).", nExp, length(pann)),
         call. = FALSE)
  }
  calls <- rep("Singlet", length(pann))
  calls[order(pann, decreasing = TRUE)[seq_len(nExp)]] <- "Doublet"
  calls
}

# Bimodality coefficient, with DoubletFinder's sample skewness and kurtosis
# written exactly as DoubletFinder writes them, so results match to the bit.
bimodality_coefficient <- function(x) {
  n <- length(x)
  S <- (1 / n) * sum((x - mean(x))^3) / (((1 / n) * sum((x - mean(x))^2))^1.5)
  G <- S * (sqrt(n * (n - 1))) / (n - 2)
  K <- (1 / n) * sum((x - mean(x))^4) / (((1 / n) * sum((x - mean(x))^2))^2) - 3
  K <- ((n - 1) * ((n + 1) * K - 3 * (n - 1)) / ((n - 2) * (n - 3))) + 3
  ((G^2) + 1) / (K + ((3 * ((n - 1)^2)) / ((n - 2) * (n - 3))))
}

# Area under the ROC curve: probability a positive scores above a negative,
# counting ties as half.
auc <- function(score, positive) {
  r <- rank(score)
  n_pos <- sum(positive)
  n_neg <- sum(!positive)
  (sum(r[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
