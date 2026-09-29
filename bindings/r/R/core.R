#' Simulate artificial doublets
#'
#' Builds `n` artificial doublets by averaging pairs of cells drawn
#' independently and with replacement, as DoubletFinder does. The count
#' matrix is read directly from R memory; it is not copied.
#'
#' @param counts Genes x cells count matrix, ideally a `dgCMatrix`. Other
#'   matrix types are converted to one first.
#' @param n Number of doublets to simulate.
#' @param seed Integer seed. The default draws one from R's RNG, so
#'   `set.seed()` makes results reproducible.
#' @param threads Number of threads to use. See [doubletrs_threads()].
#'
#' @return A genes x `n` `dgCMatrix` of doublets, with an attribute
#'   `"parents"`: a data frame of the two parent cell indices (1-based
#'   columns of `counts`) for each doublet.
#' @export
#' @examples
#' counts <- abs(Matrix::rsparsematrix(100, 50, density = 0.1))
#' doublets <- simulate_doublets(counts, n = 20, seed = 1)
#' dim(doublets)
#' head(attr(doublets, "parents"))
simulate_doublets <- function(counts, n, seed = NULL, threads = doubletrs_threads()) {
  counts <- as_dgc(counts)
  n <- as_count(n, "n")
  seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as_seed(seed)

  out <- rs_simulate_doublets(
    counts@i, counts@p, counts@x,
    nrow(counts), ncol(counts), n, seed, as_threads(threads)
  )

  doublets <- new(
    "dgCMatrix",
    i = out$i, p = out$p, x = out$x,
    Dim = c(nrow(counts), n),
    Dimnames = list(rownames(counts), paste0("doublet_", seq_len(n)))
  )
  attr(doublets, "parents") <- data.frame(cell_a = out$cell_a, cell_b = out$cell_b)
  doublets
}

#' Find nearest neighbours of real cells
#'
#' Exact k-nearest-neighbour search in PC space, giving the same neighbours
#' as DoubletFinder.
#'
#' @param pcs Numeric matrix, points x PCs, with the real cells in the
#'   first `n_real` rows and artificial doublets after them.
#' @param n_real Number of real cells.
#' @param k Number of neighbours per cell (excluding the cell itself).
#' @inheritParams simulate_doublets
#'
#' @return Integer matrix, `n_real` x `k`: 1-based row indices of each real
#'   cell's neighbours in `pcs`, nearest first.
#' @export
#' @examples
#' pcs <- matrix(rnorm(200), ncol = 2)
#' nn <- find_neighbors(pcs, n_real = 80, k = 5)
#' dim(nn)
find_neighbors <- function(pcs, n_real, k, threads = doubletrs_threads()) {
  pcs <- as_pcs(pcs)
  n_real <- as_count(n_real, "n_real")
  nn <- rs_find_neighbors(pcs, n_real, as_count(k, "k"), as_threads(threads))
  rownames(nn) <- rownames(pcs)[seq_len(n_real)]
  nn
}

#' Compute pANN scores
#'
#' The proportion of artificial nearest neighbours: for each real cell, the
#' fraction of its `k` nearest neighbours that are artificial doublets
#' (rows after `n_real` in `pcs`).
#'
#' @inheritParams find_neighbors
#' @return Numeric vector of length `n_real`, named by `rownames(pcs)`.
#' @export
#' @examples
#' pcs <- matrix(rnorm(200), ncol = 2)
#' pann <- compute_pann(pcs, n_real = 80, k = 5)
#' summary(pann)
compute_pann <- function(pcs, n_real, k, threads = doubletrs_threads()) {
  pcs <- as_pcs(pcs)
  n_real <- as_count(n_real, "n_real")
  pann <- rs_compute_pann(pcs, n_real, as_count(k, "k"), as_threads(threads))
  names(pann) <- rownames(pcs)[seq_len(n_real)]
  pann
}

#' Number of threads used by doubletrs
#'
#' The default for every `threads` argument. It is
#' `getOption("doubletrs.threads")` if set, otherwise the number of CPU
#' cores. During `R CMD check --as-cran` it is at most 2, as CRAN policy
#' requires.
#'
#' @return A single positive integer.
#' @export
#' @examples
#' doubletrs_threads()
#' old <- options(doubletrs.threads = 1)
#' doubletrs_threads()
#' options(old)
doubletrs_threads <- function() {
  opt <- getOption("doubletrs.threads")
  if (!is.null(opt)) {
    return(as_threads(opt))
  }
  n <- parallel::detectCores()
  if (is.na(n) || n < 1L) {
    n <- 1L
  }
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && !identical(tolower(limit), "false")) {
    n <- min(n, 2L)
  }
  as.integer(n)
}
