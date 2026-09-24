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
#'
#' @return A genes x `n` `dgCMatrix` of doublets, with an attribute
#'   `"parents"`: a data frame of the two parent cell indices (1-based
#'   columns of `counts`) for each doublet.
#' @export
#' @examples
#' counts <- Matrix::rsparsematrix(100, 50, density = 0.1)
#' doublets <- simulate_doublets(abs(counts), n = 20)
#' dim(doublets)
simulate_doublets <- function(counts, n, seed = NULL) {
  counts <- as_dgc(counts)
  n <- as_count(n, "n")
  seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)

  out <- rs_simulate_doublets(
    counts@i, counts@p, counts@x,
    nrow(counts), ncol(counts), n, seed
  )

  doublets <- methods::new(
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
#' @param pcs Numeric matrix, points x PCs, with the real cells in the
#'   first `n_real` rows and artificial doublets after them.
#' @param n_real Number of real cells.
#' @param k Number of neighbours per cell (excluding the cell itself).
#' @param method `"exact"` for brute-force search (matches DoubletFinder) or
#'   `"hnsw"` for fast approximate search on large datasets.
#' @param hnsw Tuning parameters for `method = "hnsw"`, see [hnsw_params()].
#'
#' @return Integer matrix, `n_real` x `k`: 1-based row indices of each real
#'   cell's neighbours in `pcs`, nearest first.
#' @export
find_neighbors <- function(pcs, n_real, k, method = c("exact", "hnsw"),
                           hnsw = hnsw_params()) {
  args <- knn_args(pcs, n_real, k, method, hnsw)
  nn <- do.call(rs_find_neighbors, args)
  rownames(nn) <- rownames(pcs)[seq_len(args$n_real)]
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
compute_pann <- function(pcs, n_real, k, method = c("exact", "hnsw"),
                         hnsw = hnsw_params()) {
  args <- knn_args(pcs, n_real, k, method, hnsw)
  pann <- do.call(rs_compute_pann, args)
  names(pann) <- rownames(pcs)[seq_len(args$n_real)]
  pann
}

#' HNSW tuning parameters
#'
#' Larger values give more accurate neighbours at the cost of speed and
#' memory.
#'
#' @param max_connections Links per node per graph layer (HNSW `M`), at
#'   most 256.
#' @param ef_construction Candidate list size while building the index.
#' @param ef_search Candidate list size while searching; raised to `k + 1`
#'   if smaller.
#' @return A list of parameters for [find_neighbors()] and [compute_pann()].
#' @export
hnsw_params <- function(max_connections = 24L, ef_construction = 200L, ef_search = 128L) {
  list(
    max_connections = as_count(max_connections, "max_connections"),
    ef_construction = as_count(ef_construction, "ef_construction"),
    ef_search = as_count(ef_search, "ef_search")
  )
}

knn_args <- function(pcs, n_real, k, method, hnsw) {
  if (!is.matrix(pcs) || !is.numeric(pcs)) {
    stop("`pcs` must be a numeric matrix (points x PCs).", call. = FALSE)
  }
  storage.mode(pcs) <- "double"
  list(
    pcs = pcs,
    n_real = as_count(n_real, "n_real"),
    k = as_count(k, "k"),
    method = match.arg(method, c("exact", "hnsw")),
    max_connections = hnsw$max_connections,
    ef_construction = hnsw$ef_construction,
    ef_search = hnsw$ef_search
  )
}

as_dgc <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (!is.matrix(x) && !methods::is(x, "Matrix")) {
    stop("`counts` must be a matrix or Matrix object (genes x cells).", call. = FALSE)
  }
  methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

as_count <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x < 0 || x != round(x) || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be a single non-negative whole number.", name), call. = FALSE)
  }
  as.integer(x)
}
