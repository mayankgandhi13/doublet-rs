as_dgc <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (!is.matrix(x) && !is(x, "Matrix")) {
    stop("`counts` must be a matrix or Matrix object (genes x cells).", call. = FALSE)
  }
  as(as(as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

as_pcs <- function(pcs) {
  if (!is.matrix(pcs) || !is.numeric(pcs) || ncol(pcs) == 0L) {
    stop("`pcs` must be a numeric matrix (points x PCs).", call. = FALSE)
  }
  storage.mode(pcs) <- "double"
  pcs
}

as_count <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < 0 ||
    x != round(x) || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be a single non-negative whole number.", name), call. = FALSE)
  }
  as.integer(x)
}

as_threads <- function(x) {
  x <- as_count(x, "threads")
  if (x < 1L) {
    stop("`threads` must be at least 1.", call. = FALSE)
  }
  x
}
