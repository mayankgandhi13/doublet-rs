toy_counts <- function(genes = 30, cells = 40, seed = 1) {
  set.seed(seed)
  m <- Matrix::rsparsematrix(genes, cells, density = 0.3, rand.x = function(n) rpois(n, 5) + 1)
  dimnames(m) <- list(paste0("g", seq_len(genes)), paste0("c", seq_len(cells)))
  m
}

test_that("doublets are the average of their parent cells", {
  counts <- toy_counts()
  d <- simulate_doublets(counts, n = 25, seed = 7)
  parents <- attr(d, "parents")

  expect_s4_class(d, "dgCMatrix")
  expect_equal(dim(d), c(nrow(counts), 25L))
  expect_equal(rownames(d), rownames(counts))
  expect_equal(nrow(parents), 25L)
  expect_true(all(unlist(parents) >= 1 & unlist(parents) <= ncol(counts)))

  expected <- (counts[, parents$cell_a] + counts[, parents$cell_b]) / 2
  expect_equal(unname(as.matrix(d)), unname(as.matrix(expected)))
})

test_that("seeds make simulation reproducible", {
  counts <- toy_counts()
  expect_identical(simulate_doublets(counts, 10, seed = 3), simulate_doublets(counts, 10, seed = 3))
  expect_false(identical(
    attr(simulate_doublets(counts, 10, seed = 3), "parents"),
    attr(simulate_doublets(counts, 10, seed = 4), "parents")
  ))

  set.seed(99); a <- simulate_doublets(counts, 10)
  set.seed(99); b <- simulate_doublets(counts, 10)
  expect_identical(a, b)
})

test_that("dense matrices are accepted", {
  counts <- toy_counts()
  expect_equal(
    simulate_doublets(as.matrix(counts), 5, seed = 1),
    simulate_doublets(counts, 5, seed = 1)
  )
})

test_that("exact neighbours match a brute-force R computation", {
  set.seed(2)
  pcs <- matrix(rnorm(60 * 5), 60, 5)
  n_real <- 40
  k <- 6
  nn <- find_neighbors(pcs, n_real, k, method = "exact")

  d <- as.matrix(dist(pcs))
  expected <- t(sapply(seq_len(n_real), function(i) order(d[i, ])[2:(k + 1)]))
  expect_equal(unname(nn), expected)
})

test_that("pANN is the share of neighbours that are artificial", {
  set.seed(3)
  pcs <- matrix(rnorm(80 * 4), 80, 4)
  rownames(pcs) <- paste0("p", 1:80)
  nn <- find_neighbors(pcs, 60, 10)
  pann <- compute_pann(pcs, 60, 10)

  expect_equal(pann, rowMeans(nn > 60))
  expect_equal(names(pann), rownames(pcs)[1:60])
})

test_that("HNSW agrees closely with exact search", {
  set.seed(4)
  pcs <- matrix(rnorm(2000 * 10), 2000, 10)
  exact <- find_neighbors(pcs, 500, 15, method = "exact")
  approx <- find_neighbors(pcs, 500, 15, method = "hnsw")
  recall <- mean(sapply(1:500, function(i) mean(approx[i, ] %in% exact[i, ])))
  expect_gt(recall, 0.95)
})

test_that("true doublets get high pANN on two-cell-type data", {
  set.seed(5)
  genes <- 20
  type_a <- c(rep(10, 10), rep(0, 10))
  type_b <- rev(type_a)
  cell <- function(profile) pmax(profile + runif(genes, -1, 1), 0)
  counts <- cbind(
    sapply(1:300, function(i) cell(type_a)),
    sapply(1:300, function(i) cell(type_b)),
    sapply(1:30, function(i) cell((type_a + type_b) / 2))
  )
  is_doublet <- rep(c(FALSE, TRUE), c(600, 30))

  sim <- simulate_doublets(counts, n = round(630 / 0.75 - 630), seed = 1)
  merged <- t(cbind(counts, as.matrix(sim)))
  k <- round(nrow(merged) * 0.02)
  pann <- compute_pann(merged, n_real = 630, k = k)

  expect_gt(mean(pann[is_doublet]), 0.9)
  expect_lt(mean(pann[!is_doublet]), 0.5)
  top <- order(pann, decreasing = TRUE)[1:30]
  expect_gte(sum(is_doublet[top]), 27)
})

test_that("bad arguments give clear errors", {
  pcs <- matrix(rnorm(20), 10, 2)
  expect_error(find_neighbors(pcs, 11, 2), "n_real")
  expect_error(find_neighbors(pcs, 5, 0), "k must be")
  expect_error(find_neighbors(pcs, 5, 2, method = "fast"), "should be one of")
  expect_error(find_neighbors(pcs, -1, 2), "non-negative")
  pcs[1, 1] <- NA
  expect_error(find_neighbors(pcs, 5, 2), "NaN")
  expect_error(find_neighbors(as.data.frame(pcs), 5, 2), "numeric matrix")
  expect_error(simulate_doublets("x", 5), "matrix")
  expect_error(hnsw_params(max_connections = 1.5), "whole number")
  expect_error(
    find_neighbors(matrix(rnorm(20), 10, 2), 5, 2, "hnsw", hnsw_params(max_connections = 300)),
    "max_connections"
  )
})
