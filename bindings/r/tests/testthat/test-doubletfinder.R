skip_if_not_installed("Seurat")
skip_if_not_installed("SeuratObject")

seu <- SeuratObject::pbmc_small

# Everything doubletFinder() does, written in plain R the way DoubletFinder
# does it (arithmetic doublets, dense distance matrix, sorting), using the
# same random draws.
reference_pann <- function(seu, PCs, pN, pK, seed) {
  set.seed(seed)
  data <- doubletrs:::seurat_counts(seu)
  n_real <- ncol(data)
  n_doublets <- round(n_real / (1 - pN) - n_real)
  cells1 <- sample(colnames(data), n_doublets, replace = TRUE)
  cells2 <- sample(colnames(data), n_doublets, replace = TRUE)
  doublets <- (data[, cells1] + data[, cells2]) / 2
  colnames(doublets) <- paste0("X", seq_len(n_doublets))
  pcs <- doubletrs:::merged_pcs(seu, cbind(data, doublets), PCs, sct = FALSE)
  k <- round(nrow(pcs) * pK)
  d <- as.matrix(stats::dist(pcs))
  vapply(seq_len(n_real), function(i) {
    neighbors <- order(d[, i])[2:(k + 1)]
    sum(neighbors > n_real) / k
  }, 0)
}

test_that("doubletFinder() matches a plain-R DoubletFinder computation exactly", {
  set.seed(42)
  out <- doubletFinder(seu, PCs = 1:5, pN = 0.25, pK = 0.1, nExp = 8)
  expect_equal(unname(out$pANN_0.25_0.1_8), reference_pann(seu, 1:5, 0.25, 0.1, seed = 42))

  calls <- out$DF.classifications_0.25_0.1_8
  expect_setequal(unique(calls), c("Singlet", "Doublet"))
  expect_equal(sum(calls == "Doublet"), 8)
  expect_true(all(out$pANN_0.25_0.1_8[calls == "Doublet"] >=
                  max(out$pANN_0.25_0.1_8[calls == "Singlet"])))
})

test_that("reuse.pANN only recomputes calls", {
  set.seed(1)
  out <- doubletFinder(seu, PCs = 1:5, pN = 0.25, pK = 0.1, nExp = 8)
  out <- doubletFinder(out, PCs = 1:5, pN = 0.25, pK = 0.1, nExp = 4,
                       reuse.pANN = "pANN_0.25_0.1_8")
  expect_equal(sum(out$DF.classifications_0.25_0.1_4 == "Doublet"), 4)
  expect_false("pANN_0.25_0.1_4" %in% colnames(out@meta.data))
})

test_that("paramSweep() covers the pN-pK grid and agrees with doubletFinder()", {
  set.seed(7)
  sweep <- paramSweep(seu, PCs = 1:5)
  # 80 cells: pK below 0.13 would give zero neighbours at pN = 0.05.
  pks <- seq(0.13, 0.3, by = 0.01)
  expect_length(sweep, 6 * length(pks))
  expect_equal(names(sweep)[1:2], c("pN_0.05_pK_0.13", "pN_0.05_pK_0.14"))
  expect_true(all(vapply(sweep, nrow, 0L) == ncol(seu)))
  expect_equal(rownames(sweep[[1]]), colnames(seu))

  # The first pN draws the same doublets as doubletFinder() with that seed.
  set.seed(7)
  single <- doubletFinder(seu, PCs = 1:5, pN = 0.05, pK = 0.2, nExp = 5)
  expect_equal(sweep[["pN_0.05_pK_0.2"]]$pANN, unname(single$pANN_0.05_0.2_5))
})

test_that("summarizeSweep() and find.pK() compute BCmvn", {
  set.seed(3)
  sweep <- paramSweep(seu, PCs = 1:5)
  stats <- summarizeSweep(sweep)
  expect_equal(nrow(stats), length(sweep))
  expect_s3_class(stats$pK, "factor")
  expect_true(all(is.finite(stats$BCreal)))

  bcmvn <- find.pK(stats, plot = FALSE)
  expect_equal(nrow(bcmvn), nlevels(stats$pK))
  first <- stats$BCreal[stats$pK == bcmvn$pK[1]]
  expect_equal(bcmvn$BCmetric[1], mean(first) / stats::var(first))

  truth <- ifelse(seq_len(ncol(seu)) %% 10 == 0, "Doublet", "Singlet")
  # Made-up labels can be perfectly separable, which glm() warns about.
  gt <- suppressWarnings(summarizeSweep(sweep[1:3], GT = TRUE, GT.calls = truth))
  expect_true(all(gt$AUC >= 0 & gt$AUC <= 1))
})

test_that("bimodality coefficient matches DoubletFinder's formula", {
  x <- c(stats::rnorm(50), stats::rnorm(50, 5))
  n <- length(x)
  skew <- (1 / n) * sum((x - mean(x))^3) / (((1 / n) * sum((x - mean(x))^2))^1.5)
  skew <- skew * (sqrt(n * (n - 1))) / (n - 2)
  kurt <- (1 / n) * sum((x - mean(x))^4) / (((1 / n) * sum((x - mean(x))^2))^2) - 3
  kurt <- ((n - 1) * ((n + 1) * kurt - 3 * (n - 1)) / ((n - 2) * (n - 3))) + 3
  expected <- ((skew^2) + 1) / (kurt + ((3 * ((n - 1)^2)) / ((n - 2) * (n - 3))))
  expect_identical(doubletrs:::bimodality_coefficient(x), expected)
})

test_that("modelHomotypic() and argument checks", {
  expect_equal(modelHomotypic(c("a", "a", "b", "b")), 0.5)
  expect_equal(modelHomotypic(rep("a", 3)), 1)
  expect_error(doubletFinder(seu, PCs = 1:5, pK = 0.1, nExp = 1000), "larger than")
  expect_error(doubletFinder(seu, PCs = 1:5, pK = 0.1, nExp = 5, annotations = "x"),
               "not supported")
})
