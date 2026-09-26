#!/usr/bin/env Rscript
# Combines per-run CSVs from compare_doubletfinder.R into summary tables.
#
#   Rscript validation/summarize.R          # from the repo root
#
# Writes validation/results/summary_exactness.csv,
# validation/results/summary_accuracy.csv, and prints Markdown tables.

dir <- "validation/results"
read_all <- function(pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  do.call(rbind, lapply(files, read.csv))
}

exact <- read_all("_seed[0-9]+_exactness\\.csv$")
acc <- read_all("_seed[0-9]+_accuracy\\.csv$")

# Exactness: did doubletrs reproduce DoubletFinder's pANN on its own PCs?
ex <- aggregate(
  cbind(runs = 1, cells_identical_pann, cells) ~ dataset,
  data = transform(exact, runs = 1), FUN = sum
)
ex$max_abs_diff <- tapply(exact$max_abs_diff, exact$dataset, max)[ex$dataset]
ex$all_calls_identical <- tapply(exact$calls_identical, exact$dataset, all)[ex$dataset]
ex$cells <- ex$cells / ex$runs
ex$cells_identical_pann <- ex$cells_identical_pann / ex$runs
write.csv(ex, file.path(dir, "summary_exactness.csv"), row.names = FALSE)

# Accuracy: mean (sd) over seeds per dataset and method.
fmt <- function(x) sprintf("%.3f (%.3f)", mean(x), if (length(x) > 1) sd(x) else 0)
methods <- c("DoubletFinder", "doubletrs (exact)", "doubletrs (hnsw)")
datasets <- sort(unique(acc$dataset))
auprc <- sapply(methods, function(m) {
  sapply(datasets, function(d) fmt(acc$auprc[acc$dataset == d & acc$method == m]))
})
mean_auprc <- aggregate(auprc ~ dataset + method, data = acc, FUN = mean)
write.csv(mean_auprc, file.path(dir, "summary_accuracy.csv"), row.names = FALSE)

md_table <- function(df) {
  header <- paste("|", paste(names(df), collapse = " | "), "|")
  rule <- paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|")
  rows <- apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))
  cat(header, rule, rows, sep = "\n")
  cat("\n")
}

cat("## Exactness (doubletrs on DoubletFinder's PCs)\n\n")
md_table(data.frame(
  Dataset = ex$dataset,
  Cells = ex$cells,
  Runs = ex$runs,
  `Identical pANN` = sprintf("%d / %d", ex$cells_identical_pann * ex$runs, ex$cells * ex$runs),
  `Max diff` = signif(ex$max_abs_diff, 2),
  `Calls identical` = ifelse(ex$all_calls_identical, "yes", "no"),
  check.names = FALSE
))

cat("## AUPRC, mean (sd) over seeds\n\n")
md_table(data.frame(Dataset = datasets, auprc, check.names = FALSE))

overall <- sapply(methods, function(m) mean(mean_auprc$auprc[mean_auprc$method == m]))
cat("Mean AUPRC across datasets:\n")
print(round(overall, 3))
