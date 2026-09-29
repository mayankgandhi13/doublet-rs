# doubletrs 0.1.0

* First release.
* DoubletFinder-compatible interface: `doubletFinder()`, `paramSweep()`,
  `summarizeSweep()`, `find.pK()` and `modelHomotypic()` take the same
  arguments and return the same outputs as the DoubletFinder functions of the
  same names. With the same seed, results are identical.
* Lower-level functions that do not need Seurat: `simulate_doublets()`,
  `find_neighbors()` and `compute_pann()`.
* Doublet simulation, nearest-neighbour search and pANN scoring run in Rust,
  in parallel; `doubletrs_threads()` sets the default thread count.
