# doubletrs

Doublet detection for single-cell RNA sequencing with the
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) method,
reimplemented in Rust to use less memory and time.

DoubletFinder builds a dense cells x cells distance matrix, which is what
makes it run out of memory on large datasets. doubletrs finds the same nearest
neighbours without ever building that matrix, and runs the heavy steps in
parallel.

## Drop-in replacement for DoubletFinder

doubletrs has DoubletFinder's functions, with the same arguments and outputs.
Replace `library(DoubletFinder)` with `library(doubletrs)`:

```r
library(doubletrs)

# seu: a Seurat object after NormalizeData, FindVariableFeatures, ScaleData, RunPCA
set.seed(1)
sweep <- paramSweep(seu, PCs = 1:10)
bcmvn <- find.pK(summarizeSweep(sweep))
pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

nExp <- round(0.075 * ncol(seu) * (1 - modelHomotypic(seu$seurat_clusters)))
seu <- doubletFinder(seu, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp)
```

Parent cells for artificial doublets are drawn with the same R random number
calls as DoubletFinder, so after the same `set.seed()` the results are
identical.

## Validation

On the 16 real datasets of the Xi & Li (2021) doublet-detection benchmark,
doubletrs reproduces DoubletFinder's pANN scores exactly for every cell, with
the same accuracy against experimentally identified doublets, and is 1.5 to
4 times faster end to end. `paramSweep()` is about 3 times faster. Details
and scripts are in the
[project repository](https://github.com/mayankgandhi13/doublet-rs#validation).

## Installation

Installing from source needs the Rust toolchain (rustc >= 1.71), available
from <https://rustup.rs>. The package shares its Rust code with the rest of
the repository, so install from a clone:

```sh
git clone https://github.com/mayankgandhi13/doublet-rs
R CMD INSTALL doublet-rs/bindings/r
```

## Reference

McGinnis CS, Murrow LM, Gartner ZJ (2019). DoubletFinder: Doublet detection
in single-cell RNA sequencing data using artificial nearest neighbors.
*Cell Systems* 8(4), 329-337. <https://doi.org/10.1016/j.cels.2019.03.003>
