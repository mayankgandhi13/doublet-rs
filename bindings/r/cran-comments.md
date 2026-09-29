## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Rust

The package uses Rust through the 'extendr' framework, following the CRAN
policy for packages using Rust:

* `SystemRequirements: Cargo (Rust's package manager), rustc >= 1.71`. The
  configure step checks for `cargo` and `rustc` and reports their versions
  in the installation log.
* All Rust dependencies are vendored in `src/rust/vendor.tar.xz` and the
  build runs offline (`cargo build --offline -j 2`); nothing is downloaded.
* The authors and licenses of the 24 bundled crates are listed in
  `inst/AUTHORS`, and "Authors of the dependency Rust crates" are listed as
  copyright holders in `Authors@R`. All are MIT, Apache-2.0, BSD-2-Clause or
  Unicode-3.0 licensed.
* `CARGO_HOME` points to a temporary directory inside `src/` during the build
  and is removed afterwards; nothing is written to the user's home directory.
* Examples, tests and the vignette use at most 2 threads.

## Test environments

* GitHub Actions: ubuntu-latest, macos-latest and windows-latest (R release),
  `R CMD check --as-cran`
* GitHub Actions: installation with the minimum supported rustc, 1.71
* Local: macOS 15 (arm64), R 4.4.2
