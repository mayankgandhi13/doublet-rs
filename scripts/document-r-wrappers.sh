#!/usr/bin/env bash
# Regenerates bindings/r/R/extendr-wrappers.R from the #[extendr] functions
# in bindings/r/src/rust. Run after changing their names or arguments; the
# R package install does not regenerate it.
set -euo pipefail
cd "$(dirname "$0")/../bindings/r/src"
cargo run --quiet --manifest-path rust/Cargo.toml --bin document
echo "updated bindings/r/R/extendr-wrappers.R"
