#!/usr/bin/env bash
# Builds a self-contained source tarball of the doubletrs R package, as
# CRAN needs it: the doublet_rs core crate is copied inside the package and
# all crates.io dependencies are vendored for an offline build.
#
#   scripts/build-r-package.sh [output-dir]     # default: dist/
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
out=$(mkdir -p "${1:-$repo/dist}" && cd "${1:-$repo/dist}" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
pkg="$work/doubletrs"

# 1. Copy the R package without build artifacts.
rsync -a \
  --exclude 'src/rust/target/' --exclude 'src/.cargo/' --exclude 'src/vendor/' \
  --exclude '*.o' --exclude '*.so' --exclude '*.dll' --exclude 'src/Makevars' \
  "$repo/bindings/r/" "$pkg/"

# 2. Bundle the core crate and point the binding crate at the copy.
mkdir -p "$pkg/src/rust/doublet_rs"
cp -R "$repo/Cargo.toml" "$repo/src" "$pkg/src/rust/doublet_rs/"
sed -i.bak 's|path = "../../../.."|path = "doublet_rs"|' "$pkg/src/rust/Cargo.toml"
rm "$pkg/src/rust/Cargo.toml.bak"
grep -q 'path = "doublet_rs"' "$pkg/src/rust/Cargo.toml"

# 3. Vendor crates.io dependencies. src/Makevars unpacks vendor.tar.xz and
#    builds offline with vendor-config.toml when NOT_CRAN is unset.
(
  cd "$pkg/src/rust"
  cargo vendor --locked --versioned-dirs vendor > vendor-config.toml
  # No extended attributes: macOS adds some that other platforms' tar warns about.
  COPYFILE_DISABLE=1 tar --no-xattrs -cJf vendor.tar.xz vendor
  rm -rf vendor
)

# 4. Build the source tarball.
(cd "$out" && COPYFILE_DISABLE=1 R CMD build ${R_BUILD_ARGS:-} "$pkg")
ls -la "$out"/doubletrs_*.tar.gz
