#!/usr/bin/env python3
"""Writes bindings/r/inst/AUTHORS: the Rust crates bundled in the doubletrs
R package, with their authors and licenses (required by CRAN).

Run after changing Rust dependencies:  python3 scripts/r-authors.py
"""

import json
import subprocess
from pathlib import Path

repo = Path(__file__).resolve().parent.parent
crate_dir = repo / "bindings/r/src/rust"
meta = json.loads(subprocess.check_output(
    ["cargo", "metadata", "--format-version", "1", "--locked"], cwd=crate_dir))

# Only crates actually in the build graph, not the whole registry index.
in_graph = {node["id"] for node in meta["resolve"]["nodes"]}
own = {"doubletrs", "doublet_rs"}
crates = sorted(
    (p for p in meta["packages"] if p["id"] in in_graph and p["name"] not in own),
    key=lambda p: (p["name"], p["version"]),
)

lines = [
    "The doubletrs package bundles the following Rust crates in",
    "src/rust/vendor.tar.xz. Each is listed with its authors and license.",
    "",
]
for p in crates:
    authors = ", ".join(a.split(" <")[0] for a in p["authors"]) or f"The {p['name']} developers"
    lines.append(f"{p['name']} {p['version']} ({p['license']}): {authors}")
    if p.get("repository"):
        lines.append(f"  {p['repository']}")

out = repo / "bindings/r/inst/AUTHORS"
out.parent.mkdir(exist_ok=True)
out.write_text("\n".join(lines) + "\n")
print(f"wrote {out} ({len(crates)} crates)")
