#!/usr/bin/env python3
"""Fetch single datasets from the Xi & Li (2021) benchmark zip on Zenodo.

The zip is 752 MB, but each dataset inside is small. This reads the zip's
central directory with HTTP range requests and downloads only the members
you ask for, in resumable chunks.

Usage:
    python3 validation/fetch_dataset.py --list
    python3 validation/fetch_dataset.py hm-6k pbmc-ch     # -> validation/data/<name>.rds
"""

import argparse
import struct
import sys
import time
import urllib.request
import zipfile
import zlib
from pathlib import Path

URL = "https://zenodo.org/api/records/4062232/files/real_datasets.zip/content"
OUT_DIR = Path(__file__).resolve().parent / "data"
CHUNK = 1 << 20  # 1 MiB per request


def get_range(start: int, end: int, retries: int = 8) -> bytes:
    """Bytes [start, end] inclusive, retrying on dropped connections."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            if len(data) == end - start + 1:
                return data
        except OSError as e:
            print(f"  retry {attempt + 1}: {e}", file=sys.stderr)
        time.sleep(3)
    raise RuntimeError(f"failed to fetch bytes {start}-{end}")


class RangeFile:
    """Minimal seekable file over HTTP range requests, enough for zipfile."""

    def __init__(self):
        req = urllib.request.Request(URL, headers={"Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            self.size = int(r.headers["Content-Range"].rsplit("/", 1)[1])
        self.pos = 0

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = [offset, self.pos + offset, self.size + offset][whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        data = get_range(self.pos, self.pos + n - 1)
        self.pos += len(data)
        return data


def extract(info: zipfile.ZipInfo, dest: Path) -> None:
    # Local header: 30 fixed bytes, then file name and extra field.
    header = get_range(info.header_offset, info.header_offset + 29)
    name_len, extra_len = struct.unpack("<HH", header[26:30])
    start = info.header_offset + 30 + name_len + extra_len
    end = start + info.compress_size

    if info.compress_type == zipfile.ZIP_DEFLATED:
        inflater = zlib.decompressobj(-15)
    elif info.compress_type == zipfile.ZIP_STORED:
        inflater = None
    else:
        raise RuntimeError(f"unsupported compression {info.compress_type}")

    tmp = dest.with_suffix(".part")
    crc = 0
    with open(tmp, "wb") as out:
        for pos in range(start, end, CHUNK):
            chunk = get_range(pos, min(pos + CHUNK, end) - 1)
            data = inflater.decompress(chunk) if inflater else chunk
            crc = zlib.crc32(data, crc)
            out.write(data)
            done = min(pos + CHUNK, end) - start
            print(f"\r  {dest.name}: {done / 1e6:.1f} / {info.compress_size / 1e6:.1f} MB",
                  end="", flush=True)
        if inflater:
            tail = inflater.flush()
            crc = zlib.crc32(tail, crc)
            out.write(tail)
    print()
    if crc != info.CRC:
        tmp.unlink()
        raise RuntimeError(f"CRC mismatch for {info.filename}")
    tmp.rename(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("names", nargs="*", help="dataset names, e.g. hm-6k")
    parser.add_argument("--list", action="store_true", help="list available datasets")
    args = parser.parse_args()

    members = {
        Path(i.filename).stem: i
        for i in zipfile.ZipFile(RangeFile()).infolist()
        if i.filename.endswith(".rds")
    }
    if args.list or not args.names:
        for name, info in sorted(members.items()):
            print(f"{name:24s} {info.compress_size / 1e6:7.1f} MB")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name in args.names:
        if name not in members:
            sys.exit(f"unknown dataset {name!r}; run with --list")
        dest = OUT_DIR / f"{name}.rds"
        if dest.exists():
            print(f"  {dest.name}: already downloaded")
            continue
        extract(members[name], dest)


if __name__ == "__main__":
    main()
