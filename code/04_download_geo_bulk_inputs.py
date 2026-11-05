#!/usr/bin/env python
"""Download public GEO inputs used for GSE57495 and GSE183795 validation."""
from __future__ import annotations

import gzip
import hashlib
import shutil
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
FILES = {
    RAW / "GSE57495" / "GSE57495_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57495/"
        "matrix/GSE57495_series_matrix.txt.gz"
    ),
    RAW / "GSE183795" / "GSE183795_normalized_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183795/"
        "suppl/GSE183795_normalized_matrix.txt.gz"
    ),
    RAW / "GSE183795" / "GSE183795_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183795/"
        "soft/GSE183795_family.soft.gz"
    ),
}


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size:
        print(f"USING {destination}")
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    temporary.replace(destination)
    print(f"DOWNLOADED {destination}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    records = ["file\turl\tsize_bytes\tsha256"]
    for destination, url in FILES.items():
        download(url, destination)
        records.append(
            f"{destination.relative_to(ROOT)}\t{url}\t"
            f"{destination.stat().st_size}\t{sha256(destination)}"
        )

    compressed = RAW / "GSE57495" / "GSE57495_series_matrix.txt.gz"
    matrix = compressed.with_suffix("")
    if not matrix.exists():
        with gzip.open(compressed, "rb") as source, matrix.open("wb") as target:
            shutil.copyfileobj(source, target)
        print(f"DECOMPRESSED {matrix}")

    (RAW / "geo_bulk_download_manifest.tsv").write_text(
        "\n".join(records) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
