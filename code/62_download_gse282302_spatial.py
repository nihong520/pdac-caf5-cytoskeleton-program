"""Download GSE282302 Visium matrices + tissue positions (skip images)."""

from __future__ import annotations

import csv
import gzip
import hashlib
import time
from pathlib import Path

import requests

OUT = Path("data/raw/GSE282302_spatial")
FILELIST = OUT / "filelist.txt"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples"


def parse_needed() -> list[tuple[str, int]]:
    rows = []
    for line in FILELIST.read_text(encoding="utf-8", errors="ignore").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        name, size = parts[1], int(parts[3])
        if name.endswith("filtered_feature_bc_matrix.h5") or (
            "tissue_positions_list.csv.gz" in name
        ):
            rows.append((name, size))
    return rows


def gsm_url(name: str) -> str:
    gsm = name.split("_", 1)[0]
    prefix = gsm[:7] + "nnn"
    return f"{BASE}/{prefix}/{gsm}/suppl/{name}"


def gzip_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
        return True
    except (OSError, EOFError):
        return False


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_one(
    session: requests.Session, name: str, expected: int, destination: Path
) -> None:
    valid = (
        destination.exists()
        and destination.stat().st_size == expected
        and (destination.suffix != ".gz" or gzip_ok(destination))
    )
    if valid:
        print(f"EXISTS\t{name}", flush=True)
        return
    if destination.exists():
        print(f"REPLACE_INVALID\t{name}", flush=True)
        destination.unlink()

    url = gsm_url(name)
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, 6):
        try:
            with session.get(url, stream=True, timeout=(60, 600)) as response:
                response.raise_for_status()
                with partial.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if partial.stat().st_size != expected:
                raise IOError(
                    f"size mismatch: got {partial.stat().st_size}, expected {expected}"
                )
            if destination.suffix == ".gz" and not gzip_ok(partial):
                raise IOError("gzip integrity failed")
            partial.replace(destination)
            print(f"DOWNLOADED\t{name}\t{expected}", flush=True)
            return
        except Exception as exc:
            if partial.exists():
                partial.unlink()
            if destination.exists():
                destination.unlink()
            if attempt == 5:
                raise
            print(f"RETRY\t{name}\t{exc}", flush=True)
            time.sleep(3 * attempt)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    needed = parse_needed()
    if len(needed) != 216:
        raise SystemExit(f"Expected 216 files, found {len(needed)}")

    manifest = []
    with requests.Session() as session:
        for index, (name, size) in enumerate(needed, start=1):
            destination = OUT / name
            print(f"[{index}/{len(needed)}]", flush=True)
            download_one(session, name, size, destination)
            manifest.append(
                {
                    "file": name,
                    "size": destination.stat().st_size,
                    "sha256": sha256_file(destination),
                    "url": gsm_url(name),
                }
            )

    with (OUT / "download_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest[0].keys())
        writer.writeheader()
        writer.writerows(manifest)
    print("GSE282302_DOWNLOAD_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
