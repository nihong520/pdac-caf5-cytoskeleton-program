"""Download tissue-only processed matrices for GSE155698 from HCA Azul."""

from __future__ import annotations

import csv
import hashlib
import json
import time
from pathlib import Path

import requests


PROJECT_ID = "bcdf233f-9246-4c0c-9843-0514120b7e3a"
API = "https://service.azul.data.humancellatlas.org/index/files"
OUT = Path("data/raw/GSE155698_tissue")


def list_files() -> list[dict]:
    filters = json.dumps({"projectId": {"is": [PROJECT_ID]}})
    url: str | None = API
    params: dict | None = {
        "filters": filters,
        "size": 100,
        "sort": "fileName",
        "order": "asc",
    }
    records: list[dict] = []
    while url:
        response = requests.get(url, params=params, timeout=120)
        response.raise_for_status()
        payload = response.json()
        for hit in payload.get("hits", []):
            for info in hit.get("files", []):
                name = str(info.get("name", ""))
                if "_TISSUE_" not in name:
                    continue
                if not name.endswith((".mtx.gz", ".tsv.gz")):
                    continue
                records.append(
                    {
                        "name": name,
                        "size": int(info["size"]),
                        "sha256": info["sha256"],
                        "url": info["azul_url"],
                        "uuid": info["uuid"],
                        "version": info["version"],
                    }
                )
        url = payload.get("pagination", {}).get("next")
        params = None
    records.sort(key=lambda x: x["name"])
    return records


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_one(session: requests.Session, record: dict) -> None:
    destination = OUT / record["name"]
    if destination.exists() and destination.stat().st_size == record["size"]:
        if sha256_file(destination) == record["sha256"]:
            print(f"VERIFIED\t{record['name']}", flush=True)
            return
        destination.unlink()

    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, 5):
        try:
            with session.get(record["url"], stream=True, timeout=(60, 300)) as response:
                response.raise_for_status()
                with partial.open("wb") as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if partial.stat().st_size != record["size"]:
                raise IOError(
                    f"size mismatch {partial.stat().st_size} != {record['size']}"
                )
            if sha256_file(partial) != record["sha256"]:
                raise IOError("sha256 mismatch")
            partial.replace(destination)
            print(f"DOWNLOADED\t{record['name']}", flush=True)
            return
        except Exception as exc:
            if partial.exists():
                partial.unlink()
            if attempt == 4:
                raise
            print(
                f"RETRY {attempt}/4\t{record['name']}\t{exc}",
                flush=True,
            )
            time.sleep(3 * attempt)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    records = list_files()
    samples = sorted(
        {
            rec["name"]
            .replace("_barcodes.tsv.gz", "")
            .replace("_features.tsv.gz", "")
            .replace("_genes.tsv.gz", "")
            .replace("_matrix.mtx.gz", "")
            for rec in records
        }
    )
    print(
        f"MANIFEST\t{len(records)} files\t{len(samples)} tissue samples\t"
        f"{sum(x['size'] for x in records) / 1024**2:.1f} MiB",
        flush=True,
    )
    with (OUT / "HCA_download_manifest.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)

    with requests.Session() as session:
        for index, record in enumerate(records, start=1):
            print(
                f"FILE {index}/{len(records)}\t{record['name']}",
                flush=True,
            )
            download_one(session, record)
    print("DOWNLOAD_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
