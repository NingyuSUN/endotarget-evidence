#!/usr/bin/env python3
"""Stage 2 Step 1: Download ONLY the 31 locked GSE179640 tissue scRNA .h5 files
after verifying receipts/protocol_freeze.json, and record SHA-256 data receipt."""
from __future__ import annotations

import concurrent.futures
import csv
import datetime as dt
import hashlib
import json
import pathlib
import sys
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from freeze_protocol import sha256_file, verify_protocol_freeze  # noqa: E402


def download_one(item: dict, out_dir: pathlib.Path) -> dict:
    url = item["supplementary_url"].replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
    fname = item["supplementary_file"]
    expected_bytes = int(item["bytes"])
    dest = out_dir / fname
    if dest.exists() and dest.stat().st_size == expected_bytes:
        return {
            "gsm": item["gsm"],
            "donor_id": item["donor_id"],
            "site_code": item["site_code"],
            "filename": fname,
            "url": url,
            "bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
            "action": "reused_verified_size",
        }
    part = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "endotarget-v2/2.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, part.open("wb") as out:
        while True:
            chunk = resp.read(512 * 1024)
            if not chunk:
                break
            out.write(chunk)
    actual_bytes = part.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(f"Size mismatch for {fname}: expected {expected_bytes}, got {actual_bytes}")
    part.rename(dest)
    return {
        "gsm": item["gsm"],
        "donor_id": item["donor_id"],
        "site_code": item["site_code"],
        "filename": fname,
        "url": url,
        "bytes": actual_bytes,
        "sha256": sha256_file(dest),
        "action": "downloaded",
    }


def main():
    freeze = verify_protocol_freeze(ROOT)
    freeze_sha = sha256_file(ROOT / "receipts/protocol_freeze.json")

    manifest_path = ROOT / "manifests/m0_GSE179640_samples.tsv"
    with manifest_path.open("r", encoding="utf-8", newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t") if r["included_discovery_scrna"] == "True"]

    assert len(rows) == 31, f"Expected 31 discovery scRNA rows, got {len(rows)}"
    out_dir = ROOT / "data/raw/scrna/GSE179640"
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(download_one, r, out_dir) for r in rows]
        for fut in concurrent.futures.as_completed(futs):
            records.append(fut.result())

    records.sort(key=lambda x: x["gsm"])
    receipt = {
        "schema_version": 2,
        "stage": "STAGE2_DISCOVERY_DATA_DOWNLOAD",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "protocol_freeze_sha256": freeze_sha,
        "cohort": "GSE179640",
        "files_count": len(records),
        "total_bytes": sum(r["bytes"] for r in records),
        "external_cohorts_downloaded": False,
        "files": {r["filename"]: r for r in records},
    }
    rec_path = ROOT / "receipts/discovery_GSE179640_data.json"
    rec_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "OK",
        "files_count": len(records),
        "total_mb": round(receipt["total_bytes"] / 1e6, 2),
        "receipt_sha256": sha256_file(rec_path),
    }, indent=2))


if __name__ == "__main__":
    main()
