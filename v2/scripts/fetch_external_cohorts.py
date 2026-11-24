#!/usr/bin/env python3
"""Stage 3 Step 1: Verify BOTH protocol_freeze.json AND prediction_seal.json,
then download ONLY the 29 locked GSE213216 scRNA libraries and GSE141549 expression workbook + GPL10558.annot.gz,
recording SHA-256 receipts in receipts/external_cohorts_data.json."""
from __future__ import annotations

import concurrent.futures
import csv
import datetime as dt
import json
import pathlib
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from freeze_protocol import sha256_file, verify_protocol_freeze  # noqa: E402


def verify_prediction_seal(root: pathlib.Path = ROOT) -> dict:
    seal_path = root / "receipts/prediction_seal.json"
    if not seal_path.exists():
        raise FileNotFoundError(f"Missing {seal_path}")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    freeze_sha = sha256_file(root / "receipts/protocol_freeze.json")
    if seal["protocol_freeze_sha256"] != freeze_sha:
        raise ValueError("prediction_seal.json protocol_freeze_sha256 mismatch!")
    for relpath, info in seal["sealed_artifacts"].items():
        p = root / relpath
        if not p.exists():
            raise FileNotFoundError(f"Sealed Stage 2 artifact missing: {relpath}")
        actual = sha256_file(p)
        if actual != info["sha256"]:
            raise ValueError(f"SHA-256 mismatch on sealed Stage 2 artifact {relpath}: expected {info['sha256']}, got {actual}")
    return seal


def download_url(url: str, dest: pathlib.Path, expected_bytes: int | None = None, retries: int = 4) -> dict:
    url_https = url.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and (expected_bytes is None or dest.stat().st_size == expected_bytes):
        return {
            "filename": dest.name,
            "url": url_https,
            "bytes": dest.stat().st_size,
            "sha256": sha256_file(dest),
            "action": "reused_verified",
        }
    part = dest.with_suffix(dest.suffix + ".part")
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            if part.exists():
                part.unlink()
            req = urllib.request.Request(url_https, headers={"User-Agent": "endotarget-v2/2.0"})
            with urllib.request.urlopen(req, timeout=300) as resp, part.open("wb") as out:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            sz = part.stat().st_size
            if expected_bytes is not None and sz != expected_bytes:
                raise ValueError(f"Size mismatch for {dest.name}: expected {expected_bytes}, got {sz}")
            part.rename(dest)
            return {
                "filename": dest.name,
                "url": url_https,
                "bytes": sz,
                "sha256": sha256_file(dest),
                "action": "downloaded",
            }
        except Exception as exc:
            last_err = exc
            time.sleep(2.0 * attempt)
    raise RuntimeError(f"Failed downloading {url_https} after {retries} attempts: {last_err}")


def main():
    verify_protocol_freeze(ROOT)
    verify_prediction_seal(ROOT)
    freeze_sha = sha256_file(ROOT / "receipts/protocol_freeze.json")
    seal_sha = sha256_file(ROOT / "receipts/prediction_seal.json")

    with (ROOT / "manifests/m0_GSE213216_locked_selection.tsv").open("r", encoding="utf-8") as f:
        sel_rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(sel_rows) == 29, f"Expected 29 locked GSE213216 rows, got {len(sel_rows)}"

    # Parse exact per-file expected sizes from GSE213216/filelist.txt
    fl_sizes = {}
    for line in (ROOT / "data/raw/metadata/GSE213216/filelist.txt").read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 5 and parts[0] == "File":
            fl_sizes[parts[1].strip()] = int(parts[3].strip())

    tasks = []
    sc_dir = ROOT / "data/raw/scrna/GSE213216"
    for r in sel_rows:
        urls = [u.strip() for u in r["supplementary_urls"].split(";") if u.strip()]
        for u in urls:
            fname = u.split("/")[-1]
            tasks.append(("GSE213216", r["gsm"], u, sc_dir / fname, fl_sizes.get(fname)))

    bulk_dir = ROOT / "data/raw/bulk/GSE141549"
    tasks.append((
        "GSE141549",
        "SERIES_MATRIX_XLSX",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/GSE141549_batchCorrectednormalizedArrayscombined.xlsx",
        bulk_dir / "GSE141549_batchCorrectednormalizedArrayscombined.xlsx",
        None,
    ))
    tasks.append((
        "GSE141549",
        "GPL10558_ANNOT",
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
        bulk_dir / "GPL10558.annot.gz",
        None,
    ))

    results_sc = []
    results_bulk = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        fut_map = {
            ex.submit(download_url, url, dest, exp_sz): (cohort, tag)
            for cohort, tag, url, dest, exp_sz in tasks
        }
        for fut in concurrent.futures.as_completed(fut_map):
            cohort, tag = fut_map[fut]
            rec = fut.result()
            rec["tag"] = tag
            if cohort == "GSE213216":
                results_sc.append(rec)
            else:
                results_bulk.append(rec)

    results_sc.sort(key=lambda x: x["filename"])
    results_bulk.sort(key=lambda x: x["filename"])

    receipt = {
        "schema_version": 2,
        "stage": "STAGE3_EXTERNAL_COHORTS_DOWNLOAD",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "protocol_freeze_sha256": freeze_sha,
        "prediction_seal_sha256": seal_sha,
        "GSE213216_files_count": len(results_sc),
        "GSE213216_total_bytes": sum(r["bytes"] for r in results_sc),
        "GSE141549_files_count": len(results_bulk),
        "GSE141549_total_bytes": sum(r["bytes"] for r in results_bulk),
        "GSE213216_files": {r["filename"]: r for r in results_sc},
        "GSE141549_files": {r["filename"]: r for r in results_bulk},
    }
    rec_path = ROOT / "receipts/external_cohorts_data.json"
    rec_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "OK",
        "gse213216_files": len(results_sc),
        "gse213216_gb": round(receipt["GSE213216_total_bytes"] / 1e9, 3),
        "gse141549_mb": round(receipt["GSE141549_total_bytes"] / 1e6, 2),
        "receipt_sha256": sha256_file(rec_path),
    }, indent=2))


if __name__ == "__main__":
    main()
