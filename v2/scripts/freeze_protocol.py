#!/usr/bin/env python3
"""Stage 1 Freeze: Hash configs/protocol_v2.yaml, all M0 manifests, STRING graph summary,
and M0 receipt into receipts/protocol_freeze.json before any expression matrix download."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import pathlib
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_protocol_freeze(root: pathlib.Path = ROOT) -> dict:
    freeze_path = root / "receipts/protocol_freeze.json"
    if not freeze_path.exists():
        raise FileNotFoundError(f"Missing {freeze_path}")
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    for relpath, info in freeze["locked_artifacts"].items():
        p = root / relpath
        if not p.exists():
            raise FileNotFoundError(f"Locked artifact missing: {relpath}")
        actual = sha256_file(p)
        if actual != info["sha256"]:
            raise ValueError(f"SHA-256 mismatch for locked artifact {relpath}: expected {info['sha256']}, got {actual}")
    return freeze


def main():
    cfg_path = ROOT / "configs/protocol_v2.yaml"
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))

    # Verify zero expression matrices exist anywhere in data/raw/scrna or data/raw/bulk
    for sub in ("scrna", "bulk"):
        d = ROOT / "data/raw" / sub
        if d.exists() and any(d.rglob("*")):
            raise RuntimeError(f"Zero-leakage violation: expression directory {d} is non-empty before protocol freeze!")

    locked_relpaths = [
        "configs/protocol_v2.yaml",
        "manifests/m0_GSE179640_samples.tsv",
        "manifests/m0_GSE179640_donor_design.tsv",
        "manifests/m0_GSE213216_samples.tsv",
        "manifests/m0_GSE213216_locked_selection.tsv",
        "manifests/m0_GSE141549_samples.tsv",
        "manifests/m0_GSE141549_subjects.tsv",
        "manifests/human_genetics_prior_v2.tsv",
        "manifests/canonical_lineage_markers.tsv",
        "manifests/curated_niche_lr_pairs.tsv",
        "manifests/string_network_summary.json",
        "receipts/m0_inputs_and_priors.json",
    ]
    locked_artifacts = {}
    for rel in locked_relpaths:
        p = ROOT / rel
        if not p.exists():
            raise FileNotFoundError(f"Required M0 artifact missing before freeze: {rel}")
        locked_artifacts[rel] = {
            "sha256": sha256_file(p),
            "bytes": p.stat().st_size,
        }

    net_summary = json.loads((ROOT / "manifests/string_network_summary.json").read_text(encoding="utf-8"))
    net_path = ROOT / net_summary["output_relpath"]
    if sha256_file(net_path) != net_summary["output_sha256"]:
        raise ValueError("STRING interim network hash mismatch!")

    freeze = {
        "schema_version": 2,
        "phase": "STAGE1_PROTOCOL_V2_FREEZE",
        "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "protocol_title": cfg["title"],
        "random_seed": cfg["runtime_policy"]["random_seed"],
        "max_threads": cfg["runtime_policy"]["max_threads"],
        "attestations": {
            "discovery_expression_downloaded": False,
            "external_expression_downloaded": False,
            "heca_annotations_used": False,
            "parameters_locked_before_expression": True,
        },
        "locked_artifacts": locked_artifacts,
        "interim_string_network": {
            "relpath": net_summary["output_relpath"],
            "sha256": net_summary["output_sha256"],
            "bytes": net_summary["output_bytes"],
        },
    }
    out_path = ROOT / "receipts/protocol_freeze.json"
    out_path.write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")
    # Self-verify immediately
    verify_protocol_freeze(ROOT)
    print(json.dumps({"status": "FROZEN", "freeze_sha256": sha256_file(out_path), "locked_count": len(locked_artifacts)}, indent=2))


if __name__ == "__main__":
    main()
