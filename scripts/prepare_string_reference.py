#!/usr/bin/env python3
"""Pin human STRING inputs and filter permitted channels; no disease input."""
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
import urllib.request

class DeadlineExceeded(RuntimeError):
    pass

def deadline_handler(signum, frame):
    raise DeadlineExceeded("bounded reference download deadline exceeded")

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def fetch(record, existing, cache_only=False):
    url = record["url"]
    expected_size = int(record["headers"]["Content-Range"].split("/")[-1])
    target = Path("data/reference/string") / url.rsplit("/", 1)[1]
    if target.exists():
        prior = existing.get(url)
        if not prior or prior["sha256"] != digest(target) or target.stat().st_size != expected_size:
            raise ValueError("existing full input lacks matching prior hash; preserve and inspect")
        return dict(prior, current_action="verified_cache_reuse")
    if cache_only:
        raise ValueError("pinned cached source missing; no network download in cached replay")
    if expected_size > 128 * 1024**2:
        raise ValueError("reference exceeds per-file size budget")
    st = os.statvfs(".")
    if st.f_bavail * st.f_frsize - expected_size < 20 * 1024**3:
        raise ValueError("20GiB free-disk reserve would be violated")
    partial = target.with_suffix(target.suffix + ".part")
    if partial.exists():
        raise ValueError("previous partial download exists; preserve for diagnosis instead of silently replacing")
    started = time.time()
    signal.alarm(900)
    req = urllib.request.Request(url, headers={"If-Match": record["headers"]["ETag"], "Accept-Encoding": "identity", "User-Agent": "EndoTarget-Evidence/0.1 reference-download"})
    total = 0
    h = hashlib.sha256()
    try:
        with urllib.request.urlopen(req, timeout=30) as response, partial.open("xb") as out:
            if response.status != 200:
                raise ValueError("expected full HTTP200 response")
            if response.headers.get("ETag") != record["headers"]["ETag"]:
                raise ValueError("source ETag changed after schema probe")
            if int(response.headers.get("Content-Length", -1)) != expected_size:
                raise ValueError("source length changed after schema probe")
            while True:
                chunk = response.read(256 * 1024)
                if not chunk:
                    break
                total += len(chunk)
                if total > expected_size:
                    raise ValueError("body larger than pinned source size")
                h.update(chunk)
                out.write(chunk)
                if total % (8 * 1024**2) == 0:
                    print(json.dumps({"event": "download_progress", "file": target.name, "bytes": total, "expected": expected_size}), flush=True)
        if total != expected_size:
            raise ValueError("incomplete source body")
        with partial.open("rb") as f:
            if hashlib.sha256(f.read(65536)).hexdigest() != record["sha256_prefix"]:
                raise ValueError("prefix differs from pinned schema probe")
        # Traversal checks gzip CRC and decompression integrity without materializing it.
        with gzip.open(partial, "rb") as f:
            for _ in iter(lambda: f.read(1024 * 1024), b""):
                pass
        partial.rename(target)
    finally:
        signal.alarm(0)
    return {"url": url, "path": str(target), "bytes": total, "sha256": h.hexdigest(), "etag": record["headers"]["ETag"], "last_modified": record["headers"]["Last-Modified"], "download_seconds": time.time() - started, "retrieved_at": datetime.datetime.now(datetime.timezone.utc).isoformat(), "gzip_crc": "PASS", "current_action": "downloaded"}

def edge_key(a, b):
    return (a, b) if a < b else (b, a)

def validate_filter_policy(cfg):
    """Reject unsupported semantics rather than silently ignoring policy keys."""
    expected = {
        "version": "12.5", "species_taxid": 9606,
        "network_kind": "functional_association_not_physical_binding_or_causal_regulation",
        "source_score_scale": 1000, "minimum_max_direct_score": 700,
        "edge_weight": "max(experiments,database)/1000",
        "undirected_duplicate_policy": "maximum_weight", "self_loop_policy": "exclude",
    }
    for key, value in expected.items():
        if cfg.get(key) != value:
            raise ValueError("unsupported frozen filter semantics: " + key)
    permitted = cfg.get("permitted_direct_channels", [])
    if len(permitted) != 2 or set(permitted) != {"experiments", "database"}:
        raise ValueError("only direct experiments/database channels permitted")
    excluded = set(cfg.get("excluded_channels", []))
    required_exclusions = {
        "experiments_transferred", "database_transferred", "textmining",
        "textmining_transferred", "coexpression", "coexpression_transferred",
        "neighborhood", "neighborhood_transferred", "fusion", "cooccurence", "homology",
    }
    if not required_exclusions <= excluded or excluded & set(permitted):
        raise ValueError("excluded/permitted evidence policy conflicts")


def filter_edges(lines, proteins, cfg):
    """Pure streaming filter: no IO, download, disease data or expression input."""
    validate_filter_policy(cfg)
    lines = iter(lines)
    header = next(lines).split()
    expected = ["experiments", "database", "experiments_transferred", "database_transferred", "textmining", "combined_score"]
    if header[:2] != ["protein1", "protein2"] or len(set(header)) != len(header) or not all(x in header for x in expected):
        raise ValueError("missing unique endpoint/direct/transferred score columns")
    channel_indices = [header.index(x) for x in cfg["permitted_direct_channels"]]
    edges = {}
    stats = {"source_rows": 0, "below_direct_score_threshold": 0, "self_loops_excluded": 0, "qualifying_directed_rows": 0}
    species_prefix = str(cfg["species_taxid"]) + "."
    for line in lines:
        fields = line.split()
        if len(fields) != len(header):
            raise ValueError("network column count mismatch")
        stats["source_rows"] += 1
        a, b = fields[:2]
        direct = [int(fields[j]) for j in channel_indices]
        if not all(0 <= score <= cfg["source_score_scale"] for score in direct):
            raise ValueError("direct channel out of expected scale")
        score = max(direct)
        if score < cfg["minimum_max_direct_score"]:
            stats["below_direct_score_threshold"] += 1
            continue
        if a not in proteins or b not in proteins or not a.startswith(species_prefix) or not b.startswith(species_prefix):
            raise ValueError("qualifying edge has unknown/nonhuman protein")
        if a == b:
            stats["self_loops_excluded"] += 1
            continue
        stats["qualifying_directed_rows"] += 1
        key = edge_key(a, b)
        edges[key] = max(score, edges.get(key, 0))
        if len(edges) > 2000000:
            raise ValueError("sparse reference exceeds two-million-edge resource cap")
    return edges, stats


def main():
    # This audited entry point replays pinned cached inputs only. Original download
    # and runtime receipts are immutable; a separate cached-replay receipt is emitted.
    signal.signal(signal.SIGALRM, deadline_handler)
    start = time.time()
    cfg = json.loads(Path("config/network_reference.json").read_text())
    validate_filter_policy(cfg)
    records = json.loads(Path(cfg["source_probe"]).read_text())["records"]
    if len(records) != 2 or any(r["version"] != cfg["version"] or r["species_taxid"] != cfg["species_taxid"] for r in records):
        raise ValueError("pinned source probe version/species differs from configuration")
    original_receipt_path = Path("receipts/string_reference_download.json")
    original_receipt_hash = digest(original_receipt_path)
    original = json.loads(original_receipt_path.read_text())
    existing = {r["url"]: r for r in original["sources"]}
    downloads = [fetch(rec, existing, cache_only=True) for rec in sorted(records, key=lambda r: ".info." not in r["url"])]
    info = next(Path(r["path"]) for r in downloads if ".info." in r["url"])
    network = next(Path(r["path"]) for r in downloads if ".links.full." in r["url"])
    with gzip.open(info, "rt") as f:
        proteins = {r["#string_protein_id"] for r in csv.DictReader(f, delimiter="\t")}
    with gzip.open(network, "rt") as f:
        edges, stats = filter_edges(f, proteins, cfg)
    out = Path("data/reference/string/direct_evidence_protein_edges.tsv.gz")
    partial = out.with_suffix(out.suffix + ".cached_replay.partial")
    with partial.open("xb") as binary, gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as gz:
        gz.write(b"protein1\tprotein2\tweight\n")
        for (a, b), score in sorted(edges.items()):
            gz.write(f"{a}\t{b}\t{score / cfg['source_score_scale']:.3f}\n".encode())
    expected_hash = "803cb3f406f33a81857c0759878485c41f018a5442b493d14f598b5b71808e09"
    replay_hash = digest(partial)
    if replay_hash != expected_hash or original["summary"]["output_sha256"] != expected_hash:
        raise ValueError("cached replay differs from pinned original edge hash; partial preserved")
    if out.exists():
        if digest(out) != expected_hash:
            raise ValueError("original output changed; original and replay preserved")
        partial.unlink()  # Only this newly generated, verified-identical replay copy.
    else:
        partial.rename(out)
    stats.update(proteins_in_info=len(proteins), undirected_edges=len(edges), nodes_in_filtered_reference=len({x for e in edges for x in e}), output_path=str(out), output_sha256=replay_hash, wall_seconds=time.time()-start, peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == "darwin" else 1024))
    if any(stats[k] != original["summary"][k] for k in ["source_rows", "below_direct_score_threshold", "self_loops_excluded", "qualifying_directed_rows", "proteins_in_info", "undirected_edges", "nodes_in_filtered_reference"]):
        raise ValueError("cached replay summary differs from original")
    if digest(original_receipt_path) != original_receipt_hash:
        raise ValueError("original download/runtime receipt changed")
    final = {"status": "CACHED_REFERENCE_REPLAY_PASS_NOT_A_GENE_UNIVERSE_OR_RANKING", "configuration": cfg, "configuration_sha256": digest(Path("config/network_reference.json")), "code_sha256": digest(Path(__file__)), "original_receipt_path": str(original_receipt_path), "original_receipt_sha256_preserved": original_receipt_hash, "sources": downloads, "summary": stats, "network_downloads_this_round": 0, "disease_input_used": False}
    Path("receipts/string_reference_cached_replay.json").write_text(json.dumps(final, indent=2)+"\n")
    print(json.dumps(stats, indent=2), flush=True)

if __name__ == "__main__":
    main()
