#!/usr/bin/env python3
"""Download the bulk expression matrices and platform annotations after the freeze, with a data receipt.

Verifies that the freeze receipt exists and has the right phase, downloads every file listed
under ``downloads.files`` of config/bulk_replication.json into data/bulk_replication/ (skipping
files already present with the same size), hashes them and writes receipts/bulk_replication_data.json
(URL, HTTP status, bytes, SHA-256, UTC per file). No expression value is parsed here. Standard
library only; intended for the analysis host.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time, urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
CONFIG = 'config/bulk_replication.json'
UA = {'User-Agent': 'endotarget-evidence bulk replication fetch (research)'}


def download(url, dest, attempts=4):
    last = None
    for _ in range(attempts):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=600) as r, dest.open('wb') as f:
                status = r.status
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk: break
                    f.write(chunk)
            return status
        except Exception as e:
            last = e; time.sleep(5)
    raise last


def main(root):
    root = Path(root); cfg = json.loads((root / CONFIG).read_text(encoding='utf-8'))
    fz = root / cfg['outputs']['freeze_receipt']
    if not fz.exists(): sys.exit('freeze receipt missing: freeze before downloading')
    freeze = json.loads(fz.read_text(encoding='utf-8'))
    if freeze.get('phase') != cfg['freeze']['phase']: sys.exit('freeze receipt phase mismatch')
    out = root / cfg['downloads']['directory']; out.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, url in cfg['downloads']['files'].items():
        dest = out / name
        status = 'cached' if dest.exists() and dest.stat().st_size > 0 else download(url, dest)
        data = dest.read_bytes()
        records[name] = {'url': url, 'http_status': status, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(), 'utc': datetime.now(timezone.utc).isoformat()}
        print(f'{name}: {status} {len(data)} bytes {records[name]["sha256"][:16]}')
    hg = cfg['downloads']['hgnc']; hp = root / hg['path']
    hdata = hp.read_bytes().replace(b'\r\n', b'\n'); hsha = hashlib.sha256(hdata).hexdigest()
    if hsha != hg['sha256_lf']: sys.exit('HGNC file differs from the pinned hash')
    receipt = {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(), 'freeze_sha256': hashlib.sha256(fz.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
               'files': records, 'hgnc': {'path': hg['path'], 'sha256_lf': hsha}, 'note': 'downloaded after the freeze; no expression value parsed by this script'}
    (root / cfg['downloads']['data_receipt']).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'files': len(records), 'receipt': cfg['downloads']['data_receipt']}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    main(Path(a.root).resolve() if a.root else ROOT)
