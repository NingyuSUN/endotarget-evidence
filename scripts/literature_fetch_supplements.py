#!/usr/bin/env python3
"""Download the supplementary-file bundle of open-access PMC articles with receipts.

Usage: literature_fetch_supplements.py PMID [PMID ...]
Uses the Europe PMC supplementaryFiles endpoint (a zip per article), saves it under
reference/papers/PMID<id>/supplementary/, lists and extracts its members into that folder
(no execution; archives are data), and appends a receipt row (route europepmc_supplementary)
to manifests/literature_receipts.tsv with URL, HTTP status, bytes, SHA-256 and time. Standard
library only.
"""
from pathlib import Path
import csv, hashlib, io, json, sys, urllib.error, urllib.request, zipfile
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / 'manifests/literature_candidates.tsv'
RECEIPTS = ROOT / 'manifests/literature_receipts.tsv'
FIELDS = ['pmid', 'pmcid', 'doi', 'route', 'url', 'http_status', 'content_type', 'bytes', 'sha256', 'utc', 'saved_as', 'status']
UA = {'User-Agent': 'endotarget-evidence literature fetch (research; contact via GitHub NingyuSUN/endotarget-evidence)'}


def main(pmids):
    cands = {r['pmid']: r for r in csv.DictReader(CANDIDATES.open(encoding='utf-8'), delimiter='\t')}
    rows = []
    for pmid in pmids:
        c = cands.get(pmid, {}); pmcid = c.get('pmcid', '')
        if not pmcid: print(f'PMID{pmid}: no PMCID, skipped'); continue
        url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles'
        req = urllib.request.Request(url, headers=UA)
        try:
            with urllib.request.urlopen(req, timeout=300) as r: status, ctype, body = r.status, r.headers.get('Content-Type', ''), r.read()
        except urllib.error.HTTPError as e: status, ctype, body = e.code, '', b''
        except Exception as e: status, ctype, body = 0, str(e)[:60], b''
        folder = ROOT / 'reference/papers' / f'PMID{pmid}' / 'supplementary'; folder.mkdir(parents=True, exist_ok=True)
        saved = ''; members = []
        ok = status == 200 and body[:2] == b'PK'
        if ok:
            (folder / 'supplementary.zip').write_bytes(body); saved = f'reference/papers/PMID{pmid}/supplementary/supplementary.zip'
            with zipfile.ZipFile(io.BytesIO(body)) as z:
                for info in z.infolist():
                    name = Path(info.filename).name
                    if not name or info.is_dir() or '..' in info.filename: continue
                    (folder / name).write_bytes(z.read(info)); members.append({'name': name, 'bytes': info.file_size})
            (folder / 'members.json').write_text(json.dumps(members, indent=1), encoding='utf-8')
        rows.append({'pmid': pmid, 'pmcid': pmcid, 'doi': c.get('doi', ''), 'route': 'europepmc_supplementary', 'url': url, 'http_status': status,
                     'content_type': ctype[:60], 'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest() if body else '',
                     'utc': datetime.now(timezone.utc).isoformat(), 'saved_as': saved, 'status': 'saved' if ok else 'not_available'})
        print(f'PMID{pmid} {pmcid}: {status} {len(body)} bytes; members: ' + ', '.join(f"{m['name']} ({m['bytes']})" for m in members)[:400])
    existing = list(csv.DictReader(RECEIPTS.open(encoding='utf-8'), delimiter='\t')) if RECEIPTS.exists() else []
    with RECEIPTS.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(existing + rows)


if __name__ == '__main__':
    main(sys.argv[1:])
