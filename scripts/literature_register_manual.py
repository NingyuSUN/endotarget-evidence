#!/usr/bin/env python3
"""Register manually downloaded papers: copy each PDF into reference/papers/PMID<id>/ and append a receipt.

Usage: literature_register_manual.py <download_folder>
The folder holds one subfolder per paper named PMID<id> with a main.pdf (and optionally the
downloader's own receipt JSON, which is copied alongside and whose URL fields are recorded).
The PDF is copied byte for byte, hashed, checked against SHA256SUMS.txt in the folder when
present, and a row (route manual_download) is appended to manifests/literature_receipts.tsv.
Nothing is executed from the downloaded folder. Standard library only.
"""
from pathlib import Path
import csv, hashlib, json, re, shutil, sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT / 'manifests/literature_receipts.tsv'
CANDIDATES = ROOT / 'manifests/literature_candidates.tsv'
FIELDS = ['pmid', 'pmcid', 'doi', 'route', 'url', 'http_status', 'content_type', 'bytes', 'sha256', 'utc', 'saved_as', 'status']


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(folder):
    folder = Path(folder)
    sums = {}
    if (folder / 'SHA256SUMS.txt').exists():
        for line in (folder / 'SHA256SUMS.txt').read_text(encoding='utf-8', errors='replace').splitlines():
            parts = line.strip().split(None, 1)
            if len(parts) == 2: sums[parts[1].replace('\\', '/')] = parts[0]
    cands = {r['pmid']: r for r in csv.DictReader(CANDIDATES.open(encoding='utf-8'), delimiter='\t')}
    rows = []
    for sub in sorted(folder.glob('PMID*')):
        pmid = sub.name[4:]; pdf = sub / 'main.pdf'
        if not pdf.exists(): print(f'{sub.name}: no main.pdf, skipped'); continue
        data = pdf.read_bytes()
        if data[:5] != b'%PDF-': print(f'{sub.name}: not a PDF, skipped'); continue
        digest = hashlib.sha256(data).hexdigest(); listed = sums.get(f'{sub.name}/main.pdf')
        check = 'matches SHA256SUMS' if listed == digest else ('NOT in SHA256SUMS' if listed is None else 'DIFFERS from SHA256SUMS')
        dest = ROOT / 'reference/papers' / sub.name; dest.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(pdf, dest / 'manual.pdf')
        url = ''
        for extra in ('download_receipt.json', 'citation_metadata.json'):
            if (sub / extra).exists():
                shutil.copyfile(sub / extra, dest / f'downloader_{extra}')
                if extra == 'download_receipt.json' and not url:
                    try:
                        obj = json.loads((sub / extra).read_text(encoding='utf-8', errors='replace'))
                        url = next((str(v) for k, v in obj.items() if 'url' in k.lower() and isinstance(v, str)), '')
                    except Exception: url = ''
        c = cands.get(pmid, {})
        rows.append({'pmid': pmid, 'pmcid': c.get('pmcid', ''), 'doi': c.get('doi', ''), 'route': 'manual_download', 'url': url, 'http_status': '',
                     'content_type': 'application/pdf', 'bytes': len(data), 'sha256': digest, 'utc': datetime.now(timezone.utc).isoformat(),
                     'saved_as': f'reference/papers/{sub.name}/manual.pdf', 'status': 'saved'})
        print(f'{sub.name}: {len(data):>9} bytes  {digest[:16]}  {check}  url={url[:60]}')
    existing = list(csv.DictReader(RECEIPTS.open(encoding='utf-8'), delimiter='\t')) if RECEIPTS.exists() else []
    with RECEIPTS.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(existing + rows)
    print(f'registered {len(rows)} papers')


if __name__ == '__main__':
    main(sys.argv[1])
