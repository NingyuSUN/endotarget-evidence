#!/usr/bin/env python3
"""Download free full texts of the literature candidates into reference/papers/PMID<id>/ with receipts.

Routes, in order: Europe PMC full-text XML (open-access subset), Europe PMC PDF render for the
same PMC articles, then any free publisher or repository PDF that Europe PMC lists for the
record. Every attempt is recorded (URL, HTTP status, bytes, SHA-256, UTC) in
manifests/literature_receipts.tsv; copies stay in the git-ignored reference folder. Nothing
is extracted here. Standard library only.
"""
from pathlib import Path
import csv, hashlib, json, sys, time, urllib.error, urllib.parse, urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / 'manifests/literature_candidates.tsv'
RECEIPTS = ROOT / 'manifests/literature_receipts.tsv'
PAPERS = ROOT / 'reference/papers'
UA = {'User-Agent': 'endotarget-evidence literature fetch (research; contact via GitHub NingyuSUN/endotarget-evidence)'}
FIELDS = ['pmid', 'pmcid', 'doi', 'route', 'url', 'http_status', 'content_type', 'bytes', 'sha256', 'utc', 'saved_as', 'status']


def get(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r: return r.status, r.headers.get('Content-Type', ''), r.read()
    except urllib.error.HTTPError as e: return e.code, e.headers.get('Content-Type', '') if e.headers else '', b''
    except Exception as e: return 0, str(e)[:80], b''


def core_record(pmid):
    q = urllib.parse.quote(f'EXT_ID:{pmid} AND SRC:MED')
    status, ctype, body = get(f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&resultType=core')
    if status != 200: return None
    res = json.loads(body.decode('utf-8')).get('resultList', {}).get('result', [])
    return res[0] if res else None


def record(rows, row): rows.append(row); print(f"  {row['route']:<22} {row['http_status']:>3} {row['bytes']:>9} {row['status']}")


def main(only=None):
    cands = list(csv.DictReader(CANDIDATES.open(encoding='utf-8'), delimiter='\t'))
    rows = []
    for c in cands:
        pmid = c['pmid']
        if only and pmid not in only: continue
        folder = PAPERS / f'PMID{pmid}'; folder.mkdir(parents=True, exist_ok=True)
        print(f'PMID{pmid} {c["first_author"]} {c["year"]}')
        rec = core_record(pmid); pmcid = c['pmcid'] or (rec or {}).get('pmcid', '')
        got_text = False
        def attempt(route, url, kind):
            nonlocal got_text
            status, ctype, body = get(url); utc = datetime.now(timezone.utc).isoformat()
            ok = status == 200 and len(body) > 2000 and ((kind == 'xml' and b'<article' in body[:5000]) or (kind == 'pdf' and body[:5] == b'%PDF-'))
            saved = ''
            if ok:
                name = 'full.xml' if kind == 'xml' else ('europepmc.pdf' if route.startswith('europepmc') else 'publisher.pdf')
                (folder / name).write_bytes(body); saved = str((folder / name).relative_to(ROOT)).replace('\\', '/'); got_text = True
            record(rows, {'pmid': pmid, 'pmcid': pmcid, 'doi': c['doi'], 'route': route, 'url': url, 'http_status': status, 'content_type': ctype[:60],
                          'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest() if body else '', 'utc': utc, 'saved_as': saved,
                          'status': 'saved' if ok else ('not_available' if status in (404, 403, 0) else 'rejected')})
            time.sleep(0.5)
        if pmcid:
            attempt('europepmc_xml', f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML', 'xml')
            attempt('europepmc_pdf', f'https://europepmc.org/articles/{pmcid}?pdf=render', 'pdf')
            if not got_text:   # author manuscripts outside the open-access subset are often served by NCBI efetch
                attempt('ncbi_efetch_xml', f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid[3:]}&rettype=xml&retmode=xml', 'xml')
        if not got_text and rec:
            for ft in rec.get('fullTextUrlList', {}).get('fullTextUrl', []):
                if ft.get('documentStyle') == 'pdf' and ft.get('availability', '').lower().startswith(('open', 'free')):
                    attempt(f"free_pdf:{ft.get('site', '')[:12]}", ft['url'], 'pdf')
                    if got_text: break
        if not got_text:
            record(rows, {'pmid': pmid, 'pmcid': pmcid, 'doi': c['doi'], 'route': 'none', 'url': '', 'http_status': '', 'content_type': '', 'bytes': 0, 'sha256': '',
                          'utc': datetime.now(timezone.utc).isoformat(), 'saved_as': '', 'status': 'needs_manual_download'})
        (folder / 'europepmc_record.json').write_text(json.dumps(rec, indent=1) if rec else 'null\n', encoding='utf-8')
    existing = list(csv.DictReader(RECEIPTS.open(encoding='utf-8'), delimiter='\t')) if RECEIPTS.exists() else []
    with RECEIPTS.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(existing + rows)
    saved = {r['pmid'] for r in rows if r['status'] == 'saved'}; manual = sorted({r['pmid'] for r in rows} - saved)
    print(json.dumps({'saved_pmids': len(saved), 'needs_manual_download': manual}))


if __name__ == '__main__':
    main(set(sys.argv[1:]) or None)
