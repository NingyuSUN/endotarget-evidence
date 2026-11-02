#!/usr/bin/env python3
"""Build and publish the gene sets of docs/GENETIC_SET_TEST_PROTOCOL.md before the freeze.

Reads the pinned Rahmioglu 2023 supplement (local, hash-checked against manifests), fetches
the Open Targets Platform endometriosis associations (saved raw with URL, time and hash),
maps symbols to Ensembl IDs with the pinned HGNC complete set, and writes the published set
table and a build receipt. No test statistic is read. Standard library plus openpyxl.
"""
from pathlib import Path
import argparse, csv, hashlib, json, re, sys, urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
SUPP = 'data/genetics/41588_2023_1323_MOESM4_ESM.xlsx'; SUPP_SHA = '01cc0de4f3cb9aeec523ba312a17b931606d627e17e83581700b559741c7585e'
HGNC = 'data/reference/gene_annotation/hgnc_complete_set_2026-10-06.tsv'; HGNC_SHA = '8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d'
OT_URL = 'https://api.platform.opentargets.org/api/v4/graphql'; OT_DISEASE = 'MONDO_0005133'
OUT_SETS = 'reports/genetic_set_test/gene_sets.tsv'; OUT_OT = 'reports/genetic_set_test/opentargets_endometriosis_targets.tsv'
OUT_UNMAPPED = 'reports/genetic_set_test/unmapped_symbols.tsv'; RECEIPT = 'receipts/genetic_set_build.json'
RAW_OT = 'data/genetics/opentargets_endometriosis_associations.json'


def sha(path, normalize=False):
    b = Path(path).read_bytes()
    if normalize: b = b.replace(b'\r\n', b'\n')
    return hashlib.sha256(b).hexdigest()


def hgnc_maps(path):
    approved, prev, alias, prev_amb, alias_amb = {}, {}, {}, set(), set()
    with Path(path).open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            eg = r['ensembl_gene_id'].strip()
            if not eg or r['status'] != 'Approved': continue
            approved[r['symbol']] = eg
            for col, store, amb in (('prev_symbol', prev, prev_amb), ('alias_symbol', alias, alias_amb)):
                for s in filter(None, (x.strip() for x in r[col].split('|'))):
                    if s in store and store[s] != eg: amb.add(s)
                    store[s] = eg
    return approved, {k: v for k, v in prev.items() if k not in prev_amb}, {k: v for k, v in alias.items() if k not in alias_amb}


def map_symbol(sym, maps):
    approved, prev, alias = maps
    if sym in approved: return approved[sym], 'approved'
    if sym in prev and sym not in approved: return prev[sym], 'previous'
    if sym in alias and sym not in approved: return alias[sym], 'alias'
    return None, 'unmapped'


def clean(sym):
    return str(sym).replace('\xa0', ' ').strip()


def supplement_sets(path):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    def sheet(prefix): return wb[next(n for n in wb.sheetnames if n.startswith(prefix + '.'))]
    def column(ws, header_row, header):
        rows = ws.iter_rows(values_only=True); hdr = None
        for i, r in enumerate(rows, 1):
            if i < header_row: continue
            if i == header_row: hdr = list(r); idx = hdr.index(header); continue
            v = r[idx] if idx < len(r) else None
            if v is None: continue
            # cells may list several genes separated by commas; footnote cells contain spaces and are not symbols
            for token in clean(v).split(','):
                token = token.strip()
                if token != '.' and re.fullmatch(r'[A-Za-z][A-Za-z0-9._-]{0,24}', token): yield token
    nearest = sorted({s for s in column(sheet('Supp11'), 2, 'Nearest gene')})
    smr = set()
    for pref, head in (('Supp14', 'Gene'), ('Supp15', 'Gene'), ('Supp16', 'Gene'), ('Supp17', 'Closest Gene')):
        smr.update(column(sheet(pref), 2, head))
    return {'rahmioglu_nearest': nearest, 'rahmioglu_smr': sorted(smr)}


def fetch_opentargets(raw_path):
    rows, index, count = [], 0, None
    while True:
        q = {'query': '{ disease(efoId: "%s") { id name associatedTargets(page: {size: 500, index: %d}, enableIndirect: false) { count rows { target { id approvedSymbol } score datatypeScores { id score } } } } }' % (OT_DISEASE, index)}
        req = urllib.request.Request(OT_URL, data=json.dumps(q).encode(), headers={'Content-Type': 'application/json', 'User-Agent': 'endotarget-evidence/1.0'})
        with urllib.request.urlopen(req, timeout=120) as resp: data = json.loads(resp.read().decode())
        at = data['data']['disease']['associatedTargets']; count = at['count']; rows += at['rows']
        if len(rows) >= count or not at['rows']: break
        index += 1
    meta_q = {'query': '{ meta { apiVersion { x y z } dataVersion { year month iteration } } }'}
    req = urllib.request.Request(OT_URL, data=json.dumps(meta_q).encode(), headers={'Content-Type': 'application/json', 'User-Agent': 'endotarget-evidence/1.0'})
    with urllib.request.urlopen(req, timeout=60) as resp: meta = json.loads(resp.read().decode())['data']['meta']
    payload = {'url': OT_URL, 'disease': OT_DISEASE, 'fetched_utc': datetime.now(timezone.utc).isoformat(), 'count': count, 'meta': meta, 'rows': rows}
    Path(raw_path).parent.mkdir(parents=True, exist_ok=True)
    Path(raw_path).write_text(json.dumps(payload) + '\n', encoding='utf-8')
    return payload


def main(root):
    root = Path(root)
    if sha(root / SUPP) != SUPP_SHA: sys.exit('supplement hash differs from the manifest pin')
    if sha(root / HGNC, normalize=True) != HGNC_SHA and sha(root / HGNC) != HGNC_SHA: sys.exit('HGNC hash differs from the download receipt pin')
    maps = hgnc_maps(root / HGNC)
    sets = supplement_sets(root / SUPP)
    ot = fetch_opentargets(root / RAW_OT)
    ot_rows = []
    for r in ot['rows']:
        ga = next((d['score'] for d in r['datatypeScores'] if d['id'] == 'genetic_association'), None)
        ot_rows.append({'ensembl_gene_id': r['target']['id'], 'symbol': r['target']['approvedSymbol'], 'overall_score': r['score'], 'genetic_association_score': ga})
    sets['opentargets_genetic_any'] = sorted(r['ensembl_gene_id'] for r in ot_rows if r['genetic_association_score'] is not None and r['genetic_association_score'] > 0)
    sets['opentargets_genetic_0.5'] = sorted(r['ensembl_gene_id'] for r in ot_rows if r['genetic_association_score'] is not None and r['genetic_association_score'] >= 0.5)
    out_rows, unmapped = [], []
    for name in ('rahmioglu_nearest', 'rahmioglu_smr'):
        for sym in sets[name]:
            eg, how = map_symbol(sym, maps)
            if eg: out_rows.append({'set': name, 'ensembl_gene_id': eg, 'symbol': sym, 'mapping': how})
            else: unmapped.append({'set': name, 'symbol': sym})
    for name in ('opentargets_genetic_any', 'opentargets_genetic_0.5'):
        for eg in sets[name]:
            sym = next(r['symbol'] for r in ot_rows if r['ensembl_gene_id'] == eg)
            out_rows.append({'set': name, 'ensembl_gene_id': eg, 'symbol': sym, 'mapping': 'ensembl_from_source'})
    (root / OUT_SETS).parent.mkdir(parents=True, exist_ok=True)
    for rel, rows, fields in ((OUT_SETS, out_rows, ['set', 'ensembl_gene_id', 'symbol', 'mapping']),
                              (OUT_OT, sorted(ot_rows, key=lambda r: -(r['overall_score'] or 0)), ['ensembl_gene_id', 'symbol', 'overall_score', 'genetic_association_score']),
                              (OUT_UNMAPPED, unmapped, ['set', 'symbol'])):
        with (root / rel).open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    sizes = {name: {'source_symbols_or_ids': len(sets[name]), 'mapped_unique_ensembl': len({r['ensembl_gene_id'] for r in out_rows if r['set'] == name})} for name in sets}
    receipt = {'schema_version': 1, 'created_utc': datetime.now(timezone.utc).isoformat(), 'purpose': 'gene sets for the genetic set-level test, built before the freeze; no test statistic read',
               'supplement': {'path': SUPP, 'sha256': SUPP_SHA}, 'hgnc': {'path': HGNC, 'sha256_lf_normalized': HGNC_SHA},
               'opentargets': {'url': OT_URL, 'disease': OT_DISEASE, 'fetched_utc': ot['fetched_utc'], 'count': ot['count'], 'meta': ot['meta'], 'raw_saved': RAW_OT, 'raw_sha256': sha(root / RAW_OT)},
               'set_sizes': sizes, 'unmapped_symbols': len(unmapped),
               'outputs_sha256_lf': {rel: sha(root / rel, normalize=True) for rel in (OUT_SETS, OUT_OT, OUT_UNMAPPED)}}
    (root / RECEIPT).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'set_sizes': sizes, 'unmapped': len(unmapped), 'ot_count': ot['count'], 'ot_meta': ot['meta']}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    main(Path(a.root).resolve() if a.root else ROOT)
