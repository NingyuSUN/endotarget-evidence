#!/usr/bin/env python3
"""Build and publish the gene sets of docs/PUBLISHED_SIGNATURE_PROTOCOL.md before the freeze.

Reads the locally downloaded supplementary workbooks of Poli-Neto et al. 2020 (Scientific
Reports, Datasets 2A and 2D: proliferative-phase differentially expressed genes, stage I-II and
stage III-IV endometriosis versus healthy controls) and Bunis et al. 2022 (Frontiers in
Immunology, Table 1 sheet ``DvC_PEsamples``: proliferative-phase disease versus control), takes
each gene's published direction, maps symbols to Ensembl IDs with the pinned HGNC complete set
and writes the published set table, the unmapped and conflicting symbols and a build receipt.
No single-cell statistic is read. Standard library only (minimal xlsx reader); intended for the
analysis host where ``reference/papers/`` holds the downloaded supplements.
"""
from pathlib import Path
import argparse, csv, hashlib, io, json, re, sys, zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
HGNC = 'data/reference/gene_annotation/hgnc_complete_set_2026-10-06.tsv'
HGNC_SHA = '8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d'
RECEIPTS = 'manifests/literature_receipts.tsv'
SOURCES = {
    'polineto2020_proliferative_stage_I_II': {
        'pmid': '31941945', 'file': 'reference/papers/PMID31941945/supplementary/41598_2019_57207_MOESM3_ESM.xlsx',
        'sheet_prefix': 'DEG stage I-II endometriosis - proliferative', 'layout': 'polineto',
        'description': 'Poli-Neto et al. 2020, Dataset 2A: differentially expressed genes, eutopic endometrium, stage I-II endometriosis versus healthy controls, proliferative phase (limma meta-analysis of GSE4888, GSE6364, GSE7305, GSE51981; adjusted p < 0.05 and twofold); logFC = endometriosis minus healthy'},
    'polineto2020_proliferative_stage_III_IV': {
        'pmid': '31941945', 'file': 'reference/papers/PMID31941945/supplementary/41598_2019_57207_MOESM6_ESM.xlsx',
        'sheet_prefix': 'DEG stage III-IV endometriosis - proliferative', 'layout': 'polineto',
        'description': 'Poli-Neto et al. 2020, Dataset 2D: as above for stage III-IV endometriosis'},
    'bunis2022_proliferative_disease_vs_control': {
        'pmid': '35069565', 'file': 'reference/papers/PMID35069565/supplementary/Table_1.xlsx',
        'sheet_prefix': 'DvC_PEsamples', 'layout': 'bunis',
        'description': 'Bunis et al. 2022, Supplementary Table 1 sheet DvC_PEsamples: proliferative-phase endometrium, disease versus control (GSE51981 reanalysis, limma on Brainarray Entrez probe sets, adjusted p < 0.05 and twofold); the published logFC is control over disease, so the direction is its negative'},
}
OUT_DIR = 'reports/published_signature_test'
OUT_SETS = OUT_DIR + '/signature_sets.tsv'; OUT_UNMAPPED = OUT_DIR + '/unmapped_or_conflicting.tsv'
RECEIPT = 'receipts/published_signature_build.json'
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
RID = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'


def sha(path, normalize=False):
    b = Path(path).read_bytes()
    if normalize: b = b.replace(b'\r\n', b'\n')
    return hashlib.sha256(b).hexdigest()


def read_sheet(path, sheet_prefix):
    """Rows (lists of strings) of the first worksheet whose name starts with ``sheet_prefix``."""
    with zipfile.ZipFile(path) as z:
        shared = []
        if 'xl/sharedStrings.xml' in z.namelist():
            for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', NS):
                shared.append(''.join(t.text or '' for t in si.iter('{%s}t' % NS['m'])))
        wb = ET.fromstring(z.read('xl/workbook.xml')); rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        for sh in wb.find('m:sheets', NS):
            if not sh.get('name').startswith(sheet_prefix): continue
            target = next(r.get('Target') for r in rels if r.get('Id') == sh.get(RID))
            target = target if target.startswith('xl/') else 'xl/' + target.lstrip('/')
            rows = []
            for _, el in ET.iterparse(io.BytesIO(z.read(target))):
                if el.tag != '{%s}row' % NS['m']: continue
                cells = []
                for c in el.findall('m:c', NS):
                    v = c.find('m:v', NS); t = c.get('t')
                    if v is None: cells.append(''); continue
                    cells.append(shared[int(v.text)] if t == 's' else (v.text or ''))
                rows.append(cells); el.clear()
            return sh.get('name'), rows
    sys.exit(f'sheet starting with {sheet_prefix!r} not found in {path}')


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
    if sym in prev: return prev[sym], 'previous'
    if sym in alias: return alias[sym], 'alias'
    return None, 'unmapped'


def symbol_rows(layout, rows):
    """Yield (symbol, published direction, logFC as published, adjusted p) from a sheet."""
    if layout == 'polineto':
        # title row, then a header row ('logFC', 'AveExpr', ...) shifted one column: the symbol is the unlabeled first cell
        hdr_i = next(i for i, r in enumerate(rows) if r and r[0] == 'logFC')
        for r in rows[hdr_i + 1:]:
            if len(r) < 6 or not r[0]: continue
            lfc = float(r[1]); yield r[0].strip(), (1 if lfc > 0 else -1), lfc, float(r[5])
    elif layout == 'bunis':
        # header row names seven columns; data rows carry the probe-set ID in an unlabeled first cell, so every index shifts by one
        hdr = rows[0]; i_lfc = hdr.index('logFC (Control/Disease)') + 1; i_sym = hdr.index('gene_symbol') + 1; i_q = hdr.index('adj.P.Val') + 1
        for r in rows[1:]:
            if len(r) <= i_sym or not r[i_sym] or not r[0].endswith('_at'): continue
            lfc = float(r[i_lfc]); yield r[i_sym].strip(), (-1 if lfc > 0 else 1), lfc, float(r[i_q])
    else:
        raise ValueError(layout)


def main(root):
    root = Path(root)
    if sha(root / HGNC, normalize=True) != HGNC_SHA: sys.exit('HGNC file differs from the pinned hash')
    maps = hgnc_maps(root / HGNC)
    receipt_rows = {}
    with (root / RECEIPTS).open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['route'] == 'europepmc_supplementary' and r['status'] == 'saved': receipt_rows[r['pmid']] = {k: r[k] for k in ('url', 'http_status', 'bytes', 'sha256', 'utc')}
    sets, problems, sources = {}, [], {}
    for name, src in SOURCES.items():
        path = root / src['file']
        if not path.exists(): sys.exit(f'supplement missing: {src["file"]}')
        sheet, rows = read_sheet(path, src['sheet_prefix'])
        per_symbol = {}
        n_rows = 0
        for sym, d, lfc, q in symbol_rows(src['layout'], rows):
            n_rows += 1; per_symbol.setdefault(sym, []).append((d, lfc, q))
        members = {}
        for sym, recs in per_symbol.items():
            dirs = {d for d, _, _ in recs}
            if len(dirs) > 1:
                problems.append({'set': name, 'symbol': sym, 'ensembl_gene_id': '', 'reason': 'conflicting_directions_across_probes', 'n_rows': len(recs)}); continue
            eg, how = map_symbol(sym, maps)
            if eg is None:
                problems.append({'set': name, 'symbol': sym, 'ensembl_gene_id': '', 'reason': 'unmapped_symbol', 'n_rows': len(recs)}); continue
            if eg in members and members[eg]['direction'] != recs[0][0]:
                problems.append({'set': name, 'symbol': sym, 'ensembl_gene_id': eg, 'reason': 'conflicting_directions_across_symbols_of_one_gene', 'n_rows': len(recs)})
                members[eg]['conflict'] = True; continue
            if eg in members: continue
            best = max(recs, key=lambda x: abs(x[1]))
            members[eg] = {'symbol': sym, 'direction': recs[0][0], 'published_logFC_disease_minus_control': (best[1] if src['layout'] == 'polineto' else -best[1]), 'published_adj_p': min(x[2] for x in recs), 'mapping': how, 'conflict': False}
        for eg in [g for g, m in members.items() if m['conflict']]: del members[eg]
        sets[name] = members
        sources[name] = {**src, 'sheet': sheet, 'sha256': sha(path), 'bytes': path.stat().st_size, 'rows_read': n_rows, 'symbols': len(per_symbol), 'genes_published': len(members),
                         'direction_split': {'up_in_disease': sum(m['direction'] == 1 for m in members.values()), 'down_in_disease': sum(m['direction'] == -1 for m in members.values())},
                         'europepmc_receipt': receipt_rows.get(src['pmid'])}
    a, b = sets['polineto2020_proliferative_stage_I_II'], sets['bunis2022_proliferative_disease_vs_control']
    shared = sorted(set(a) & set(b)); agree = [g for g in shared if a[g]['direction'] == b[g]['direction']]
    sets['consensus_proliferative'] = {g: {**a[g], 'published_logFC_disease_minus_control': (a[g]['published_logFC_disease_minus_control'] + b[g]['published_logFC_disease_minus_control']) / 2, 'published_adj_p': max(a[g]['published_adj_p'], b[g]['published_adj_p']), 'mapping': a[g]['mapping'] + '+' + b[g]['mapping']} for g in agree}
    sources['consensus_proliferative'] = {'description': 'genes in both the Poli-Neto stage I-II and the Bunis proliferative sets with the same published direction', 'genes_published': len(agree),
                                          'shared_genes': len(shared), 'same_direction_share_among_shared': (len(agree) / len(shared) if shared else None),
                                          'direction_split': {'up_in_disease': sum(a[g]['direction'] == 1 for g in agree), 'down_in_disease': sum(a[g]['direction'] == -1 for g in agree)}}
    outdir = root / OUT_DIR; outdir.mkdir(parents=True, exist_ok=True)
    with (root / OUT_SETS).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['set', 'ensembl_gene_id', 'symbol', 'direction', 'published_logFC_disease_minus_control', 'published_adj_p', 'mapping'], delimiter='\t', lineterminator='\n'); w.writeheader()
        for name, members in sets.items():
            for eg in sorted(members):
                m = members[eg]; w.writerow({'set': name, 'ensembl_gene_id': eg, 'symbol': m['symbol'], 'direction': m['direction'], 'published_logFC_disease_minus_control': f"{m['published_logFC_disease_minus_control']:.6g}", 'published_adj_p': f"{m['published_adj_p']:.3g}", 'mapping': m['mapping']})
    with (root / OUT_UNMAPPED).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['set', 'symbol', 'ensembl_gene_id', 'reason', 'n_rows'], delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(problems)
    receipt = {'schema_version': 1, 'phase': 'PUBLISHED_SIGNATURE_SETS_BUILT', 'created_utc': datetime.now(timezone.utc).isoformat(),
               'hgnc': {'path': HGNC, 'sha256_lf': HGNC_SHA}, 'sources': sources, 'problems': len(problems),
               'outputs': {OUT_SETS: sha(root / OUT_SETS, normalize=True), OUT_UNMAPPED: sha(root / OUT_UNMAPPED, normalize=True)},
               'note': 'built from locally downloaded supplements with receipts in manifests/literature_receipts.tsv; no single-cell statistic was read'}
    (root / RECEIPT).parent.mkdir(parents=True, exist_ok=True)
    (root / RECEIPT).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({n: {'genes': s['genes_published'], 'split': s['direction_split']} for n, s in sources.items()}))
    print(json.dumps({'problems': len(problems), 'consensus_shared': sources['consensus_proliferative']['shared_genes'], 'same_direction_share': sources['consensus_proliferative']['same_direction_share_among_shared']}))
    return receipt


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    main(Path(a.root).resolve() if a.root else ROOT)
