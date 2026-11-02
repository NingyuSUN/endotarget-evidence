#!/usr/bin/env python3
"""Metadata-only audit of public bulk eutopic-endometrium series for a program-level replication.

For each GEO series given (default: GSE51981, GSE6364, GSE141549, GSE25628, GSE7307), downloads
the SOFT family file and the supplementary file listing, parses per-sample characteristics
(tissue, cycle phase, disease / stage, hormonal treatment where stated), counts the strata a
proliferative disease contrast and a secretory x disease interaction would need, and records
series-matrix sizes. No expression value is read. Writes manifests/bulk_replication_<GSE>_samples.tsv,
reports/bulk_replication_feasibility.json and a receipt. Standard library only; intended to run on
the analysis host.
"""
from pathlib import Path
import csv, gzip, hashlib, io, json, re, sys, urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ['GSE51981', 'GSE6364', 'GSE141549', 'GSE25628']
UA = {'User-Agent': 'endotarget-evidence metadata audit (research)'}
OUT_DIR = ROOT / 'data/metadata/bulk_replication'


def get(url, timeout=300, attempts=4):
    last = None
    for i in range(attempts):
        req = urllib.request.Request(url, headers=UA)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r: return r.status, r.read()
        except Exception as e:   # GEO FTP mirrors drop long transfers now and then; retry a few times
            last = e
    raise last


def soft_url(gse):
    stem = gse[:-3] + 'nnn'
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{gse}/soft/{gse}_family.soft.gz'


def matrix_listing(gse):
    stem = gse[:-3] + 'nnn'
    status, body = get(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{gse}/matrix/')
    return re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"[^<]*</a>\s+\S+\s+\S+\s+(\S+)', body.decode('utf-8', 'replace'))


def parse_soft(data):
    samples, series = [], {}
    cur = None
    for line in gzip.decompress(data).decode('utf-8', 'replace').splitlines():
        if line.startswith('^SAMPLE'):
            cur = {'gsm': line.split('=', 1)[1].strip(), 'characteristics': []}; samples.append(cur)
        elif line.startswith('^SERIES'): cur = None
        elif cur is None and line.startswith('!Series_'):
            k, _, v = line[1:].partition(' = '); series.setdefault(k, []).append(v)
        elif cur is not None and line.startswith('!Sample_'):
            k, _, v = line[1:].partition(' = ')
            if k == 'Sample_characteristics_ch1': cur['characteristics'].append(v)
            elif k in ('Sample_title', 'Sample_source_name_ch1', 'Sample_platform_id', 'Sample_type', 'Sample_description'): cur.setdefault(k, v)
    return series, samples


def classify(sample):
    text = ' | '.join(sample['characteristics'] + [sample.get('Sample_title', ''), sample.get('Sample_source_name_ch1', ''), sample.get('Sample_description', '')]).lower()
    phase = ('menstrual' if re.search(r'menstrual|\bme\b', text) else 'proliferative' if re.search(r'proliferative|\bpe\b|\bp\b', text) else
             'early_secretory' if re.search(r'early[- ]secretory|\bese\b', text) else 'mid_secretory' if re.search(r'mid[- ]?secretory|\bmse\b', text) else
             'late_secretory' if re.search(r'late[- ]secretory|\blse\b', text) else 'secretory' if 'secretory' in text else '')
    disease = ('endometriosis' if re.search(r'endometriosis|\bendo\b|\be\b(?!\w)', text) and not re.search(r'no endometriosis|without endometriosis|non-?endometriosis|ne\.', text) else
               'control' if re.search(r'control|healthy|normal|no endometriosis|without endometriosis|non-?endometriosis|\bne\b', text) else '')
    tissue = ('eutopic_endometrium' if re.search(r'eutopic|endometrium', text) and not re.search(r'lesion|ectopic|peritone|ovar|endometrioma', text) else
              'lesion' if re.search(r'lesion|ectopic|endometrioma|peritone', text) else 'other')
    stage = (re.search(r'stage\s*([iv]+|[1-4])', text) or [None, ''])[1]
    hormone = 'stated' if re.search(r'hormon|contracept|progestin|gnrh|oral', text) else ''
    return {'phase_guess': phase, 'disease_guess': disease, 'tissue_guess': tissue, 'stage_text': stage, 'hormone_mention': hormone}


def main(series_ids):
    OUT_DIR.mkdir(parents=True, exist_ok=True); (ROOT / 'manifests').mkdir(exist_ok=True)
    summary = {'created_utc': datetime.now(timezone.utc).isoformat(), 'series': {}, 'note': 'metadata only; guesses are regex classifications of GEO characteristics and must be confirmed against the primary publication before use'}
    for gse in series_ids:
        print('==', gse)
        status, data = get(soft_url(gse)); p = OUT_DIR / f'{gse}_family.soft.gz'; p.write_bytes(data)
        series, samples = parse_soft(data)
        rows = []
        for s in samples:
            c = classify(s); rows.append({'gsm': s['gsm'], 'title': s.get('Sample_title', ''), 'source': s.get('Sample_source_name_ch1', ''), 'platform': s.get('Sample_platform_id', ''),
                                          'characteristics': ' | '.join(s['characteristics']), **c})
        with (ROOT / f'manifests/bulk_replication_{gse}_samples.tsv').open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
        strata = {}
        for r in rows:
            if r['tissue_guess'] != 'eutopic_endometrium': continue
            k = f"{r['phase_guess'] or 'unknown'}|{r['disease_guess'] or 'unknown'}"; strata[k] = strata.get(k, 0) + 1
        try: matrices = [{'file': f, 'size': sz} for f, sz in matrix_listing(gse)]
        except Exception as e: matrices = [{'error': str(e)[:100]}]
        summary['series'][gse] = {'title': (series.get('Series_title') or [''])[0], 'platforms': sorted({r['platform'] for r in rows}), 'samples': len(rows),
                                  'eutopic_strata_phase|disease': dict(sorted(strata.items())), 'soft_sha256': hashlib.sha256(data).hexdigest(), 'soft_bytes': len(data),
                                  'series_matrix_files': matrices, 'pubmed_ids': series.get('Series_pubmed_id', []), 'summary': (series.get('Series_summary') or [''])[0][:600]}
        print(json.dumps({k: summary['series'][gse][k] for k in ('samples', 'platforms', 'eutopic_strata_phase|disease', 'series_matrix_files')}))
    (ROOT / 'reports/bulk_replication_feasibility.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main(sys.argv[1:] or DEFAULT)
