#!/usr/bin/env python3
"""Program-level replication in public bulk eutopic endometrium, per docs/BULK_REPLICATION_PROTOCOL.md.

Verifies the freeze and the data receipt, parses the downloaded GEO matrices and platform
annotations, maps probes to Ensembl genes (Entrez, then symbol, through the pinned HGNC set),
and per series computes E1 (secretory program disease response: gate, interaction theta,
secretory-only delta, MDE), E2 (direction transfer of the frozen B0 top-20 in proliferative
samples against expression-matched random sets) and E3 (tissue-level proliferative contrast:
moderated t, BH, MDE, agreement with the pooled single-cell t). Labels fixed by the protocol.
Runs once. ``--root`` targets another project tree (tests). Python 3.10, NumPy, SciPy.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, io, json, math, platform, re, sys, zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
import numpy as np
from scipy import stats, special

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from followon_freeze import verify, sha  # noqa: E402

CONFIG = 'config/bulk_replication.json'   # a correction run passes its own configuration with --config
SCRIPT_REL = 'scripts/bulk_replication_analysis.py'


# ----- readers --------------------------------------------------------------------
def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))


def read_series_matrix(path):
    """GEO series matrix: returns (gsm list, gene_ids(probes), values array probes x samples)."""
    ids, rows, gsms, inside = [], [], None, False
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'): inside = True; continue
            if line.startswith('!series_matrix_table_end'): break
            if not inside: continue
            parts = line.rstrip('\n').split('\t')
            if gsms is None: gsms = [p.strip('"') for p in parts[1:]]; continue
            ids.append(parts[0].strip('"'))
            rows.append([float(p) if p not in ('', 'null', 'NA', 'NaN') else np.nan for p in parts[1:]])
    return gsms, ids, np.array(rows, float)


def read_annot(path):
    """GEO .annot.gz platform annotation: probe -> (entrez ids, symbols)."""
    out = {}; header = None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith(('#', '^', '!')): continue
            parts = line.rstrip('\n').split('\t')
            if header is None:
                header = parts
                if 'ID' not in header: header = None; continue
                i_id, i_gid, i_sym = header.index('ID'), header.index('Gene ID') if 'Gene ID' in header else None, header.index('Gene symbol') if 'Gene symbol' in header else None
                continue
            if len(parts) <= i_id: continue
            gid = parts[i_gid].strip() if i_gid is not None and i_gid < len(parts) else ''
            sym = parts[i_sym].strip() if i_sym is not None and i_sym < len(parts) else ''
            out[parts[i_id]] = ([g for g in re.split(r'///|;|,', gid) if g.strip()], [s.strip() for s in re.split(r'///|;|,', sym) if s.strip()])
    return out


def hgnc_maps(path):
    by_entrez, by_symbol = {}, {}
    with Path(path).open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            eg = r['ensembl_gene_id'].strip()
            if not eg or r['status'] != 'Approved': continue
            if r['entrez_id'].strip(): by_entrez[r['entrez_id'].strip()] = eg
            by_symbol[r['symbol']] = eg
    return by_entrez, by_symbol


def probes_to_genes(probe_ids, values, annot, maps):
    """One row per Ensembl gene: the probe with the highest mean expression; probes without a unique gene dropped."""
    by_entrez, by_symbol = maps; best = {}
    means = np.nanmean(values, axis=1)
    for i, pid in enumerate(probe_ids):
        if pid not in annot: continue
        entrez, syms = annot[pid]; genes = {by_entrez[g] for g in entrez if g in by_entrez}
        if not genes: genes = {by_symbol[s] for s in syms if s in by_symbol}
        if len(genes) != 1: continue
        g = next(iter(genes))
        if g not in best or means[i] > best[g][1]: best[g] = (i, means[i])
    genes = sorted(best); return genes, values[[best[g][0] for g in genes]]


def read_xlsx_sheet(path, sheet_index=0):
    """Minimal xlsx reader (standard library): rows of strings/floats from the first worksheet."""
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        shared = []
        if 'xl/sharedStrings.xml' in z.namelist():
            for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns):
                shared.append(''.join(t.text or '' for t in si.iter('{%s}t' % ns['m'])))
        wb = ET.fromstring(z.read('xl/workbook.xml')); sheets = wb.find('m:sheets', ns)
        rels = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        rid = list(sheets)[sheet_index].get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')
        target = next(r.get('Target') for r in rels if r.get('Id') == rid)
        target = target if target.startswith('xl/') else 'xl/' + target.lstrip('/')
        rows = []
        for row in ET.iterparse(io.BytesIO(z.read(target))):
            el = row[1]
            if el.tag != '{%s}row' % ns['m']: continue
            cells = []
            for c in el.findall('m:c', ns):
                v = c.find('m:v', ns); t = c.get('t')
                if v is None: cells.append(''); continue
                cells.append(shared[int(v.text)] if t == 's' else (v.text if t in ('str', 'inlineStr') else float(v.text)))
            rows.append(cells); el.clear()
        return rows


# ----- statistics -----------------------------------------------------------------
def ols(y, X, idx, min_df):
    n, k = X.shape; rank = np.linalg.matrix_rank(X)
    if rank < k: return {'status': 'NOT_ESTIMABLE', 'reason': 'rank_deficient', 'n': n}
    df = n - k
    if df < min_df: return {'status': 'NOT_ESTIMABLE', 'reason': 'residual_df_below_minimum', 'n': n, 'df': df}
    beta, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ beta; s2 = float(r @ r) / df
    cov = s2 * np.linalg.inv(X.T @ X); est = float(beta[idx]); se = float(np.sqrt(cov[idx, idx])); t = est / se
    tq = stats.t.ppf(0.975, df)
    return {'status': 'ESTIMATED', 'n': n, 'df': df, 'estimate': est, 'se': se, 't': t, 'p': float(2 * stats.t.sf(abs(t), df)), 'ci_low': est - tq * se, 'ci_high': est + tq * se}


def nct_power(tcrit, df, ncp):
    up = stats.nct.sf(tcrit, df, ncp); up = np.where(np.isfinite(up), up, stats.norm.sf(tcrit - ncp))
    lo = stats.nct.cdf(-tcrit, df, ncp); lo = np.where(np.isfinite(lo), lo, 0.0); return up + lo


def mde_exact_vec(se, df, alpha, power):
    se, df = np.asarray(se, float), np.asarray(df, float); tcrit = stats.t.ppf(1 - alpha / 2, df)
    lo = np.zeros_like(se); hi = tcrit + stats.norm.ppf(power) + 2.0
    for _ in range(20):
        short = nct_power(tcrit, df, hi) < power
        if not short.any(): break
        hi = np.where(short, hi * 2, hi)
    for _ in range(60):
        mid = (lo + hi) / 2; above = nct_power(tcrit, df, mid) >= power; hi = np.where(above, mid, hi); lo = np.where(above, lo, mid)
    return (lo + hi) / 2 * se


def squeeze_var(s2, df):
    """limma-style empirical Bayes moderation of residual variances (Smyth 2004, fitFDist by moments)."""
    s2 = np.asarray(s2, float); ok = s2 > 0; z = np.log(s2[ok]); e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean(); evar = max(np.var(e, ddof=1) - special.polygamma(1, df / 2), 1e-8)
    x = trigamma_inverse(evar)
    d0 = 2 * x; s0 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    return (d0 * s0 + df * s2) / (d0 + df), float(d0), float(s0)


def trigamma_inverse(y):
    """Solve trigamma(x) = y by Newton's method (limma::trigammaInverse): the update is x + dif."""
    if y > 1e7: return 1 / np.sqrt(y)
    if y < 1e-6: return 1 / y
    x = 0.5 + 1 / y
    for _ in range(50):
        tri = special.polygamma(1, x); dif = tri * (1 - tri / y) / special.polygamma(2, x); x = x + dif
        if abs(-dif / x) < 1e-8: break
    return float(x)


def bh(p):
    p = np.asarray(p, float); m = p.size; order = np.argsort(p); ranked = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]; out = np.empty(m); out[order] = np.minimum(q, 1.0); return out


def auc(scores, labels):
    pos, neg = scores[labels == 1], scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0: return float('nan')
    return float(stats.mannwhitneyu(pos, neg, alternative='two-sided').statistic / (len(pos) * len(neg)))


# ----- per-series analysis --------------------------------------------------------
def gene_de(M, design, idx, min_df):
    """Moderated t for every gene (rows of M) on one design; returns dict of arrays."""
    n, k = design.shape; df = n - k
    XtX_inv = np.linalg.inv(design.T @ design); H = XtX_inv @ design.T
    B = M @ H.T; R = M - B @ design.T; s2 = (R ** 2).sum(axis=1) / df
    s2_mod, d0, s0 = squeeze_var(s2, df); se = np.sqrt(s2_mod * XtX_inv[idx, idx]); coef = B[:, idx]
    t = coef / se; dft = df + d0; p = 2 * stats.t.sf(np.abs(t), dft)
    return {'coef': coef, 'se': se, 't': t, 'p': p, 'q': bh(p), 'df_total': dft, 'd0': d0, 's0': s0, 'df_residual': df}


def analyze_series(gse, spec, genes, M, samples, programs, b0, pooled, P, rng):
    """samples: list of dicts with gsm, phase, disease, subgroup, age; M genes x samples aligned to samples."""
    out = {'series': gse, 'role': spec['role'], 'n_genes': len(genes), 'n_samples': len(samples)}
    phase = np.array([s['phase'] for s in samples]); dis = np.array([1.0 if s['disease'] == 'ENDO' else 0.0 for s in samples])
    sec = np.isin(phase, ['early_secretory', 'mid_secretory', 'secretory']).astype(float); pro = (phase == 'proliferative')
    age = np.array([s['age'] for s in samples], float) if all(s.get('age') is not None for s in samples) else None
    out['groups'] = {'proliferative_ENDO': int((pro & (dis == 1)).sum()), 'proliferative_non_ENDO': int((pro & (dis == 0)).sum()),
                     'secretory_ENDO': int(((sec == 1) & (dis == 1)).sum()), 'secretory_non_ENDO': int(((sec == 1) & (dis == 0)).sum())}
    gpos = {g: i for i, g in enumerate(genes)}
    ref = pro & (dis == 0)
    if ref.sum() < 3: out['status'] = 'NOT_EVALUABLE'; out['reason'] = 'fewer than 3 proliferative non-ENDO samples'; return out
    centered = M - M[:, ref].mean(axis=1, keepdims=True)
    # E1
    out['E1'] = {}
    for name, prog in programs.items():
        idx = [gpos[g] for g in prog if g in gpos]
        if len(idx) < 50: out['E1'][name] = {'status': 'NOT_EVALUABLE', 'genes_in_family': len(idx)}; continue
        beta = np.array([prog[genes[i]] for i in idx]); X = centered[idx, :]
        score = (beta @ X) / float(beta @ beta); mss = (np.sign(beta) @ X) / len(idx)
        cols = [np.ones(len(samples)), sec, dis] + ([age - age.mean()] if age is not None else [])
        gate = ols(score, np.column_stack(cols), 1, P['minimum_residual_df']); a = auc(score, sec.astype(int))
        gate_pass = gate['status'] == 'ESTIMATED' and gate['estimate'] > 0 and gate['p'] < P['gate_p_max'] and a >= P['gate_auc_min']
        res = {'genes_in_program': len(prog), 'genes_in_family': len(idx), 'gate': {**gate, 'auc': a, 'pass': bool(gate_pass)},
               'score_means': {k: float(score[m].mean()) for k, m in (('proliferative_non_ENDO', pro & (dis == 0)), ('proliferative_ENDO', pro & (dis == 1)), ('secretory_non_ENDO', (sec == 1) & (dis == 0)), ('secretory_ENDO', (sec == 1) & (dis == 1))) if m.any()},
               'projection_vs_mean_signed_shift_correlation': float(np.corrcoef(score, mss)[0, 1]) if np.std(mss) > 0 else None}
        if not gate_pass: res['label'] = 'NOT_INTERPRETABLE'; out['E1'][name] = res; continue
        X1 = np.column_stack(cols + [sec * dis]); m1 = ols(score, X1, X1.shape[1] - 1, P['minimum_residual_df'])
        if m1['status'] == 'ESTIMATED':
            m1['mde'] = float(mde_exact_vec(np.array([m1['se']]), np.array([m1['df']]), P['alpha'], P['mde_power'])[0])
            m1['relative_to_shift'] = {'ci_low_share': m1['ci_low'] / gate['estimate'], 'ci_high_share': m1['ci_high'] / gate['estimate']}
            label = 'attenuation_detected' if m1['p'] < P['per_series_alpha'] and m1['estimate'] < 0 else 'enhancement_detected' if m1['p'] < P['per_series_alpha'] and m1['estimate'] > 0 else 'no_detectable_difference'
        else: label = 'NOT_ESTIMABLE'
        smask = sec == 1; cols2 = [np.ones(int(smask.sum())), dis[smask]] + ([age[smask] - age[smask].mean()] if age is not None else [])
        m2 = ols(score[smask], np.column_stack(cols2), 1, P['minimum_residual_df'])
        sub = (phase[smask] == 'mid_secretory').astype(float)
        m2s = ols(score[smask], np.column_stack(cols2 + [sub]), 1, P['minimum_residual_df']) if 0 < sub.sum() < smask.sum() else {'status': 'NOT_EVALUABLE', 'reason': 'no sub-phase split'}
        res.update({'M1_theta': m1, 'M2_secretory_delta': m2, 'M2_subphase_adjusted': m2s, 'label': label}); out['E1'][name] = res
    # E2 / E3 in proliferative samples
    pm = pro; Mp = M[:, pm]; dp = dis[pm]
    out['E2'], out['E3'] = {}, {}
    if dp.sum() >= 3 and (dp == 0).sum() >= 3:
        contrasts = {'all_non_ENDO': dp == dp}
        sub_other = np.array([s.get('subgroup') == 'other_pathology' for s in samples])[pm]
        if sub_other.sum() >= 3: contrasts['other_pathology_controls_only'] = (dp == 1) | sub_other
        for cname, keep in contrasts.items():
            Xp = np.column_stack([np.ones(int(keep.sum())), dp[keep]] + ([age[pm][keep] - age[pm][keep].mean()] if age is not None else []))
            de = gene_de(Mp[:, keep], Xp, 1, P['minimum_residual_df']); nsig = int((de['q'] < P['alpha']).sum())
            mde = mde_exact_vec(de['se'], np.full(len(genes), de['df_total']), P['alpha'], P['mde_power'])
            expr = Mp[:, keep].mean(axis=1); bins = np.digitize(expr, np.quantile(expr, np.linspace(0, 1, P['expression_bins'] + 1)[1:-1]))
            present = [g for g in b0 if g in gpos]; agree = int(sum(np.sign(de['coef'][gpos[g]]) == b0[g] for g in present)); evaluable = len(present)
            miss = 20 - evaluable; bounds = [agree / 20, (agree + miss) / 20]
            mem = np.array([gpos[g] for g in present], int); draws = np.empty(P['matched_draws'])
            if evaluable:
                by_bin = {b: np.where(bins == b)[0] for b in np.unique(bins)}; need = {b: int((bins[mem] == b).sum()) for b in by_bin}
                sgn = np.array([b0[g] for g in present])
                for i in range(P['matched_draws']):
                    idx = np.concatenate([rng.choice(by_bin[b], need[b], replace=False) for b in by_bin if need[b]])
                    draws[i] = (np.sign(de['coef'][idx]) == sgn).sum()
                pct = float((draws < agree).mean() * 100)
            else: pct = float('nan')
            transfers = evaluable >= 16 and agree >= P['transfer_same_direction_min'] and pct >= P['transfer_percentile_min']
            out['E2'][cname] = {'same_direction': agree, 'evaluable': evaluable, 'selected': 20, 'missing_bounds': bounds, 'matched_percentile': pct,
                                'matched_reference_mean': float(draws.mean()) if evaluable else None, 'label': 'transfers' if transfers else 'does_not_transfer',
                                'genes': [{'gene_id': g, 'discovery_sign': int(b0[g]), 'bulk_coef': float(de['coef'][gpos[g]]), 'bulk_q': float(de['q'][gpos[g]])} for g in present]}
            shared = [g for g in genes if g in pooled]; bt = np.array([de['t'][gpos[g]] for g in shared]); st = np.array([pooled[g] for g in shared])
            out['E3'][cname] = {'n_ENDO': int(dp[keep].sum()), 'n_non_ENDO': int((dp[keep] == 0).sum()), 'df_residual': de['df_residual'], 'prior_df': de['d0'],
                                'genes_BH_q_lt_0_05': nsig, 'min_q': float(de['q'].min()), 'median_mde_log2': float(np.median(mde)), 'median_mde_fold': float(2 ** np.median(mde)),
                                'ci_within_1_log2_share': float(((de['coef'] - stats.t.ppf(0.975, de['df_total']) * de['se'] >= -1) & (de['coef'] + stats.t.ppf(0.975, de['df_total']) * de['se'] <= 1)).mean()),
                                'shared_with_pooled_single_cell': len(shared), 'spearman_t_vs_pooled': float(stats.spearmanr(bt, st).statistic) if len(shared) > 10 else None,
                                'sign_agreement_with_pooled': float((np.sign(bt) == np.sign(st)).mean()) if len(shared) else None}
            out.setdefault('_de_tables', {})[cname] = (de, mde)
    else: out['E2']['status'] = 'NOT_EVALUABLE'
    return out


# ----- main ------------------------------------------------------------------------
def main(root, config=CONFIG):
    root = Path(root); cfg, freeze, freeze_sha = verify(root, config, SCRIPT_REL); P = cfg['parameters']; out = cfg['outputs']
    ddir = root / cfg['downloads']['directory']; drec = json.loads((root / cfg['downloads']['data_receipt']).read_text(encoding='utf-8'))
    for name, rec in drec['files'].items():
        if hashlib.sha256((ddir / name).read_bytes()).hexdigest() != rec['sha256']: sys.exit(f'downloaded file differs from the data receipt: {name}')
    maps = hgnc_maps(root / cfg['downloads']['hgnc']['path'])
    programs = {}
    for key, name in (('phase_program_heca', 'heca_program'), ('phase_program_discovery', 'discovery_program')):
        programs[name] = {r['gene_id']: float(r['logFC']) for r in read_tsv(root / cfg['inputs'][key]['path']) if r['fit_estimable'] == 'TRUE' and float(r['adj.P.Val']) < P['program_q_threshold']}
    b0 = {r['gene_id']: int(float(r['discovery_sign'])) for r in read_tsv(root / cfg['inputs']['dossier_pool']['path']) if r['in_B0_top20'] == 'TRUE'}
    if len(b0) != 20: sys.exit('frozen B0 top-20 must have 20 genes')
    pooled = {r['gene_id']: float(r['t']) for r in read_tsv(root / cfg['inputs']['combined_per_gene']['path'])}
    rng = np.random.default_rng(P['seed']); results = {'analysis': cfg['analysis'], 'synthetic': bool(cfg.get('synthetic')), 'series': {}}
    used_rows = []; outdir = root / out['directory']; outdir.mkdir(parents=True, exist_ok=True)
    for gse, spec in cfg['cohorts'].items():
        srows = read_tsv(root / cfg['inputs'][f'samples_{gse}']['path'])
        try:
            if spec['expression'].endswith('.xlsx'):
                res = {'series': gse, 'role': spec['role'], 'status': 'not_analyzed', 'reason': 'xlsx cohort parsing not attempted in this frozen version'}
                results['series'][gse] = res; continue
            gsms, probes, vals = read_series_matrix(ddir / spec['expression'])
            if np.nanmax(vals) > P['log2_expected_range'][1] * 2: vals = np.log2(np.clip(vals, 0, None) + 1)
            genes, M = probes_to_genes(probes, vals, read_annot(ddir / spec['annotation']), maps)
            samples = []
            for r in srows:
                ch = r['characteristics']; fields = dict(kv.split(': ', 1) for kv in ch.split(' | ') if ': ' in kv)
                if 'phase_map' in spec:
                    ph = spec['phase_map'].get(fields.get(spec['phase_field'], ''), ''); dz = spec['disease_map'].get(fields.get(spec['disease_field'], ''), '')
                    sub = 'other_pathology' if fields.get(spec.get('subgroup_field', ''), '') == spec.get('other_pathology_value') else ('no_pathology' if dz == 'non-ENDO' else '')
                else:
                    text = (ch + ' ' + r['title'] + ' ' + r['source']).lower()
                    ph = next((k for k, rx in spec['phase_regex'].items() if re.search(rx, text, re.I)), '')
                    dz = 'non-ENDO' if re.search(spec['disease_regex']['non-ENDO'], text, re.I) else ('ENDO' if re.search(spec['disease_regex']['ENDO'], text, re.I) else ''); sub = ''
                if not ph or not dz or r['gsm'] not in gsms: continue
                samples.append({'gsm': r['gsm'], 'phase': ph, 'disease': dz, 'subgroup': sub, 'age': None})
            cols = [gsms.index(s['gsm']) for s in samples]; Ms = M[:, cols]
            keep = np.isfinite(Ms).all(axis=1); genes = [g for g, k in zip(genes, keep) if k]; Ms = Ms[keep]
            res = analyze_series(gse, spec, genes, Ms, samples, programs, b0, pooled, P, rng)
            res['excluded_samples'] = len(srows) - len(samples); res['log2_range'] = [float(np.nanmin(Ms)), float(np.nanmax(Ms))]
            for s in samples: used_rows.append({'series': gse, **s})
            for name, prog in programs.items():
                if res.get('E1', {}).get(name, {}).get('status') == 'NOT_EVALUABLE': continue
                idx = [i for i, g in enumerate(genes) if g in prog]
                if idx:
                    beta = np.array([prog[genes[i]] for i in idx]); ref = np.array([s['phase'] == 'proliferative' and s['disease'] == 'non-ENDO' for s in samples])
                    sc = (beta @ (Ms[idx] - Ms[idx][:, ref].mean(axis=1, keepdims=True))) / float(beta @ beta)
                    with (root / f"{out['scores_prefix']}{gse}_{name}.tsv").open('w', encoding='utf-8', newline='') as f:
                        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gsm', 'phase', 'disease', 'subgroup', 'projection_score'])
                        for s, v in zip(samples, sc): w.writerow([s['gsm'], s['phase'], s['disease'], s['subgroup'], repr(float(v))])
            for cname, (de, mde) in res.pop('_de_tables', {}).items():
                with gzip.open(root / f"{out['per_gene_prefix']}{gse}_{cname}.tsv.gz", 'wt', encoding='utf-8', newline='') as f:
                    w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 'coef_log2', 'se_moderated', 't', 'p', 'q_bh', 'mde_alpha_0.05'])
                    for i, g in enumerate(genes): w.writerow([g, repr(float(de['coef'][i])), repr(float(de['se'][i])), repr(float(de['t'][i])), repr(float(de['p'][i])), repr(float(de['q'][i])), repr(float(mde[i]))])
            results['series'][gse] = res
        except Exception as e:
            results['series'][gse] = {'series': gse, 'role': spec['role'], 'status': 'not_analyzed', 'reason': f'{type(e).__name__}: {str(e)[:200]}'}
    # cross-series E1 label on the HECA-derived program
    prim = cfg['cohorts']; lab = {}
    for gse in prim:
        e1 = results['series'].get(gse, {}).get('E1', {}).get('heca_program', {}); lab[gse] = (e1.get('label'), e1.get('M1_theta', {}).get('estimate'), e1.get('M1_theta', {}).get('p'))
    primary = next(g for g, s in prim.items() if s['role'] == 'primary'); pl = lab.get(primary)
    replicated = bool(pl and pl[2] is not None and pl[2] < P['replication_alpha_each'] and any(l[2] is not None and l[2] < P['replication_alpha_each'] and np.sign(l[1]) == np.sign(pl[1]) for g, l in lab.items() if g != primary))
    results['E1_across_series'] = {'primary_series': primary, 'labels': {g: l[0] for g, l in lab.items()}, 'label': 'replicated' if replicated else 'not_replicated'}
    results['software'] = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': __import__('scipy').__version__}; results['finished_utc'] = datetime.now(timezone.utc).isoformat()
    with (root / out['sample_table']).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['series', 'gsm', 'phase', 'disease', 'subgroup', 'age'], delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(used_rows)
    (root / out['summary']).write_text(json.dumps(results, indent=2, default=float) + '\n', encoding='utf-8')
    outputs = [out['summary'], out['sample_table']] + sorted(str(p.relative_to(root)).replace('\\', '/') for p in outdir.glob('*') if p.name not in ('summary.json', 'samples_used.tsv'))
    receipt = {'schema_version': 1, 'phase': cfg['results_phase'], 'created_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': bool(cfg.get('synthetic')), 'protocol': cfg['protocol'],
               'freeze_receipt': out['freeze_receipt'], 'freeze_sha256': freeze_sha, 'data_receipt': cfg['downloads']['data_receipt'], 'data_receipt_sha256': sha(root / cfg['downloads']['data_receipt']),
               'input_sha256': {s['path']: s['sha256'] for s in cfg['inputs'].values()}, 'script_sha256': sha(root / SCRIPT_REL), 'output_sha256': {rel: sha(root / rel) for rel in outputs},
               'labels': {g: {'E1_heca_program': r.get('E1', {}).get('heca_program', {}).get('label'), 'E2': {c: v.get('label') for c, v in r.get('E2', {}).items() if isinstance(v, dict)}} for g, r in results['series'].items()} | {'E1_across_series': results['E1_across_series']['label']},
               'software': results['software']}
    (root / out['results_receipt']).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return results


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); ap.add_argument('--config', default=CONFIG); a = ap.parse_args()
    r = main(Path(a.root).resolve() if a.root else ROOT, a.config)
    print(json.dumps({g: {'status': s.get('status', 'analyzed'), 'groups': s.get('groups'), 'E1_heca': {k: s.get('E1', {}).get('heca_program', {}).get(k) for k in ('label',)} | {'theta': s.get('E1', {}).get('heca_program', {}).get('M1_theta', {}).get('estimate'), 'p': s.get('E1', {}).get('heca_program', {}).get('M1_theta', {}).get('p')},
                           'E2': {c: (v.get('same_direction'), v.get('evaluable'), v.get('matched_percentile'), v.get('label')) for c, v in s.get('E2', {}).items() if isinstance(v, dict)},
                           'E3_sig': {c: v.get('genes_BH_q_lt_0_05') for c, v in s.get('E3', {}).items()}} for g, s in r['series'].items()} | {'across': r['E1_across_series']['label']}, default=str))
