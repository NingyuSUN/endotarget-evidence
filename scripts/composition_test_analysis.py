#!/usr/bin/env python3
"""Stromal-state composition and disease, per docs/COMPOSITION_TEST_PROTOCOL.md.

Reads the published composition grids and the v1 primary design matrices, verifies the
freeze (scripts/followon_freeze.py), and fits per cohort: the logit decidualized-like share
on the v1 design (primary), CLR of every author state (secondary, BH within cohort), the
other shares, and the sensitivities of the protocol. Writes tables, a summary, a results
receipt. Runs once. ``--root`` targets another project tree (tests with synthetic grids).
"""
from pathlib import Path
import argparse, csv, json, platform, sys
from datetime import datetime, timezone
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from followon_freeze import verify, sha  # noqa: E402

CONFIG = 'config/composition_test.json'
SCRIPT_REL = 'scripts/composition_test_analysis.py'


def read_tsv(path):
    with Path(path).open(newline='', encoding='utf-8') as f: return list(csv.DictReader(f, delimiter='\t'))


def ols(y, X, coef_index, min_df, rank_tol=1e-7):
    """Least squares with a t test on one coefficient; returns a dict or a NOT_ESTIMABLE status."""
    n, k = X.shape; rank = np.linalg.matrix_rank(X, tol=rank_tol)
    if rank < k: return {'status': 'NOT_ESTIMABLE', 'reason': 'rank_deficient', 'n': n, 'rank': int(rank)}
    df = n - k
    if df < min_df: return {'status': 'NOT_ESTIMABLE', 'reason': 'residual_df_below_minimum', 'n': n, 'df': df}
    beta, *_ = np.linalg.lstsq(X, y, rcond=None); resid = y - X @ beta
    s2 = float(resid @ resid) / df; cov = s2 * np.linalg.inv(X.T @ X)
    est = float(beta[coef_index]); se = float(np.sqrt(cov[coef_index, coef_index])); t = est / se
    p = float(2 * stats.t.sf(abs(t), df)); tq = stats.t.ppf(0.975, df)
    return {'status': 'ESTIMATED', 'n': n, 'df': df, 'estimate': est, 'se': se, 't': t, 'p': p, 'ci_low': est - tq * se, 'ci_high': est + tq * se}


def mde_t(se, df, alpha=0.05, power=0.8):
    tcrit = stats.t.ppf(1 - alpha / 2, df); lo, hi = 0.0, tcrit + 5
    for _ in range(80):
        mid = (lo + hi) / 2
        pw = stats.nct.sf(tcrit, df, mid) + stats.nct.cdf(-tcrit, df, mid)
        if not np.isfinite(pw): pw = stats.norm.sf(tcrit - mid)
        lo, hi = (lo, mid) if pw >= power else (mid, hi)
    return (lo + hi) / 2 * se


def bh(p):
    p = np.asarray(p, float); m = p.size; order = np.argsort(p); ranked = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]; out = np.empty(m); out[order] = np.minimum(q, 1.0); return out


def main(root):
    root = Path(root)
    cfg, freeze, freeze_sha = verify(root, CONFIG, SCRIPT_REL)
    P = cfg['parameters']; out = cfg['outputs']; cohorts = ('discovery', 'heca')
    donors = read_tsv(root / cfg['inputs']['composition_donors']['path'])
    denom = read_tsv(root / cfg['inputs']['composition_denominators']['path'])
    results = {'analysis': cfg['analysis'], 'synthetic': bool(cfg.get('synthetic')), 'cohorts': {}}
    rows_out = []
    for c in cohorts:
        design_rows = read_tsv(root / cfg['inputs'][f'design_{c}']['path'])
        ids = [r['donor_id'] for r in design_rows]; cols = [k for k in design_rows[0] if k != 'donor_id']
        X_full = np.array([[float(r[k]) for k in cols] for r in design_rows])
        ci = cols.index(P['disease_column'])
        dn = {r['donor_id']: int(r['primary_cells']) for r in denom if r['cohort'] == c}
        if set(dn) != set(ids): sys.exit(f'{c}: donors of the design and the denominators differ')
        grid = [r for r in donors if r['cohort'] == c]
        labels = sorted({r['source_subtype_label'] for r in grid})
        part = cfg['partition'][c]
        unassigned = [l for l in labels if l not in part['decidualized'] + part['non_decidualized'] + part['other']]
        if unassigned: sys.exit(f'{c}: author labels not covered by the frozen partition: {unassigned}')
        counts = {d: {l: 0 for l in labels} for d in ids}
        for r in grid:
            if r['donor_id'] in counts: counts[r['donor_id']][r['source_subtype_label']] = int(r['cells'])
        for d in ids:
            if sum(counts[d].values()) != dn[d]: sys.exit(f'{c}: state counts of {d} do not sum to the primary denominator')
        n = np.array([dn[d] for d in ids], float)
        share = {g: np.array([sum(counts[d][l] for l in part[g]) for d in ids], float) for g in ('decidualized', 'non_decidualized', 'other')}
        logit = lambda k: np.log((k + 0.5) / (n - k + 0.5))
        asin = lambda k: np.arcsin(np.sqrt((k + 0.5) / (n + 1.0)))
        disease = X_full[:, ci]
        def fit(y, X=X_full, keep=None, coef=ci):
            if keep is not None: y, X = y[keep], X[keep]
            # drop constant non-intercept columns after an exclusion (a batch with no remaining donor)
            keepcols = [j for j in range(X.shape[1]) if j == 0 or np.ptp(X[:, j]) > 0]
            r = ols(y, X[:, keepcols], keepcols.index(coef), P['minimum_residual_df'])
            if r['status'] == 'ESTIMATED': r['mde_log_odds'] = float(mde_t(r['se'], r['df'], P['alpha'], P['mde_power'])); r['odds_ratio'] = float(np.exp(r['estimate'])); r['mde_odds_ratio'] = float(np.exp(r['mde_log_odds']))
            return r
        primary = fit(logit(share['decidualized']))
        label = 'NOT_ESTIMABLE' if primary['status'] != 'ESTIMATED' else ('composition_shift_detected' if primary['p'] < P['per_cohort_alpha'] else 'no_detectable_shift')
        secondary_shares = {g: fit(logit(share[g])) for g in ('non_decidualized', 'other')}
        # CLR per author state
        K = np.array([[counts[d][l] + 0.5 for l in labels] for d in ids]); clr = np.log(K) - np.log(K).mean(axis=1, keepdims=True)
        per_state = []
        for j, l in enumerate(labels):
            r = ols(clr[:, j], X_full, ci, P['minimum_residual_df']); r['state'] = l
            r['group'] = next(g for g in ('decidualized', 'non_decidualized', 'other') if l in part[g])
            r['mean_share_endo'] = float(np.mean([counts[d][l] / dn[d] for d, dz in zip(ids, disease) if dz == 1]))
            r['mean_share_non_endo'] = float(np.mean([counts[d][l] / dn[d] for d, dz in zip(ids, disease) if dz == 0]))
            per_state.append(r)
        ps = np.array([r['p'] if r['status'] == 'ESTIMATED' else np.nan for r in per_state]); ok = ~np.isnan(ps)
        qs = np.full(ps.shape, np.nan); qs[ok] = bh(ps[ok])
        for r, q in zip(per_state, qs): r['q_bh_within_cohort'] = None if np.isnan(q) else float(q)
        reduced_cols = [0, ci, cols.index(P['age_column'])]
        sens = {'arcsine_sqrt': fit(asin(share['decidualized'])),
                'reduced_design': fit(logit(share['decidualized']), X=X_full[:, reduced_cols], coef=1),
                'depth_at_least_500': fit(logit(share['decidualized']), keep=n >= P['depth_minimum']),
                'adenomyosis_excluded': fit(logit(share['decidualized']), keep=X_full[:, cols.index(P['adenomyosis_column'])] == 0)}
        sens['depth_at_least_500']['donors_excluded'] = int((n < P['depth_minimum']).sum())
        sens['adenomyosis_excluded']['donors_excluded'] = int((X_full[:, cols.index(P['adenomyosis_column'])] == 1).sum())
        lodo = []
        for i, d in enumerate(ids):
            r = fit(logit(share['decidualized']), keep=np.arange(len(ids)) != i); lodo.append({'left_out': d, **{k: r.get(k) for k in ('status', 'estimate', 'p')}})
        est_l = [r['estimate'] for r in lodo if r['status'] == 'ESTIMATED']
        results['cohorts'][c] = {'donors': len(ids), 'disease_counts': {'ENDO': int(disease.sum()), 'non-ENDO': int((disease == 0).sum())}, 'author_states': labels,
                                 'partition_sizes': {g: len(part[g]) for g in part},
                                 'decidualized_share_mean': {'ENDO': float(np.mean(share['decidualized'][disease == 1] / n[disease == 1])), 'non-ENDO': float(np.mean(share['decidualized'][disease == 0] / n[disease == 0]))},
                                 'primary': primary, 'label': label, 'secondary_shares': secondary_shares, 'per_state_clr': per_state, 'sensitivities': sens,
                                 'leave_one_donor_out': {'fits': lodo, 'estimate_range': [min(est_l), max(est_l)] if est_l else None,
                                                         'sign_constant': bool(est_l) and (all(e > 0 for e in est_l) or all(e < 0 for e in est_l)),
                                                         'min_p': min(r['p'] for r in lodo if r['status'] == 'ESTIMATED') if est_l else None}}
        for d, i in zip(ids, range(len(ids))):
            rows_out.append({'cohort': c, 'donor_id': d, 'disease': 'ENDO' if disease[i] == 1 else 'non-ENDO', 'primary_cells': int(n[i]),
                             'decidualized_cells': int(share['decidualized'][i]), 'non_decidualized_cells': int(share['non_decidualized'][i]), 'other_cells': int(share['other'][i]),
                             'decidualized_share': share['decidualized'][i] / n[i], 'logit_decidualized': float(logit(share['decidualized'])[i])})
    both = [results['cohorts'][c]['primary'] for c in cohorts]
    pair = 'replicated_shift' if all(b['status'] == 'ESTIMATED' and b['p'] < P['replication_alpha_each'] for b in both) and np.sign(both[0]['estimate']) == np.sign(both[1]['estimate']) else 'no_replicated_shift'
    results['pair_label'] = pair
    results['software'] = {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': __import__('scipy').__version__}
    results['finished_utc'] = datetime.now(timezone.utc).isoformat()
    outdir = root / out['directory']; outdir.mkdir(parents=True, exist_ok=True)
    with (root / out['donor_shares']).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows_out)
    with (root / out['per_state']).open('w', encoding='utf-8', newline='') as f:
        fields = ['cohort', 'state', 'group', 'status', 'estimate', 'se', 't', 'p', 'q_bh_within_cohort', 'ci_low', 'ci_high', 'mean_share_endo', 'mean_share_non_endo']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore'); w.writeheader()
        for c in cohorts:
            for r in results['cohorts'][c]['per_state_clr']: w.writerow({'cohort': c, **r})
    (root / out['summary']).write_text(json.dumps(results, indent=2, default=float) + '\n', encoding='utf-8')
    receipt = {'schema_version': 1, 'phase': cfg['results_phase'], 'created_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': bool(cfg.get('synthetic')),
               'protocol': cfg['protocol'], 'freeze_receipt': out['freeze_receipt'], 'freeze_sha256': freeze_sha,
               'input_sha256': {s['path']: s['sha256'] for s in cfg['inputs'].values()}, 'script_sha256': sha(root / SCRIPT_REL),
               'output_sha256': {rel: sha(root / rel) for rel in (out['donor_shares'], out['per_state'], out['summary'])},
               'labels': {c: results['cohorts'][c]['label'] for c in cohorts} | {'pair': pair}, 'software': results['software']}
    (root / out['results_receipt']).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return results


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    r = main(Path(a.root).resolve() if a.root else ROOT)
    print(json.dumps({'pair_label': r['pair_label'], **{c: {'label': r['cohorts'][c]['label'], 'estimate': r['cohorts'][c]['primary'].get('estimate'), 'p': r['cohorts'][c]['primary'].get('p')} for c in r['cohorts']}}))
