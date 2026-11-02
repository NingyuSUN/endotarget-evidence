#!/usr/bin/env python3
"""Two-cohort combined bound on the proliferative-phase stromal disease contrast.

Implements docs/COMBINED_BOUND_PROTOCOL.md on the two published per-gene tables of the
MDE analysis (saved v1 fit quantities). Fixed-effect inverse-variance combination with
Welch-Satterthwaite degrees of freedom, BH over the combined family, exact non-central t
minimum detectable effects, interval-exclusion counts, Cochran Q, Stouffer and
normal-approximation sensitivities, fixed labels. Verifies the freeze receipt and the
input pins first and refuses to run otherwise. ``--root`` points at another project tree
(tests use a temporary copy with synthetic tables marked synthetic in the configuration).
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, html, json, math, platform, sys
from datetime import datetime, timezone
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
FREEZE_PHASE = 'COMBINED_BOUND_RULES_FROZEN_PRE_ANALYSIS'
RESULT_PHASE = 'COMBINED_BOUND_COMPLETED'


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def jload(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def read_table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))


# ----- statistics -----------------------------------------------------------------
def ivw(b, s):
    """Fixed-effect inverse-variance combination of columns of estimates b and SEs s (n x k)."""
    w = 1.0 / s ** 2
    beta = (w * b).sum(axis=1) / w.sum(axis=1)
    se = 1.0 / np.sqrt(w.sum(axis=1))
    q = (w * (b - beta[:, None]) ** 2).sum(axis=1)
    return beta, se, q


def satterthwaite_df(s, df):
    """Welch-Satterthwaite df of the inverse-variance weighted combination."""
    w = 1.0 / s ** 2; a = w / w.sum(axis=1, keepdims=True)
    num = (a ** 2 * s ** 2).sum(axis=1) ** 2
    den = ((a ** 2 * s ** 2) ** 2 / df).sum(axis=1)
    return num / den


def bh(p):
    """Benjamini-Hochberg adjusted p values (monotone, capped at 1)."""
    p = np.asarray(p, float); m = p.size; order = np.argsort(p); ranked = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]; out = np.empty(m); out[order] = np.minimum(q, 1.0); return out


def nct_power(tcrit, df, ncp):
    """Two-sided power of a t test at critical value tcrit under non-centrality ncp (arrays allowed).

    scipy's non-central t can return NaN far in the tails; those cases are recovered by the
    symmetry P(T > t | ncp) = P(T < -t | -ncp) and, failing that, by their limiting values.
    """
    tcrit, df, ncp = np.broadcast_arrays(np.asarray(tcrit, float), np.asarray(df, float), np.asarray(ncp, float))
    upper = stats.nct.sf(tcrit, df, ncp)
    bad = ~np.isfinite(upper)
    if bad.any(): upper = np.where(bad, stats.nct.cdf(-tcrit, df, -ncp), upper)
    bad = ~np.isfinite(upper)
    if bad.any(): upper = np.where(bad, stats.norm.sf(tcrit - ncp), upper)
    lower = stats.nct.cdf(-tcrit, df, ncp)
    lower = np.where(np.isfinite(lower), lower, 0.0)
    return upper + lower


def mde_exact_vec(se, df, alpha, power, iterations=80):
    """Smallest |effect| detected with the given power by a two-sided t test, solved by bisection on the
    non-centrality parameter for every gene at once (exact non-central t; precision far below 1e-9)."""
    se, df = np.asarray(se, float), np.asarray(df, float)
    tcrit = stats.t.ppf(1 - alpha / 2, df)
    lo = np.zeros_like(se); hi = tcrit + stats.norm.ppf(power) + 2.0
    for _ in range(20):
        short = nct_power(tcrit, df, hi) < power
        if not short.any(): break
        hi = np.where(short, hi * 2, hi)
    for _ in range(iterations):
        mid = (lo + hi) / 2
        above = nct_power(tcrit, df, mid) >= power
        hi = np.where(above, mid, hi); lo = np.where(above, lo, mid)
    return (lo + hi) / 2 * se


def mde_exact(se, df, alpha, power):
    return float(mde_exact_vec(np.array([se]), np.array([df]), alpha, power)[0])


def mde_approx(se, df, alpha, power):
    return (stats.t.ppf(1 - alpha / 2, df) + stats.t.ppf(power, df)) * se


def stouffer(b, s, df):
    """Weighted Stouffer combination: signed z from each cohort's two-sided t p value, weights 1/se."""
    p_each = 2 * stats.t.sf(np.abs(b / s), df)
    z_each = np.sign(b) * stats.norm.isf(p_each / 2)
    w = 1.0 / s
    z = (w * z_each).sum(axis=1) / np.sqrt((w ** 2).sum(axis=1))
    return z, 2 * stats.norm.sf(np.abs(z))


# ----- figure (standard library SVG) ----------------------------------------------
def ecdf_svg(series, title, note):
    """series: list of (label, color, sorted values in log2). ECDF of per-gene MDE."""
    W, H = 900, 520; left, right, top, bottom = 90, 860, 90, 440
    xmax = 3.0
    def X(v): return left + min(v, xmax) / xmax * (right - left)
    def Y(f): return bottom - f * (bottom - top)
    esc = lambda t: html.escape(str(t), quote=True)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img">',
             '<rect width="100%" height="100%" fill="white"/><g font-family="DejaVu Sans, Arial, sans-serif" fill="#243247">',
             f'<text x="35" y="40" font-size="20">{esc(title)}</text>', f'<text x="35" y="64" font-size="12" fill="#5B6B7C">{esc(note)}</text>']
    for fr in (0, 0.25, 0.5, 0.75, 1.0):
        parts += [f'<line x1="{left}" y1="{Y(fr):.1f}" x2="{right}" y2="{Y(fr):.1f}" stroke="#dce2e9"/>',
                  f'<text x="{left - 10}" y="{Y(fr) + 4:.1f}" font-size="11" text-anchor="end" fill="#5B6B7C">{int(fr * 100)}%</text>']
    for v, lab in ((0.585, '1.5×'), (1.0, '2×'), (2.0, '4×'), (3.0, '8×')):
        parts += [f'<line x1="{X(v):.1f}" y1="{top}" x2="{X(v):.1f}" y2="{bottom}" stroke="#dce2e9" stroke-dasharray="3 3"/>',
                  f'<text x="{X(v):.1f}" y="{bottom + 16}" font-size="11" text-anchor="middle" fill="#5B6B7C">{v:g} log2 ({lab})</text>']
    parts += [f'<text x="{(left + right) / 2:.1f}" y="{bottom + 36}" font-size="12" text-anchor="middle">minimum detectable |log2 fold change| at 80% power</text>',
              f'<text x="{left - 55}" y="{(top + bottom) / 2:.1f}" font-size="12" text-anchor="middle" transform="rotate(-90 {left - 55} {(top + bottom) / 2:.1f})">share of genes</text>']
    for i, (label, color, vals) in enumerate(series):
        n = len(vals); pts = ' '.join(f'{X(v):.2f},{Y((k + 1) / n):.2f}' for k, v in enumerate(vals) if v <= xmax)
        parts += [f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"/>',
                  f'<rect x="{left + 15}" y="{top + 12 + i * 20}" width="12" height="12" fill="{color}"/>',
                  f'<text x="{left + 32}" y="{top + 23 + i * 20}" font-size="12">{esc(label)}</text>']
    parts += ['</g></svg>']
    return '\n'.join(parts) + '\n'


# ----- main ------------------------------------------------------------------------
def main(root):
    root = Path(root)
    cfg = jload(root / 'config/combined_bound.json')
    out_paths = cfg['outputs']
    freeze_path = root / out_paths['freeze_receipt']
    if not freeze_path.exists(): sys.exit('freeze receipt missing: run scripts/combined_bound_freeze.py first')
    freeze = jload(freeze_path)
    if freeze.get('phase') != FREEZE_PHASE: sys.exit('freeze receipt has the wrong phase')
    if bool(freeze.get('synthetic')) != bool(cfg.get('synthetic')): sys.exit('synthetic flag differs between freeze and configuration')
    for rel in ('config/combined_bound.json', cfg['protocol'], 'scripts/combined_bound_analysis.py'):
        if freeze['frozen_sha256'].get(rel) != sha(root / rel): sys.exit(f'frozen file changed since the freeze: {rel}')
    if (root / out_paths['results_receipt']).exists(): sys.exit('results receipt exists; this analysis runs once')
    cohorts = ('discovery', 'heca'); tables = {}
    for c in cohorts:
        spec = cfg['inputs'][c]; p = root / spec['path']
        if sha(p) != spec['sha256']: sys.exit(f'input hash mismatch: {spec["path"]}')
        if freeze['pinned_inputs'][spec['path']] != spec['sha256']: sys.exit(f'freeze pin differs from configuration: {spec["path"]}')
        rows = read_table(p)
        if len(rows) != spec['expected_rows']: sys.exit(f'row count {len(rows)} differs from expected {spec["expected_rows"]}')
        tables[c] = {r[cfg['columns']['gene_id']]: r for r in rows}
        if len(tables[c]) != len(rows): sys.exit('duplicate gene_id')
    col = cfg['columns']
    est = {c: {g for g, r in t.items() if r[col['estimable']] == 'TRUE'} for c, t in tables.items()}
    family = sorted(est['discovery'] & est['heca'])
    only = {c: sorted(est[c] - est[o]) for c, o in (('discovery', 'heca'), ('heca', 'discovery'))}
    m = len(family)
    if m < 2: sys.exit('combined family too small')
    b = np.array([[float(tables[c][g][col['estimate']]) for c in cohorts] for g in family])
    s = np.array([[float(tables[c][g][col['se']]) for c in cohorts] for g in family])
    df = np.array([[float(tables[c][g][col['df']]) for c in cohorts] for g in family])
    if not (np.isfinite(b).all() and np.isfinite(s).all() and (s > 0).all() and (df > 0).all()): sys.exit('non-finite or non-positive inputs')
    P = cfg['primary']; alpha = P['alpha']
    beta, se, q_het = ivw(b, s); nu = satterthwaite_df(s, df)
    t = beta / se; p = 2 * stats.t.sf(np.abs(t), nu); q = bh(p)
    tq = stats.t.ppf(1 - alpha / 2, nu); tq_bonf = stats.t.ppf(1 - alpha / (2 * m), nu)
    ci = (beta - tq * se, beta + tq * se); ci_b = (beta - tq_bonf * se, beta + tq_bonf * se)
    q_p = stats.chi2.sf(q_het, 1); i2 = np.maximum(0.0, (q_het - 1) / np.where(q_het > 0, q_het, np.inf))
    z_st, p_st = stouffer(b, s, df); q_st = bh(p_st)
    p_norm = 2 * stats.norm.sf(np.abs(t)); q_norm = bh(p_norm)
    scenarios = {'alpha_0.05': alpha, **{f'bh_k{k}': alpha * k / m for k in P['bh_boundaries_k']}}
    mde = {name: mde_exact_vec(se, nu, a, P['mde_power']) for name, a in scenarios.items()}
    mde_ap = {name: mde_approx(se, nu, a, P['mde_power']) for name, a in scenarios.items()}
    # single-cohort MDE at the same alphas (combined-family boundaries), for the comparison
    single = {c: {name: mde_exact_vec(s[:, j], df[:, j], a, P['mde_power']) for name, a in scenarios.items()} for j, c in enumerate(cohorts)}
    thr = P['exclusion_thresholds_log2']
    def within(lo, hi, d): return int(((lo >= -d) & (hi <= d)).sum())
    exclusion = {f'within_{d}': {'per_gene': within(ci[0], ci[1], d), 'family_bonferroni': within(ci_b[0], ci_b[1], d),
                                 'single_cohort_per_gene': {c: within(b[:, j] - stats.t.ppf(1 - alpha / 2, df[:, j]) * s[:, j],
                                                                      b[:, j] + stats.t.ppf(1 - alpha / 2, df[:, j]) * s[:, j], d) for j, c in enumerate(cohorts)}}
                 for d in thr}
    same_sign = int((np.sign(b[:, 0]) == np.sign(b[:, 1])).sum())
    pearson = float(np.corrcoef(b[:, 0], b[:, 1])[0, 1]); spearman = float(stats.spearmanr(b[:, 0], b[:, 1]).statistic)
    n_sig = int((q < alpha).sum()); het_flag = q_p < P['heterogeneity_flag_p']
    label = 'no_combined_discovery' if n_sig == 0 else 'combined_discovery'
    # outputs
    outdir = root / out_paths['directory']; outdir.mkdir(parents=True, exist_ok=True)
    header = ['gene_id', 'logFC_discovery', 'SE_discovery', 'df_discovery', 'logFC_heca', 'SE_heca', 'df_heca', 'beta_ivw', 'se_ivw', 'df_satterthwaite',
              't', 'p', 'q_bh', 'ci_low', 'ci_high', 'ci_low_bonferroni', 'ci_high_bonferroni', 'cochran_Q', 'Q_p', 'I2', 'same_sign',
              'stouffer_z', 'stouffer_p', 'stouffer_q_bh', 'normal_p', 'normal_q_bh'] + [f'mde_{n}' for n in scenarios] + [f'mde_approx_{n}' for n in scenarios]
    fmt = lambda x: repr(float(x))
    with gzip.open(root / out_paths['per_gene'], 'wt', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(header)
        for i, g in enumerate(family):
            w.writerow([g, fmt(b[i, 0]), fmt(s[i, 0]), fmt(df[i, 0]), fmt(b[i, 1]), fmt(s[i, 1]), fmt(df[i, 1]), fmt(beta[i]), fmt(se[i]), fmt(nu[i]),
                        fmt(t[i]), fmt(p[i]), fmt(q[i]), fmt(ci[0][i]), fmt(ci[1][i]), fmt(ci_b[0][i]), fmt(ci_b[1][i]), fmt(q_het[i]), fmt(q_p[i]), fmt(i2[i]),
                        'TRUE' if np.sign(b[i, 0]) == np.sign(b[i, 1]) else 'FALSE', fmt(z_st[i]), fmt(p_st[i]), fmt(q_st[i]), fmt(p_norm[i]), fmt(q_norm[i])]
                       + [fmt(mde[n][i]) for n in scenarios] + [fmt(mde_ap[n][i]) for n in scenarios])
    with (root / out_paths['cohort_only']).open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 'cohort', 'logFC', 'moderated_SE', 'df_total'])
        for c in cohorts:
            for g in only[c]: r = tables[c][g]; w.writerow([g, c, r[col['estimate']], r[col['se']], r[col['df']]])
    top = np.argsort(-np.abs(t))[:20]
    with (root / out_paths['top20']).open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['rank', 'gene_id', 'beta_ivw', 'se_ivw', 't', 'p', 'q_bh', 'logFC_discovery', 'logFC_heca', 'Q_p'])
        for r_, i in enumerate(top, 1): w.writerow([r_, family[i], fmt(beta[i]), fmt(se[i]), fmt(t[i]), fmt(p[i]), fmt(q[i]), fmt(b[i, 0]), fmt(b[i, 1]), fmt(q_p[i])])
    qs = lambda v: {k: float(np.quantile(v, x)) for k, x in (('q0.05', .05), ('q0.25', .25), ('q0.5', .5), ('q0.75', .75), ('q0.95', .95))}
    summary = {
        'analysis': cfg['analysis'], 'synthetic': bool(cfg.get('synthetic')), 'label': label, 'family_size': m,
        'cohort_only': {c: len(only[c]) for c in cohorts}, 'estimable_per_cohort': {c: len(est[c]) for c in cohorts},
        'primary': {'genes_BH_q_lt_0_05': n_sig, 'min_p': float(p.min()), 'min_q': float(q.min()), 'max_abs_t': float(np.abs(t).max()),
                    'df_satterthwaite_range': [float(nu.min()), float(nu.max())], 'se_quantiles': qs(se),
                    'significant_genes': [{'gene_id': family[i], 'beta': float(beta[i]), 'ci': [float(ci[0][i]), float(ci[1][i])], 'q': float(q[i]),
                                           'logFC_discovery': float(b[i, 0]), 'logFC_heca': float(b[i, 1]), 'Q_p': float(q_p[i])} for i in np.where(q < alpha)[0]]},
        'heterogeneity': {'genes_Q_p_lt_0_05': int(het_flag.sum()), 'share': float(het_flag.mean()), 'I2_median': float(np.median(i2)),
                          'primary_significant_among_Q_p_ge_0_05': int(((q < alpha) & ~het_flag).sum())},
        'agreement': {'same_sign': same_sign, 'same_sign_share': same_sign / m, 'pearson_logFC': pearson, 'spearman_logFC': spearman},
        'scenarios': {name: {'alpha': a, 'mde_log2_quantiles': qs(mde[name]), 'mde_fold_median': float(2 ** np.median(mde[name])),
                             'share_le_1': float((mde[name] <= 1).mean()), 'share_le_0.585': float((mde[name] <= 0.585).mean()),
                             'approx_max_abs_diff': float(np.abs(mde[name] - mde_ap[name]).max()),
                             'single_cohort_mde_log2_median': {c: float(np.median(single[c][name])) for c in cohorts}} for name, a in scenarios.items()},
        'interval_exclusion': {'thresholds_log2': thr, 'counts': exclusion, 'ci_level': P['ci_level'], 'bonferroni_alpha': alpha / m},
        'sensitivities': {'stouffer_weighted_z': {'genes_BH_q_lt_0_05': int((q_st < alpha).sum()), 'min_p': float(p_st.min())},
                          'normal_approximation': {'genes_BH_q_lt_0_05': int((q_norm < alpha).sum()), 'min_p': float(p_norm.min())}},
        'software': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': __import__('scipy').__version__},
        'finished_utc': datetime.now(timezone.utc).isoformat()}
    (root / out_paths['summary']).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    series = [('combined, α 0.05', '#228833', np.sort(mde['alpha_0.05'])),
              ('discovery alone, α 0.05', '#4477AA', np.sort(single['discovery']['alpha_0.05'])),
              ('HECA alone, α 0.05', '#EE7733', np.sort(single['heca']['alpha_0.05'])),
              ('combined, BH boundary for 100 discoveries', '#AA3377', np.sort(mde['bh_k100'])),
              ('combined, BH boundary for one discovery', '#000000', np.sort(mde['bh_k1']))]
    (root / out_paths['figure']).write_text(ecdf_svg(series, 'Per-gene minimum detectable effect: two cohorts combined versus each alone',
        f'{m:,} genes in both v1 disease fits; fixed-effect inverse-variance combination, Satterthwaite df; 80% power. Rendered from the per-gene table; no model rerun.'), encoding='utf-8')
    outputs = {k: out_paths[k] for k in ('per_gene', 'cohort_only', 'top20', 'summary', 'figure')}
    receipt = {'schema_version': 1, 'phase': RESULT_PHASE, 'created_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': bool(cfg.get('synthetic')),
               'protocol': cfg['protocol'], 'freeze_receipt': out_paths['freeze_receipt'], 'freeze_sha256': sha(freeze_path),
               'input_sha256': {cfg['inputs'][c]['path']: cfg['inputs'][c]['sha256'] for c in cohorts},
               'script_sha256': sha(Path(__file__)), 'output_sha256': {rel: sha(root / rel) for rel in outputs.values()},
               'label': label, 'family_size': m, 'genes_BH_q_lt_0_05': n_sig, 'software': summary['software']}
    (root / out_paths['results_receipt']).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', help='project root (default: this repository)')
    args = ap.parse_args()
    s = main(Path(args.root).resolve() if args.root else ROOT)
    print(json.dumps({k: s[k] for k in ('label', 'family_size')} | {'genes_BH_q_lt_0_05': s['primary']['genes_BH_q_lt_0_05'],
                     'mde_fold_median_alpha_0.05': s['scenarios']['alpha_0.05']['mde_fold_median'], 'software': s['software']}))
