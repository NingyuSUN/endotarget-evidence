#!/usr/bin/env python3
"""Published bulk disease signatures tested in the single-cell stromal results (docs/PUBLISHED_SIGNATURE_PROTOCOL.md).

T1: unsigned enrichment of each signature's genes in the pooled two-cohort |t| (competitive Wilcoxon and
an expression-matched reference). T2: directional concordance between the published direction and the
pooled single-cell statistic (signed mean t against expression-matched random sets carrying the same
direction labels; sign agreement with an exact binomial test). T3: donor-level signature score on each
cohort's pseudobulk log-CPM, fitted on the v1 primary design, with 95% interval and 80%-power MDE.
Runs once after the freeze; refuses if the results receipt exists. Python 3.10, NumPy, SciPy.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, platform, sys
from datetime import datetime, timezone
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from followon_freeze import verify, sha  # noqa: E402

CONFIG = 'config/published_signature_test.json'
SCRIPT_REL = 'scripts/published_signature_test_analysis.py'
DISEASE = 'ENDO_minus_non_ENDO'


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))


def ols(y, X, idx, min_df):
    n, k = X.shape
    if np.linalg.matrix_rank(X) < k: return {'status': 'NOT_ESTIMABLE', 'reason': 'rank_deficient', 'n': n}
    df = n - k
    if df < min_df: return {'status': 'NOT_ESTIMABLE', 'reason': 'residual_df_below_minimum', 'n': n, 'df': df}
    beta, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ beta; s2 = float(r @ r) / df
    cov = s2 * np.linalg.inv(X.T @ X); est = float(beta[idx]); se = float(np.sqrt(cov[idx, idx])); t = est / se if se > 0 else float('nan')
    tq = stats.t.ppf(0.975, df)
    return {'status': 'ESTIMATED', 'n': n, 'df': df, 'estimate': est, 'se': se, 't': float(t), 'p': float(2 * stats.t.sf(abs(t), df)), 'ci_low': est - tq * se, 'ci_high': est + tq * se}


def nct_power(tcrit, df, ncp):
    up = stats.nct.sf(tcrit, df, ncp); up = up if np.isfinite(up) else stats.norm.sf(tcrit - ncp)
    lo = stats.nct.cdf(-tcrit, df, ncp); lo = lo if np.isfinite(lo) else 0.0
    return up + lo


def mde_exact(se, df, alpha, power):
    """Smallest |effect| detected with the given power by a two-sided t test at the given SE and df (exact non-central t)."""
    tcrit = stats.t.ppf(1 - alpha / 2, df); lo, hi = 0.0, tcrit + stats.norm.ppf(power) + 2.0
    for _ in range(20):
        if nct_power(tcrit, df, hi) >= power: break
        hi *= 2
    for _ in range(60):
        mid = (lo + hi) / 2
        if nct_power(tcrit, df, mid) >= power: hi = mid
        else: lo = mid
    return (lo + hi) / 2 * se


def auc(score, is_case):
    a, b = score[is_case], score[~is_case]
    if a.size == 0 or b.size == 0: return float('nan')
    return float(((a[:, None] > b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()) / (a.size * b.size))


def matched_draws(bins, members, draws, rng):
    """Index matrix (draws x m): random family genes with the member set's expression-bin composition.

    Columns are ordered bin by bin, so a direction vector built with ``member_order`` aligns with them."""
    by_bin = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    order = np.concatenate([members[bins[members] == b] for b in by_bin])
    need = [int((bins[members] == b).sum()) for b in by_bin]
    idx = np.empty((draws, members.size), int)
    for i in range(draws):
        idx[i] = np.concatenate([rng.choice(by_bin[b], n, replace=False) for b, n in zip(by_bin, need) if n])
    return idx, order


def load_cohort(root, cfg, name):
    c = cfg['cohorts'][name]; inp = cfg['inputs']; host = cfg['host_inputs']
    passing = {r['gene_id'] for r in read_tsv(root / inp[c['expression_filter']]['path']) if r['expression_pass'] == 'TRUE'}
    lib = {r['donor_id']: float(r['effective_library_size']) for r in read_tsv(root / inp[c['normalization']]['path'])}
    design_rows = read_tsv(root / inp[c['design']]['path'])
    donors = [r['donor_id'] for r in design_rows]; cols = [k for k in design_rows[0] if k != 'donor_id']
    X = np.array([[float(r[k]) for k in cols] for r in design_rows])
    cp = root / host[c['counts']]['path']
    if hashlib.sha256(cp.read_bytes()).hexdigest() != host[c['counts']]['sha256']: sys.exit(f'host input hash differs from the pin: {host[c["counts"]]["path"]}')
    with gzip.open(cp, 'rt', encoding='utf-8', newline='') as f:
        rd = csv.reader(f, delimiter='\t'); header = next(rd); col = {h: i for i, h in enumerate(header)}
        missing = [d for d in donors if d not in col]
        if missing: sys.exit(f'{name}: design donors absent from the count matrix: {missing[:3]}')
        ci = [col[d] for d in donors]; genes, M = [], []
        for row in rd:
            if row[0] in passing: genes.append(row[0]); M.append([float(row[i]) for i in ci])
    M = np.array(M); eff = np.array([lib[d] for d in donors])
    if (eff <= 0).any(): sys.exit(f'{name}: non-positive effective library size')
    logcpm = np.log2((M + 0.5) / (eff[None, :] + 1.0) * 1e6)
    return {'name': name, 'genes': genes, 'pos': {g: i for i, g in enumerate(genes)}, 'logcpm': logcpm, 'donors': donors, 'design_cols': cols, 'X': X,
            'n_passing': len(passing), 'n_endo': int(X[:, cols.index(DISEASE)].sum()), 'n_non_endo': int(len(donors) - X[:, cols.index(DISEASE)].sum())}


def fit_score(score, coh, P):
    X, cols = coh['X'], coh['design_cols']
    keep = [j for j, k in enumerate(cols) if k == 'intercept' or k == DISEASE or len(np.unique(X[:, j])) > 1]
    Xk = X[:, keep]; idx = [cols[j] for j in keep].index(DISEASE)
    fit = ols(score, Xk, idx, P['minimum_residual_df'])
    fit['columns'] = [cols[j] for j in keep]; fit['dropped_constant_columns'] = [cols[j] for j in range(len(cols)) if j not in keep]
    if fit['status'] == 'ESTIMATED': fit['mde'] = float(mde_exact(fit['se'], fit['df'], P['alpha'], P['mde_power']))
    return fit, Xk, idx


def score_cohort(coh, members, P):
    """Donor-level signature scores and the disease coefficient on the v1 primary design."""
    idx_g = [coh['pos'][g] for g in members if g in coh['pos']]
    res = {'genes_in_cohort_family': len(idx_g), 'cohort_family_size': len(coh['genes']), 'n_endo': coh['n_endo'], 'n_non_endo': coh['n_non_endo']}
    if len(idx_g) < P['minimum_family_genes']: res['status'] = 'NOT_EVALUABLE'; res['label'] = 'not_labelled'; return res, None
    d = np.array([members[coh['genes'][i]]['direction'] for i in idx_g], float); b = np.array([members[coh['genes'][i]]['logfc'] for i in idx_g], float)
    Xc = coh['logcpm'][idx_g, :]; Xc = Xc - Xc.mean(axis=1, keepdims=True)
    score = (d @ Xc) / len(idx_g); proj = (b @ Xc) / float(b @ b)
    fit, Xk, di = fit_score(score, coh, P); pfit, _, _ = fit_score(proj, coh, P)
    is_case = Xk[:, di] == 1
    res.update({'status': fit['status'], 'primary_sign_weighted': fit, 'sensitivity_logfc_weighted': pfit, 'score_auc_for_disease': auc(score, is_case),
                'score_means': {'ENDO': float(score[is_case].mean()), 'non_ENDO': float(score[~is_case].mean())}})
    if fit['status'] == 'ESTIMATED':
        lodo = []
        for k in range(len(score)):
            m = np.ones(len(score), bool); m[k] = False
            keep = [j for j in range(Xk.shape[1]) if j == di or len(np.unique(Xk[m, j])) > 1]
            f = ols(score[m], Xk[m][:, keep], keep.index(di), P['minimum_residual_df'])
            if f['status'] == 'ESTIMATED': lodo.append(f['estimate'])
        res['leave_one_donor_out_range'] = [float(min(lodo)), float(max(lodo))] if lodo else None
        res['label'] = ('difference_detected_toward_published_direction' if fit['estimate'] > 0 else 'difference_detected_against_published_direction') if fit['p'] < P['per_cohort_alpha'] else 'no_detectable_difference'
    else:
        res['label'] = 'not_labelled'
    rows = [{'donor_id': dn, 'score_sign_weighted': float(score[j]), 'score_logfc_weighted': float(proj[j]), 'disease': 'ENDO' if is_case[j] else 'non-ENDO'} for j, dn in enumerate(coh['donors'])]
    return res, rows


def main(root):
    root = Path(root)
    cfg, freeze, freeze_sha = verify(root, CONFIG, SCRIPT_REL)
    P = cfg['parameters']; out = cfg['outputs']
    per_gene = read_tsv(root / cfg['inputs']['combined_per_gene']['path'])
    bins_map = {r['gene_id']: int(r['expression_bin']) for r in read_tsv(root / cfg['inputs']['discovery_mde_table']['path'])}
    genes = [r['gene_id'] for r in per_gene]; pos = {g: i for i, g in enumerate(genes)}
    missing_bin = [g for g in genes if g not in bins_map]
    if missing_bin: sys.exit(f'{len(missing_bin)} family genes lack an expression bin')
    bins = np.array([bins_map[g] for g in genes])
    t_pooled = np.array([float(r['t']) for r in per_gene]); beta = np.array([float(r['beta_ivw']) for r in per_gene])
    stat = {'pooled': np.abs(t_pooled), 'discovery': np.abs(np.array([float(r['logFC_discovery']) / float(r['SE_discovery']) for r in per_gene])),
            'heca': np.abs(np.array([float(r['logFC_heca']) / float(r['SE_heca']) for r in per_gene]))}
    sets = {}
    for r in read_tsv(root / cfg['inputs']['signature_sets']['path']):
        sets.setdefault(r['set'], {})[r['ensembl_gene_id']] = {'direction': int(r['direction']), 'logfc': float(r['published_logFC_disease_minus_control']), 'symbol': r['symbol']}
    cohorts = {name: load_cohort(root, cfg, name) for name in cfg['cohorts']}
    rng = np.random.default_rng(P['seed']); results = {}; member_rows = []; donor_rows = []
    for name in cfg['set_order']:
        full = sets.get(name, {}); members = np.array(sorted(pos[g] for g in full if g in pos), int)
        res = {'set_size_published': len(full), 'set_size_in_family': int(members.size), 'direction_split_in_family': {'up_in_disease': int(sum(full[genes[i]]['direction'] == 1 for i in members)), 'down_in_disease': int(sum(full[genes[i]]['direction'] == -1 for i in members))}}
        if members.size >= P['minimum_family_genes']:
            mask = np.zeros(len(genes), bool); mask[members] = True
            idx, order = matched_draws(bins, members, P['matched_draws'], rng)
            d = np.array([full[genes[i]]['direction'] for i in order], float)
            # T1 unsigned enrichment
            res['T1_enrichment'] = {}
            for which in ('pooled', 'discovery', 'heca'):
                s = stat[which]; w = stats.mannwhitneyu(s[mask], s[~mask], alternative='greater'); ref = s[idx].mean(axis=1); obs = float(s[mask].mean())
                res['T1_enrichment'][which] = {'wilcoxon_p_one_sided': float(w.pvalue), 'mean_abs_t_set': obs, 'mean_abs_t_rest': float(s[~mask].mean()),
                                               'matched_reference_mean': float(ref.mean()), 'matched_percentile': float((ref < obs).mean() * 100), 'matched_p_one_sided': float(((ref >= obs).sum() + 1) / (P['matched_draws'] + 1))}
            # T2 directional concordance on the pooled statistic
            signed = float((d * t_pooled[order]).mean()); ref_signed = (t_pooled[idx] * d[None, :]).mean(axis=1)
            agree = (np.sign(beta[order]) == d); share = float(agree.mean()); ref_share = (np.sign(beta[idx]) == d[None, :]).mean(axis=1)
            pct = float((ref_signed < signed).mean() * 100)
            binom_greater = stats.binomtest(int(agree.sum()), int(agree.size), 0.5, alternative='greater').pvalue
            binom_less = stats.binomtest(int(agree.sum()), int(agree.size), 0.5, alternative='less').pvalue
            res['T2_direction'] = {'signed_mean_t': signed, 'matched_reference_mean': float(ref_signed.mean()), 'matched_percentile': pct,
                                   'matched_p_two_sided': float(min(1.0, 2 * min(((ref_signed >= signed).sum() + 1), ((ref_signed <= signed).sum() + 1)) / (P['matched_draws'] + 1))),
                                   'sign_agreement_share': share, 'sign_agreement_n': int(agree.size), 'sign_agreement_matched_mean': float(ref_share.mean()),
                                   'sign_agreement_matched_percentile': float((ref_share < share).mean() * 100),
                                   'binomial_p_greater': float(binom_greater), 'binomial_p_less': float(binom_less),
                                   'spearman_published_logfc_vs_pooled_t': float(stats.spearmanr(np.array([full[genes[i]]['logfc'] for i in order]), t_pooled[order]).statistic)}
            pr = res['T1_enrichment']['pooled']
            res['label_T1'] = 'enriched' if pr['wilcoxon_p_one_sided'] < P['label_alpha'] and pr['matched_percentile'] >= P['label_percentile'] else 'no_enrichment'
            if pct >= P['label_percentile'] and share > 0.5 and binom_greater < P['label_alpha']: res['label_T2'] = 'concordant'
            elif pct <= 100 - P['label_percentile'] and share < 0.5 and binom_less < P['label_alpha']: res['label_T2'] = 'discordant'
            else: res['label_T2'] = 'no_detectable_concordance'
            for i in members:
                r = per_gene[i]; m = full[genes[i]]
                member_rows.append({'set': name, 'gene_id': genes[i], 'symbol': m['symbol'], 'published_direction': m['direction'], 'published_logFC': m['logfc'], 'beta_ivw': r['beta_ivw'], 'ci_low': r['ci_low'], 'ci_high': r['ci_high'], 'q_bh': r['q_bh'],
                                    't_pooled': float(t_pooled[i]), 'same_direction_as_published': int(np.sign(beta[i]) == m['direction']), 'logFC_discovery': r['logFC_discovery'], 'logFC_heca': r['logFC_heca'], 'expression_bin': int(bins[i])})
        else:
            res['label_T1'] = res['label_T2'] = 'not_labelled'
        # T3 donor-level scores
        res['T3_donor_score'] = {}
        for cname, coh in cohorts.items():
            r3, rows = score_cohort(coh, full, P); res['T3_donor_score'][cname] = r3
            if rows: donor_rows.extend([{'set': name, 'cohort': cname, **row} for row in rows])
        fits = [res['T3_donor_score'][c].get('primary_sign_weighted') for c in cfg['cohorts']]
        ok = [f for f in fits if f and f.get('status') == 'ESTIMATED']
        res['label_T3_pair'] = 'replicated' if len(ok) == len(fits) and len(ok) >= 2 and all(f['p'] < P['replication_alpha_each'] for f in ok) and len({np.sign(f['estimate']) for f in ok}) == 1 else ('not_labelled' if len(ok) < len(fits) else 'not_replicated')
        results[name] = res
    summary = {'analysis': cfg['analysis'], 'synthetic': bool(cfg.get('synthetic')), 'family_size': len(genes), 'matched_draws': P['matched_draws'], 'seed': P['seed'],
               'cohorts': {n: {'donors': len(c['donors']), 'n_endo': c['n_endo'], 'n_non_endo': c['n_non_endo'], 'expression_passing_genes': c['n_passing'], 'genes_in_matrix': len(c['genes'])} for n, c in cohorts.items()},
               'sets': results, 'software': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': __import__('scipy').__version__}, 'finished_utc': datetime.now(timezone.utc).isoformat()}
    outdir = root / out['directory']; outdir.mkdir(parents=True, exist_ok=True)
    with (root / out['set_results']).open('w', encoding='utf-8', newline='') as f:
        fields = ['set', 'set_size_published', 'set_size_in_family', 'test', 'statistic_or_cohort', 'value', 'matched_reference_mean', 'matched_percentile', 'p', 'label']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader()
        for name, res in results.items():
            base = {'set': name, 'set_size_published': res['set_size_published'], 'set_size_in_family': res['set_size_in_family']}
            for which, t in res.get('T1_enrichment', {}).items():
                w.writerow({**base, 'test': 'T1_mean_abs_t', 'statistic_or_cohort': which, 'value': t['mean_abs_t_set'], 'matched_reference_mean': t['matched_reference_mean'], 'matched_percentile': t['matched_percentile'], 'p': t['wilcoxon_p_one_sided'], 'label': res['label_T1']})
            if 'T2_direction' in res:
                t = res['T2_direction']
                w.writerow({**base, 'test': 'T2_signed_mean_t', 'statistic_or_cohort': 'pooled', 'value': t['signed_mean_t'], 'matched_reference_mean': t['matched_reference_mean'], 'matched_percentile': t['matched_percentile'], 'p': t['matched_p_two_sided'], 'label': res['label_T2']})
                w.writerow({**base, 'test': 'T2_sign_agreement', 'statistic_or_cohort': 'pooled', 'value': t['sign_agreement_share'], 'matched_reference_mean': t['sign_agreement_matched_mean'], 'matched_percentile': t['sign_agreement_matched_percentile'], 'p': t['binomial_p_greater'], 'label': res['label_T2']})
            for cname, r3 in res['T3_donor_score'].items():
                f3 = r3.get('primary_sign_weighted') or {}
                w.writerow({**base, 'test': 'T3_disease_coefficient', 'statistic_or_cohort': cname, 'value': f3.get('estimate', ''), 'matched_reference_mean': '', 'matched_percentile': '', 'p': f3.get('p', ''), 'label': r3['label']})
    with (root / out['set_members']).open('w', encoding='utf-8', newline='') as f:
        fields = ['set', 'gene_id', 'symbol', 'published_direction', 'published_logFC', 'beta_ivw', 'ci_low', 'ci_high', 'q_bh', 't_pooled', 'same_direction_as_published', 'logFC_discovery', 'logFC_heca', 'expression_bin']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(member_rows)
    with (root / out['donor_scores']).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['set', 'cohort', 'donor_id', 'disease', 'score_sign_weighted', 'score_logfc_weighted'], delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(donor_rows)
    (root / out['summary']).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    receipt = {'schema_version': 1, 'phase': cfg['results_phase'], 'created_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': bool(cfg.get('synthetic')),
               'protocol': cfg['protocol'], 'freeze_receipt': out['freeze_receipt'], 'freeze_sha256': freeze_sha,
               'input_sha256': {s['path']: s['sha256'] for s in cfg['inputs'].values()}, 'host_input_sha256': {s['path']: s['sha256'] for s in cfg['host_inputs'].values()},
               'script_sha256': sha(root / SCRIPT_REL), 'output_sha256': {rel: sha(root / rel) for rel in (out['set_results'], out['set_members'], out['donor_scores'], out['summary'])},
               'labels': {name: {'T1': res['label_T1'], 'T2': res['label_T2'], 'T3': {c: r['label'] for c, r in res['T3_donor_score'].items()}, 'T3_pair': res['label_T3_pair']} for name, res in results.items()},
               'software': summary['software']}
    (root / out['results_receipt']).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    s = main(Path(a.root).resolve() if a.root else ROOT)
    print(json.dumps({n: {'in_family': r['set_size_in_family'], 'T1': r['label_T1'], 'T2': r['label_T2'], 'T2_pct': (r.get('T2_direction') or {}).get('matched_percentile'), 'share': (r.get('T2_direction') or {}).get('sign_agreement_share'),
                          'T3': {c: (x.get('primary_sign_weighted') or {}).get('estimate') for c, x in r['T3_donor_score'].items()}, 'T3_labels': {c: x['label'] for c, x in r['T3_donor_score'].items()}, 'pair': r['label_T3_pair']} for n, r in s['sets'].items()}, indent=1))
