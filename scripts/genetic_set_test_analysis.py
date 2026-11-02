#!/usr/bin/env python3
"""Genetic support and the stromal disease contrast, per docs/GENETIC_SET_TEST_PROTOCOL.md.

Reads the published gene sets, the combined-bound per-gene table (pooled and per-cohort t
statistics) and the discovery MDE table (expression quintile), verifies the freeze, and for
each set computes the one-sided Wilcoxon rank-sum test of |t| (set versus the rest of the
family) and the percentile of the set's mean |t| among expression-matched random sets, for
the pooled statistic (primary) and each cohort (secondary), plus the sign-agreement share.
Writes tables, summary, receipt. Runs once. ``--root`` targets another project tree (tests).
"""
from pathlib import Path
import argparse, csv, gzip, json, platform, sys
from datetime import datetime, timezone
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from followon_freeze import verify, sha  # noqa: E402

CONFIG = 'config/genetic_set_test.json'
SCRIPT_REL = 'scripts/genetic_set_test_analysis.py'


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))


def matched_reference(stat, bins, members, draws, rng):
    """Mean |t| of random sets drawn with the member set's expression-bin composition."""
    by_bin = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    need = {b: int((bins[members]).tolist().count(b)) for b in by_bin}
    means = np.empty(draws)
    for i in range(draws):
        idx = np.concatenate([rng.choice(by_bin[b], need[b], replace=False) for b in by_bin if need[b]])
        means[i] = stat[idx].mean()
    return means


def main(root):
    root = Path(root)
    cfg, freeze, freeze_sha = verify(root, CONFIG, SCRIPT_REL)
    P = cfg['parameters']; out = cfg['outputs']
    per_gene = read_tsv(root / cfg['inputs']['combined_per_gene']['path'])
    mde_disc = read_tsv(root / cfg['inputs']['discovery_mde_table']['path'])
    sets_rows = read_tsv(root / cfg['inputs']['gene_sets']['path'])
    genes = [r['gene_id'] for r in per_gene]; pos = {g: i for i, g in enumerate(genes)}
    bins_map = {r['gene_id']: int(r['expression_bin']) for r in mde_disc}
    missing_bin = [g for g in genes if g not in bins_map]
    if missing_bin: sys.exit(f'{len(missing_bin)} family genes lack an expression bin')
    bins = np.array([bins_map[g] for g in genes])
    stat = {'pooled': np.abs(np.array([float(r['t']) for r in per_gene])),
            'discovery': np.abs(np.array([float(r['logFC_discovery']) / float(r['SE_discovery']) for r in per_gene])),
            'heca': np.abs(np.array([float(r['logFC_heca']) / float(r['SE_heca']) for r in per_gene]))}
    same_sign = np.array([r['same_sign'] == 'TRUE' for r in per_gene])
    sets = {}
    for r in sets_rows: sets.setdefault(r['set'], set()).add(r['ensembl_gene_id'])
    rng = np.random.default_rng(P['seed']); results = {}; member_rows = []
    for name in cfg['set_order']:
        full = sets.get(name, set()); members = np.array(sorted(pos[g] for g in full if g in pos), int)
        res = {'set_size_source': len(full), 'set_size_in_family': int(members.size), 'tests': {}}
        if members.size >= P['minimum_family_genes']:
            mask = np.zeros(len(genes), bool); mask[members] = True
            for which in ('pooled', 'discovery', 'heca'):
                s = stat[which]
                w = stats.mannwhitneyu(s[mask], s[~mask], alternative='greater')
                ref = matched_reference(s, bins, members, P['matched_draws'], rng)
                obs = float(s[mask].mean())
                res['tests'][which] = {'wilcoxon_p_one_sided': float(w.pvalue), 'mean_abs_t_set': obs, 'mean_abs_t_rest': float(s[~mask].mean()),
                                       'matched_reference_mean': float(ref.mean()), 'matched_percentile': float((ref < obs).mean() * 100),
                                       'matched_p_one_sided': float(((ref >= obs).sum() + 1) / (P['matched_draws'] + 1))}
            ref_sign = matched_reference(same_sign.astype(float), bins, members, P['matched_draws'], rng)
            obs_sign = float(same_sign[mask].mean())
            res['sign_agreement'] = {'set_share': obs_sign, 'matched_reference_mean': float(ref_sign.mean()), 'matched_percentile': float((ref_sign < obs_sign).mean() * 100)}
            pr = res['tests']['pooled']
            res['label'] = 'enriched' if pr['wilcoxon_p_one_sided'] < P['label_alpha'] and pr['matched_percentile'] >= P['label_percentile'] else 'no_enrichment'
            for i in members:
                r = per_gene[i]
                member_rows.append({'set': name, 'gene_id': genes[i], 'beta_ivw': r['beta_ivw'], 'ci_low': r['ci_low'], 'ci_high': r['ci_high'], 'q_bh': r['q_bh'],
                                    'abs_t_pooled': float(stat['pooled'][i]), 'logFC_discovery': r['logFC_discovery'], 'logFC_heca': r['logFC_heca'], 'same_sign': r['same_sign'], 'expression_bin': int(bins[i])})
        else:
            res['label'] = 'not_labelled'
        results[name] = res
    summary = {'analysis': cfg['analysis'], 'synthetic': bool(cfg.get('synthetic')), 'family_size': len(genes), 'matched_draws': P['matched_draws'], 'seed': P['seed'],
               'sets': results, 'software': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': __import__('scipy').__version__},
               'finished_utc': datetime.now(timezone.utc).isoformat()}
    outdir = root / out['directory']; outdir.mkdir(parents=True, exist_ok=True)
    with (root / out['set_results']).open('w', encoding='utf-8', newline='') as f:
        fields = ['set', 'set_size_source', 'set_size_in_family', 'statistic', 'wilcoxon_p_one_sided', 'mean_abs_t_set', 'mean_abs_t_rest', 'matched_reference_mean', 'matched_percentile', 'matched_p_one_sided', 'label']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader()
        for name, res in results.items():
            for which, t in res.get('tests', {}).items():
                w.writerow({'set': name, 'set_size_source': res['set_size_source'], 'set_size_in_family': res['set_size_in_family'], 'statistic': which, **t, 'label': res['label']})
    with (root / out['set_members']).open('w', encoding='utf-8', newline='') as f:
        fields = ['set', 'gene_id', 'beta_ivw', 'ci_low', 'ci_high', 'q_bh', 'abs_t_pooled', 'logFC_discovery', 'logFC_heca', 'same_sign', 'expression_bin']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(member_rows)
    (root / out['summary']).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    receipt = {'schema_version': 1, 'phase': cfg['results_phase'], 'created_utc': datetime.now(timezone.utc).isoformat(), 'synthetic': bool(cfg.get('synthetic')),
               'protocol': cfg['protocol'], 'freeze_receipt': out['freeze_receipt'], 'freeze_sha256': freeze_sha,
               'input_sha256': {s['path']: s['sha256'] for s in cfg['inputs'].values()}, 'script_sha256': sha(root / SCRIPT_REL),
               'output_sha256': {rel: sha(root / rel) for rel in (out['set_results'], out['set_members'], out['summary'])},
               'labels': {name: res['label'] for name, res in results.items()}, 'software': summary['software']}
    (root / out['results_receipt']).write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    return summary


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--root'); a = ap.parse_args()
    s = main(Path(a.root).resolve() if a.root else ROOT)
    print(json.dumps({name: {'n': r['set_size_in_family'], 'label': r['label'], 'pooled': r.get('tests', {}).get('pooled', {}).get('wilcoxon_p_one_sided')} for name, r in s['sets'].items()}))
