#!/usr/bin/env python3
"""Freeze discovery-only expression/degree matched gene sets before external effects."""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import math
import pathlib
import platform
import re
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
METHODS = ('B0', 'B1', 'G1')
GENE = re.compile(r'^ENSG[0-9]{11}$')
NULL_CONFIG = {
    'n_sets': 200, 'k': 20, 'seed': 42, 'rng': 'numpy_PCG64',
    'method_order': list(METHODS),
    'expression': 'rowMeans_edgeR_cpm_log_TRUE_prior_count_0.5_normalized_lib_sizes_TRUE',
    'degree': 'weighted_degree_in_frozen_discovery_induced_U',
    'bins_per_axis': 5, 'quantile_type': 7,
    'collapse_duplicate_quantile_boundaries': True,
    'boundary_assignment': 'internal_cutpoints_equal_values_go_lower',
    'joint_bins': True, 'pool': 'U_excluding_own_method_top20',
    'external_mask_used_for_sampling': False,
    'discovery_sign_for_each_sampled_gene': True,
    'within_set_replacement': False, 'between_set_reuse': True,
    'bin_failure_policy': 'entire_method_null_not_estimable_no_merging_no_rescue',
    'descriptive_only': True,
}
BIN_FIELDS = ['gene_id', 'null_matching_logCPM', 'weighted_degree', 'expression_bin',
              'degree_bin', 'logFC', 'discovery_sign']
MEMBER_FIELDS = ['method', 'draw', 'position', 'gene_id', 'expression_bin', 'degree_bin',
                 'null_matching_logCPM', 'weighted_degree', 'logFC', 'discovery_sign']
DEMAND_FIELDS = ['method', 'expression_bin', 'degree_bin', 'demand', 'available', 'insufficient']


def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(1024*1024), b''): h.update(part)
    return h.hexdigest()


def read_table(path, fields):
    path = pathlib.Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        header = reader.fieldnames
        if not header or len(header) != len(set(header)) or not set(fields).issubset(header):
            raise ValueError(f'Missing or duplicate TSV columns: {path}')
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'Ragged TSV row {line}: {path}')
            yield row


def number(value, field):
    try: value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'Invalid {field}') from exc
    if not math.isfinite(value): raise ValueError(f'Nonfinite {field}')
    return value


def sign(value):
    return (value > 0) - (value < 0)


def axis_bins(values):
    """Type-7 pooled quintiles; duplicates collapse; internal equal values go lower."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or len(values) == 0 or not np.all(np.isfinite(values)):
        raise ValueError('Binning requires a nonempty finite vector')
    cuts = np.unique(np.quantile(values, [0.2, 0.4, 0.6, 0.8], method='linear'))
    return cuts, np.searchsorted(cuts, values, side='left') + 1


def prepare_universe(rows):
    genes = {}
    for row in rows:
        gid = row.get('gene_id')
        if not isinstance(gid, str) or not GENE.fullmatch(gid) or gid in genes:
            raise ValueError('Invalid or duplicate discovery U gene_id')
        item = {'gene_id': gid}
        for field in ('null_matching_logCPM', 'weighted_degree', 'logFC'):
            item[field] = number(row.get(field), f'{gid}/{field}')
        if item['weighted_degree'] <= 0: raise ValueError('Discovery U degree must be positive')
        item['discovery_sign'] = sign(item['logFC'])
        genes[gid] = item
    if len(genes) < 50: raise ValueError('Frozen discovery U must contain at least 50 genes')
    ordered = [genes[g] for g in sorted(genes)]
    expression, ebins = axis_bins([g['null_matching_logCPM'] for g in ordered])
    degree, dbins = axis_bins([g['weighted_degree'] for g in ordered])
    for row, e, d in zip(ordered, ebins, dbins):
        row.update(expression_bin=int(e), degree_bin=int(d))
    return ordered, {'expression_internal_cutpoints': expression.tolist(),
                     'degree_internal_cutpoints': degree.tolist(),
                     'quantile_probabilities': [0.2, 0.4, 0.6, 0.8],
                     'quantile_method': 'numpy_linear_R_type7',
                     'assignment': '1+count(value>unique_internal_cutpoints)',
                     'pooled_universe_genes': len(ordered)}


def sample_methods(rows, selected):
    """Fixed 200x20 sampler; membership ignores effect/sign/external availability."""
    universe, boundaries = prepare_universe(rows)
    by_id = {r['gene_id']: r for r in universe}
    if set(selected) != set(METHODS): raise ValueError('Exactly B0/B1/G1 top20 required')
    methods, memberships = {}, {}
    for method in METHODS:
        own = selected[method]
        if len(own) != 20 or len(set(own)) != 20 or not set(own).issubset(by_id):
            raise ValueError(f'{method}: top20 requires 20 unique members of frozen U')
        own = set(own)
        occupancy = collections.Counter((by_id[g]['expression_bin'], by_id[g]['degree_bin']) for g in own)
        candidates = collections.defaultdict(list)
        for row in universe:  # Already stable-ID sorted, independent of input order.
            if row['gene_id'] not in own:
                candidates[(row['expression_bin'], row['degree_bin'])].append(row['gene_id'])
        bin_receipts = [{'expression_bin': e, 'degree_bin': d, 'demand': occupancy.get((e,d), 0),
                         'available': len(candidates[(e,d)]),
                         'insufficient': len(candidates[(e,d)]) < occupancy.get((e,d), 0)}
                        for e, d in sorted(set(candidates) | set(occupancy))]
        failed = any(b['insufficient'] for b in bin_receipts)
        methods[method] = {'status': 'NOT_ESTIMABLE' if failed else 'ESTIMABLE',
                           'reason': 'INSUFFICIENT_OWN_EXCLUDED_JOINT_BIN_POOL' if failed else None,
                           'requested_sets': 200, 'sampled_sets': 0 if failed else 200,
                           'genes_per_set': 20, 'bins': bin_receipts,
                           'excluded_own_top20': sorted(own), 'candidate_pool_genes': len(universe)-20}
        members = []
        if not failed:
            rng = np.random.Generator(np.random.PCG64(42))  # Fresh stream EACH method.
            for draw in range(1, 201):
                position = 0
                for joint in sorted(occupancy):
                    chosen = rng.choice(candidates[joint], size=occupancy[joint], replace=False)
                    for gid in chosen.tolist():
                        position += 1
                        members.append({'method': method, 'draw': draw, 'position': position, **by_id[gid]})
                if position != 20: raise RuntimeError('Internal null-set size invariant failed')
        memberships[method] = members
    return {'universe_bins': universe, 'boundaries': boundaries, 'methods': methods,
            'memberships': memberships}


def load_inputs(de_path, ranking_dir, synthetic_test):
    """Validate and hash the upstream rank freeze; never silently alter its universe."""
    names = ['universe.tsv', 'induced_edges.tsv', 'provenance.json']
    names += [f'{m}{suffix}.tsv' for m in METHODS for suffix in ('', '_top20')]
    paths = {'discovery_de': de_path, **{name: ranking_dir/name for name in names}}
    hashes = {name: digest(path) for name, path in paths.items()}
    provenance = json.loads((ranking_dir/'provenance.json').read_text())
    expected_status = ('SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE' if synthetic_test
                       else 'DISCOVERY_RANKING_ARTIFACT_NOT_COMPLETE_EXTERNAL_FREEZE')
    if provenance.get('status') != expected_status or provenance.get('external_data_used') is not False:
        raise ValueError('Ranking provenance is not the required discovery-only production/synthetic state')
    if provenance.get('inputs', {}).get('discovery_de', {}).get('sha256') != hashes['discovery_de']:
        raise ValueError('Discovery DE differs from ranking source hash')
    for name in names:
        if name == 'provenance.json': continue
        if provenance.get('outputs', {}).get(name, {}).get('sha256') != hashes[name]:
            raise ValueError(f'Frozen ranking output changed: {name}')
    urows = list(read_table(ranking_dir/'universe.tsv', ['gene_id']))
    universe = [r['gene_id'] for r in urows]
    if universe != sorted(set(universe)) or len(universe) < 50:
        raise ValueError('Ranking U must be sorted, unique and have >=50 genes')
    de = {}
    for row in read_table(de_path, ['gene_id', 'logFC', 'null_matching_logCPM']):
        if not GENE.fullmatch(row['gene_id']) or row['gene_id'] in de:
            raise ValueError('Invalid/duplicate discovery DE gene_id')
        # DE outside U stays outside U; all provided fit values must still be finite.
        de[row['gene_id']] = {'gene_id': row['gene_id'],
                             'logFC': number(row['logFC'], 'discovery logFC'),
                             'null_matching_logCPM': number(row['null_matching_logCPM'], 'matching expression')}
    if not set(universe).issubset(de): raise ValueError('Discovery U gene missing in discovery DE')
    degree = dict.fromkeys(universe, 0.)
    pairs = set()
    for row in read_table(ranking_dir/'induced_edges.tsv', ['gene1', 'gene2', 'weight']):
        a, b = row['gene1'], row['gene2']; pair=tuple(sorted((a,b)))
        w = number(row['weight'], 'induced weight')
        if a not in degree or b not in degree or a == b or pair in pairs or not .7 <= w <= 1:
            raise ValueError('Invalid induced graph endpoints/duplicate/self-loop/weight')
        pairs.add(pair); degree[a] += w; degree[b] += w
    selected = {}
    for method in METHODS:
        fields = ['rank', 'gene_id', 'method', 'logFC', 'discovery_sign', 'weighted_degree']
        full = list(read_table(ranking_dir/f'{method}.tsv', fields))
        top = list(read_table(ranking_dir/f'{method}_top20.tsv', fields))
        if len(full) != len(universe) or set(r['gene_id'] for r in full) != set(universe):
            raise ValueError(f'{method} does not preserve exact frozen U')
        if len(top) != 20 or top != full[:20]: raise ValueError(f'{method} top20 differs from rank prefix')
        for i, row in enumerate(full, 1):
            gid = row['gene_id']
            if row['rank'] != str(i) or row['method'] != method:
                raise ValueError(f'{method}: invalid rank/method')
            if number(row['logFC'], 'rank logFC') != de[gid]['logFC'] or row['discovery_sign'] != str(sign(de[gid]['logFC'])):
                raise ValueError(f'{method}: discovery coefficient/sign changed')
            observed = number(row['weighted_degree'], 'rank weighted degree')
            if not math.isclose(observed, degree[gid], rel_tol=1e-13, abs_tol=1e-13):
                raise ValueError(f'{method}: degree differs from frozen induced graph')
            # Use original ranking float64 covariate; graph summation order may differ.
            if method == 'B0': de[gid]['weighted_degree'] = observed
            elif observed != de[gid]['weighted_degree']:
                raise ValueError('Methods have different weighted degree covariates')
        selected[method] = [r['gene_id'] for r in top]
    return [de[g] for g in universe], selected, paths, hashes


def write_table(path, rows, fields):
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({f: format(row[f], '.17g') if isinstance(row[f], float) else row[f] for f in fields})


def execute(de_path, ranking_dir, output_dir, config_path, *, synthetic_test=False):
    de_path, ranking_dir, output_dir, config_path = map(pathlib.Path, (de_path, ranking_dir, output_dir, config_path))
    partial = output_dir.with_name(output_dir.name+'.partial')
    if output_dir.exists() or partial.exists(): raise ValueError('Existing output/partial preserved; use a new directory')
    config_hash = digest(config_path); code_hash = digest(__file__)
    config = json.loads(config_path.read_text())
    if config.get('null') != NULL_CONFIG: raise ValueError('Null configuration differs from frozen contract')
    rows, selected, paths, hashes = load_inputs(de_path, ranking_dir, synthetic_test)
    paths.update(configuration=config_path, implementation=pathlib.Path(__file__))
    hashes.update(configuration=config_hash, implementation=code_hash)
    result = sample_methods(rows, selected)
    for name, path in paths.items():
        if digest(path) != hashes[name]: raise ValueError(f'Input changed during sampling: {name}')
    partial.mkdir(parents=True, exist_ok=False)
    write_table(partial/'universe_bins.tsv', result['universe_bins'], BIN_FIELDS)
    write_table(partial/'all_U_directions.tsv', result['universe_bins'], ['gene_id', 'logFC', 'discovery_sign'])
    write_table(partial/'null_memberships.tsv', [row for m in METHODS for row in result['memberships'][m]], MEMBER_FIELDS)
    (partial/'boundaries.json').write_text(json.dumps(result['boundaries'], indent=2, allow_nan=False)+'\n')
    demands = []
    for method in METHODS:
        write_table(partial/f'{method}_null_membership.tsv', result['memberships'][method], MEMBER_FIELDS)
        demands += [{'method': method, **b} for b in result['methods'][method]['bins']]
    write_table(partial/'bin_demands.tsv', demands, DEMAND_FIELDS)
    outputs = {p.name: {'sha256': digest(p), 'bytes': p.stat().st_size} for p in sorted(partial.iterdir())}
    status = {
        'schema_version': 1,
        'status': ('SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE' if synthetic_test
                   else 'DISCOVERY_NULL_MEMBERSHIPS_FROZEN_NOT_COMPLETE_EXTERNAL_FREEZE'),
        'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'inputs': {name: {'path': str(path.resolve()), 'sha256': hashes[name]} for name, path in paths.items()},
        'software': {'python': platform.python_version(), 'numpy': np.__version__, 'executable': sys.executable},
        'algorithm': NULL_CONFIG, 'methods': result['methods'], 'outputs': outputs,
        'universe_genes': len(result['universe_bins']), 'external_data_used': False,
        'complete_external_freeze': False, 'scientific_validation': False,
        'interpretation': 'descriptive_expression_degree_sensitivity_not_permutation_p_value_or_causal_control',
    }
    (partial/'status.json').write_text(json.dumps(status, indent=2, allow_nan=False)+'\n')
    partial.rename(output_dir)
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--de', required=True, type=pathlib.Path, help='Same discovery fitted DE as rankings; includes null_matching_logCPM')
    parser.add_argument('--ranking-dir', required=True, type=pathlib.Path)
    parser.add_argument('--output-dir', required=True, type=pathlib.Path)
    parser.add_argument('--config', type=pathlib.Path, default=ROOT/'config/statistical_analysis.json')
    parser.add_argument('--synthetic-test', action='store_true', help='Synthetic ranking provenance required; never a scientific freeze')
    args = parser.parse_args()
    try:
        status = execute(args.de, args.ranking_dir, args.output_dir, args.config, synthetic_test=args.synthetic_test)
    except (ValueError, OSError, json.JSONDecodeError, OverflowError) as exc:
        parser.exit(1, f'Null sampling failed: {exc}\n')
    print(json.dumps({'status': status['status'], 'output_dir': str(args.output_dir),
                      'methods': {m: {k: status['methods'][m][k] for k in ['status','sampled_sets','reason']} for m in METHODS}}, indent=2))


if __name__ == '__main__': main()
