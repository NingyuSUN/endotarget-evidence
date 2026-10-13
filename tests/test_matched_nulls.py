"""Prespecified null sampler hand cases; no real effects/external inputs."""
import collections
import copy
import csv
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import sample_matched_nulls as nulls


def gid(i):
    return 'ENSG' + str(i).zfill(11)


def rows(n=80):
    return [{'gene_id': gid(i), 'null_matching_logCPM': 3., 'weighted_degree': 1.4,
             'logFC': [0., -1., 1., 1e-300][i % 4]} for i in range(1, n + 1)]


def selections(data):
    return {m: [r['gene_id'] for r in data[j:j+20]]
            for m, j in [('B0', 0), ('B1', 20), ('G1', 40)]}


class MatchingTests(unittest.TestCase):
    def test_type7_hand_boundaries_and_equal_goes_lower(self):
        cuts, bins = nulls.axis_bins([0, 1, 2, 3, 4, 5])
        self.assertEqual(cuts.tolist(), [1., 2., 3., 4.])
        self.assertEqual(bins.tolist(), [1, 1, 2, 3, 4, 5])
        cuts, _ = nulls.axis_bins([0, 10])
        np.testing.assert_allclose(cuts, [2, 4, 6, 8])

    def test_duplicate_boundaries_identical_values_not_split(self):
        cuts, bins = nulls.axis_bins([2] * 50)
        self.assertEqual(cuts.tolist(), [2.])
        self.assertEqual(bins.tolist(), [1] * 50)
        cuts, bins = nulls.axis_bins([0] * 8 + [1] * 2)
        np.testing.assert_allclose(cuts, [0., .2], rtol=0, atol=1e-15)
        self.assertEqual(bins.tolist(), [1] * 8 + [3] * 2)

    def test_exact_200_sets_20_unique_own_exclusions_own_direction(self):
        data = rows()
        selected = selections(data)
        result = nulls.sample_methods(data, selected)
        reference = {r['gene_id']: r for r in data}
        for method in nulls.METHODS:
            status = result['methods'][method]
            self.assertEqual(status['status'], 'ESTIMABLE')
            self.assertEqual(status['sampled_sets'], 200)
            members = result['memberships'][method]
            self.assertEqual(len(members), 4000)
            for draw in range(1, 201):
                sample = [r for r in members if r['draw'] == draw]
                self.assertEqual(len(set(r['gene_id'] for r in sample)), 20)
                self.assertFalse(set(r['gene_id'] for r in sample) & set(selected[method]))
                self.assertEqual({(r['expression_bin'], r['degree_bin']) for r in sample}, {(1, 1)})
                for row in sample:
                    x = reference[row['gene_id']]['logFC']
                    self.assertEqual(row['discovery_sign'], (x > 0) - (x < 0))
            # Other methods' targets remain in eligible candidate pool.
            self.assertTrue(set(r['gene_id'] for r in members) & set(selected['B1' if method == 'B0' else 'B0']))
        self.assertEqual({r['discovery_sign'] for r in result['universe_bins']}, {-1, 0, 1})

    def test_joint_occupancies_every_draw_not_only_margins(self):
        data = []
        for expression in range(5):
            for degree in range(5):
                for _ in range(20):
                    data.append({'gene_id': gid(len(data)+1), 'null_matching_logCPM': expression,
                                 'weighted_degree': degree+1., 'logFC': -1.})
        chosen = [data[j*20]['gene_id'] for j in range(20)]
        selected = {m: chosen for m in nulls.METHODS}
        result = nulls.sample_methods(data, selected)
        bins = {r['gene_id']: (r['expression_bin'], r['degree_bin']) for r in result['universe_bins']}
        demand = collections.Counter(bins[g] for g in chosen)
        for method in nulls.METHODS:
            for draw in range(1, 201):
                self.assertEqual(collections.Counter((r['expression_bin'], r['degree_bin'])
                    for r in result['memberships'][method] if r['draw'] == draw), demand)

    def test_rng_separately_reinitialized_stable_sorted_candidates(self):
        data = rows()
        selected = {m: [r['gene_id'] for r in data[:20]] for m in nulls.METHODS}
        first = nulls.sample_methods(data, selected)
        second = nulls.sample_methods(list(reversed(data)), selected)
        self.assertEqual(first, second)
        for m in ['B1', 'G1']:
            self.assertEqual([r['gene_id'] for r in first['memberships']['B0']],
                             [r['gene_id'] for r in first['memberships'][m]])
        # Independent executable RNG hand-check, same stable pool and PCG64 choice.
        expected = np.random.Generator(np.random.PCG64(42)).choice(
            sorted(r['gene_id'] for r in data[20:]), size=20, replace=False).tolist()
        self.assertEqual([r['gene_id'] for r in first['memberships']['B0'][:20]], expected)

    def test_no_effect_sign_or_t_or_external_fields_select_memberships(self):
        data = rows()
        selected = selections(data)
        original = nulls.sample_methods(data, selected)
        changed = copy.deepcopy(data)
        for r in changed:
            r.update(logFC=-r['logFC'], t=999, **{'P.Value': .0001,
                     'external_available': False, 'external_logFC': 999})
        updated = nulls.sample_methods(changed, selected)
        for m in nulls.METHODS:
            self.assertEqual([r['gene_id'] for r in original['memberships'][m]],
                             [r['gene_id'] for r in updated['memberships'][m]])

    def test_one_scarce_bin_fails_entire_method_no_retry_or_refill(self):
        data = rows()
        # 20 high-degree genes form the final quintile and B0 selects them all.
        for r in data[-20:]: r['weighted_degree'] = 10.
        selected = selections(data)
        selected['B0'] = [r['gene_id'] for r in data[-20:]]
        result = nulls.sample_methods(data, selected)
        status = result['methods']['B0']
        self.assertEqual(status['status'], 'NOT_ESTIMABLE')
        self.assertEqual(status['sampled_sets'], 0)
        self.assertEqual(result['memberships']['B0'], [])
        scarce = [b for b in status['bins'] if b['insufficient']]
        self.assertEqual([(b['demand'], b['available']) for b in scarce], [(20, 0)])
        self.assertEqual(result['methods']['B1']['status'], 'ESTIMABLE')

    def test_invalid_inputs_fail_instead_of_repair(self):
        data = rows()
        bad = copy.deepcopy(data); bad[1]['gene_id'] = bad[0]['gene_id']
        with self.assertRaises(ValueError): nulls.sample_methods(bad, selections(data))
        for field, value in [('null_matching_logCPM', float('nan')), ('weighted_degree', 0), ('logFC', float('inf'))]:
            bad = copy.deepcopy(data); bad[0][field] = value
            with self.assertRaises(ValueError): nulls.sample_methods(bad, selections(data))
        for selected in [{m: [gid(999)]*20 for m in nulls.METHODS},
                         {m: [gid(i) for i in range(1,20)] for m in nulls.METHODS}]:
            with self.assertRaises(ValueError): nulls.sample_methods(data, selected)


class OutputTests(unittest.TestCase):
    def fixture(self, directory):
        import rank_targets
        de = directory / 'de.tsv'
        with de.open('w') as h:
            h.write('gene_id\tlogFC\tt\tP.Value\tadj.P.Val\tAveExpr\tnull_matching_logCPM\n')
            for i in range(1,81):
                h.write(f'{gid(i)}\t{[0,-1,1][i%3]}\t{81-i}\t0.1\t0.2\t2\t3\n')
        edges = directory / 'edges.tsv'
        with edges.open('w') as h:
            h.write('gene1\tgene2\tweight\n')
            for i in range(1,81): h.write(f'{gid(i)}\t{gid(i%80+1)}\t0.7\n')
        rankings = directory / 'ranks'
        rank_targets.execute(de, edges, rankings, ROOT/'config/ranking_algorithm.json', synthetic_test=True)
        return de, rankings

    def test_cli_archives_all_bins_directions_hashes_and_no_external(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'work') as t:
            directory = pathlib.Path(t)
            de, rankings = self.fixture(directory)
            out = directory/'nulls'
            cmd = [sys.executable, str(ROOT/'scripts/sample_matched_nulls.py'), '--de', str(de),
                   '--ranking-dir', str(rankings), '--output-dir', str(out), '--synthetic-test']
            done = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            status = json.loads((out/'status.json').read_text())
            self.assertEqual(status['status'], 'SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE')
            self.assertFalse(status['complete_external_freeze'])
            self.assertFalse(status['external_data_used'])
            self.assertEqual(len(list(nulls.read_table(out/'universe_bins.tsv', ['gene_id']))), 80)
            self.assertEqual(len(list(nulls.read_table(out/'all_U_directions.tsv', ['gene_id']))), 80)
            self.assertEqual(len(list(nulls.read_table(out/'null_memberships.tsv', ['gene_id']))), 12000)
            for name, receipt in status['outputs'].items():
                self.assertEqual(receipt['sha256'], nulls.digest(out/name))
            self.assertEqual(status['methods']['B0']['sampled_sets'], 200)
            self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)

    def test_production_rejects_synthetic_rank_and_changed_rank_or_de(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'work') as t:
            directory=pathlib.Path(t); de, rankings=self.fixture(directory)
            with self.assertRaises(ValueError):
                nulls.execute(de, rankings, directory/'real', ROOT/'config/statistical_analysis.json')
            with (rankings/'B0_top20.tsv').open('a') as h: h.write('broken\n')
            with self.assertRaises(ValueError):
                nulls.execute(de, rankings, directory/'changed', ROOT/'config/statistical_analysis.json', synthetic_test=True)

    def test_changed_discovery_source_fails_before_any_output(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'work') as t:
            directory=pathlib.Path(t); de, rankings=self.fixture(directory)
            de.write_text(de.read_text().replace('\t3\n', '\t4\n'))
            with self.assertRaisesRegex(ValueError, 'Discovery DE differs'):
                nulls.execute(de, rankings, directory/'changed', ROOT/'config/statistical_analysis.json', synthetic_test=True)
            self.assertFalse((directory/'changed').exists())

    def test_config_drift_rejected_and_cli_has_no_sampling_override(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'work') as t:
            directory=pathlib.Path(t); de, rankings=self.fixture(directory)
            cfg=json.loads((ROOT/'config/statistical_analysis.json').read_text())
            cfg['null']['seed']=43
            path=directory/'bad.json'; path.write_text(json.dumps(cfg))
            with self.assertRaises(ValueError):
                nulls.execute(de, rankings, directory/'bad', path, synthetic_test=True)
        helptext=subprocess.run([sys.executable,str(ROOT/'scripts/sample_matched_nulls.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(helptext.returncode,0)
        self.assertNotIn('--external',helptext.stdout)
        self.assertNotIn('--seed',helptext.stdout)
        self.assertNotIn('--n-sets',helptext.stdout)


if __name__ == '__main__': unittest.main()
