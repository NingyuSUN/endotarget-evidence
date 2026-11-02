# Synthetic boundary checks for the two-cohort combined bound: estimator identities against
# hand calculations, a planted-effect fixture end to end, and the freeze/pin refusals.
import csv, gzip, hashlib, importlib.util, json, pathlib, shutil, subprocess, sys, tempfile, unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / 'scripts/combined_bound_analysis.py'; FREEZE = ROOT / 'scripts/combined_bound_freeze.py'
spec = importlib.util.spec_from_file_location('combined_bound', ANALYSIS)
cb = importlib.util.module_from_spec(spec); spec.loader.exec_module(cb)
import numpy as np
from scipy import stats


class EstimatorIdentities(unittest.TestCase):
    def test_inverse_variance_hand_values(self):
        beta, se, q = cb.ivw(np.array([[1.0, 3.0]]), np.array([[1.0, 2.0]]))
        self.assertAlmostEqual(beta[0], (1.0 * 1 + 0.25 * 3) / 1.25); self.assertAlmostEqual(se[0], 1 / np.sqrt(1.25))
        self.assertAlmostEqual(q[0], 1.0 * (1 - beta[0]) ** 2 + 0.25 * (3 - beta[0]) ** 2)

    def test_satterthwaite_equal_precision_adds_degrees_of_freedom(self):
        nu = cb.satterthwaite_df(np.array([[0.3, 0.3]]), np.array([[12.0, 15.0]]))
        self.assertAlmostEqual(nu[0], 4 / (1 / 12 + 1 / 15))
        nu2 = cb.satterthwaite_df(np.array([[0.3, 0.3]]), np.array([[10.0, 10.0]])); self.assertAlmostEqual(nu2[0], 20.0)

    def test_bh_matches_textbook_example(self):
        # ranked: 0.01*4/1 = 0.04, 0.03*4/2 = 0.06, 0.04*4/3 = 0.0533, 0.2; the step-up minimum gives 0.0533 to both middle genes
        q = cb.bh([0.01, 0.04, 0.03, 0.2]); np.testing.assert_allclose(q, [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.2])
        self.assertTrue((np.diff(np.sort(q)) >= 0).all())

    def test_exact_mde_has_the_requested_power_and_matches_the_approximation_at_alpha_0_05(self):
        se, df = 0.2, 25.0
        d = cb.mde_exact(se, df, 0.05, 0.8); tcrit = stats.t.ppf(0.975, df)
        power = stats.nct.sf(tcrit, df, d / se) + stats.nct.cdf(-tcrit, df, d / se)
        self.assertAlmostEqual(power, 0.8, places=6)
        self.assertLess(abs(d - cb.mde_approx(se, df, 0.05, 0.8)), 0.01)
        tiny = cb.mde_exact(se, df, 1e-6, 0.8); self.assertGreater(tiny, d)
        self.assertAlmostEqual(float(cb.nct_power(stats.t.ppf(1 - 5e-7, df), df, tiny / se)), 0.8, places=6)
        vec = cb.mde_exact_vec(np.array([0.2, 0.5]), np.array([25.0, 12.6]), 0.05, 0.8)
        self.assertAlmostEqual(vec[0], d); self.assertGreater(vec[1], vec[0])

    def test_stouffer_with_equal_inputs_reduces_to_sqrt_two_scaling(self):
        z, p = cb.stouffer(np.array([[0.5, 0.5]]), np.array([[0.25, 0.25]]), np.array([[20.0, 20.0]]))
        z_each = stats.norm.isf(stats.t.sf(2.0, 20.0)); self.assertAlmostEqual(z[0], np.sqrt(2) * z_each); self.assertLess(p[0], 0.05)


def write_table(path, rows):
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 'logFC', 'moderated_SE', 'df_total', 'fit_estimable'])
        for r in rows: w.writerow(r)


class EndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'project'
        for d in ('config', 'scripts', 'tests', 'docs', 'receipts', 'reports/mde_analysis'): (self.root / d).mkdir(parents=True)
        for rel in ('scripts/combined_bound_analysis.py', 'scripts/combined_bound_freeze.py', 'tests/test_combined_bound.py', 'docs/COMBINED_BOUND_PROTOCOL.md'):
            if (ROOT / rel).exists(): shutil.copyfile(ROOT / rel, self.root / rel)
        rng = np.random.default_rng(7); n = 400
        self.planted = [f'G{i:04d}' for i in range(5)]
        rows = {'discovery': [], 'heca': []}
        for i in range(n):
            g = f'G{i:04d}'; true = 1.6 if g in self.planted else 0.0
            for c, se, df in (('discovery', 0.25, 12.7), ('heca', 0.3, 15.1)):
                if c == 'heca' and i >= n - 10: continue          # ten discovery-only genes
                rows[c].append([g, repr(float(true + rng.normal(0, se))), repr(se), repr(df), 'TRUE'])
        rows['heca'].append(['H_ONLY', '0.1', '0.3', '15.1', 'TRUE']); rows['discovery'].append(['NOFIT', '0.0', '0.3', '12.7', 'FALSE'])
        self.paths = {c: self.root / f'reports/mde_analysis/v1_{c}_mde_per_gene.tsv.gz' for c in rows}
        for c in rows: write_table(self.paths[c], rows[c])
        cfg = json.loads((ROOT / 'config/combined_bound.json').read_text(encoding='utf-8'))
        cfg['synthetic'] = True
        for c in rows:
            cfg['inputs'][c]['sha256'] = hashlib.sha256(self.paths[c].read_bytes()).hexdigest(); cfg['inputs'][c]['expected_rows'] = len(rows[c])
        self.cfg_path = self.root / 'config/combined_bound.json'; self.cfg_path.write_text(json.dumps(cfg, indent=1), encoding='utf-8'); self.cfg = cfg

    def run_py(self, script, *args):
        return subprocess.run([sys.executable, str(self.root / script), '--root', str(self.root), *args], capture_output=True, text=True)

    def test_freeze_then_analysis_recovers_planted_genes_and_binds_outputs(self):
        f = self.run_py('scripts/combined_bound_freeze.py'); self.assertEqual(f.returncode, 0, f.stderr)
        a = self.run_py('scripts/combined_bound_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        summary = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8'))
        self.assertTrue(summary['synthetic']); self.assertEqual(summary['family_size'], 390)
        self.assertEqual(summary['cohort_only'], {'discovery': 10, 'heca': 1})
        found = {g['gene_id'] for g in summary['primary']['significant_genes']}
        self.assertEqual(found, set(self.planted)); self.assertEqual(summary['label'], 'combined_discovery')
        self.assertGreater(summary['scenarios']['alpha_0.05']['mde_log2_quantiles']['q0.5'], 0)
        for name in summary['scenarios']:
            self.assertLess(summary['scenarios'][name]['mde_log2_quantiles']['q0.5'],
                            min(summary['scenarios'][name]['single_cohort_mde_log2_median'].values()), 'combined MDE must be tighter than either cohort')
        receipt = json.loads((self.root / self.cfg['outputs']['results_receipt']).read_text(encoding='utf-8'))
        for rel, h in receipt['output_sha256'].items(): self.assertEqual(h, hashlib.sha256((self.root / rel).read_bytes()).hexdigest())
        self.assertTrue(receipt['synthetic'])
        again = self.run_py('scripts/combined_bound_analysis.py'); self.assertNotEqual(again.returncode, 0, 'the analysis must run once')
        with gzip.open(self.root / self.cfg['outputs']['per_gene'], 'rt', encoding='utf-8') as fh: rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(len(rows), 390); self.assertNotIn('NOFIT', {r['gene_id'] for r in rows})

    def test_null_fixture_gets_no_combined_discovery(self):
        for c in self.paths:
            with gzip.open(self.paths[c], 'rt', encoding='utf-8') as fh: rows = list(csv.reader(fh, delimiter='\t'))
            rng = np.random.default_rng(3)
            for r in rows[1:]:
                if r[0] in self.planted: r[1] = repr(float(rng.normal(0, float(r[2]))))
            write_table(self.paths[c], rows[1:])
            self.cfg['inputs'][c]['sha256'] = hashlib.sha256(self.paths[c].read_bytes()).hexdigest()
        self.cfg_path.write_text(json.dumps(self.cfg, indent=1), encoding='utf-8')
        self.assertEqual(self.run_py('scripts/combined_bound_freeze.py').returncode, 0)
        a = self.run_py('scripts/combined_bound_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        summary = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8'))
        self.assertEqual(summary['label'], 'no_combined_discovery'); self.assertEqual(summary['primary']['genes_BH_q_lt_0_05'], 0)

    def test_freeze_refuses_a_pin_that_does_not_match_the_file(self):
        self.cfg['inputs']['heca']['sha256'] = '0' * 64; self.cfg_path.write_text(json.dumps(self.cfg), encoding='utf-8')
        self.assertNotEqual(self.run_py('scripts/combined_bound_freeze.py').returncode, 0)
        self.assertFalse((self.root / self.cfg['outputs']['freeze_receipt']).exists())

    def test_analysis_refuses_without_freeze_and_after_a_frozen_file_changes(self):
        self.assertNotEqual(self.run_py('scripts/combined_bound_analysis.py').returncode, 0, 'no freeze, no run')
        self.assertEqual(self.run_py('scripts/combined_bound_freeze.py').returncode, 0)
        p = self.root / 'scripts/combined_bound_analysis.py'; p.write_text(p.read_text(encoding='utf-8') + '\n# changed\n', encoding='utf-8')
        r = self.run_py('scripts/combined_bound_analysis.py'); self.assertNotEqual(r.returncode, 0); self.assertIn('changed since the freeze', r.stderr)
        self.assertFalse((self.root / self.cfg['outputs']['results_receipt']).exists())

    def test_synthetic_flag_must_agree_between_freeze_and_configuration(self):
        self.assertEqual(self.run_py('scripts/combined_bound_freeze.py').returncode, 0)
        self.cfg['synthetic'] = False; self.cfg_path.write_text(json.dumps(self.cfg, indent=1), encoding='utf-8')
        r = self.run_py('scripts/combined_bound_analysis.py'); self.assertNotEqual(r.returncode, 0)


if __name__ == '__main__':
    unittest.main()
