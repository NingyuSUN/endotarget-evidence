# Synthetic boundary checks for the stromal-state composition test: a planted shift is found,
# a null fixture is not, labels and partition coverage behave, and the freeze is enforced.
import csv, json, pathlib, shutil, subprocess, sys, tempfile, unittest
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
FILES = ['scripts/composition_test_analysis.py', 'scripts/followon_freeze.py', 'tests/test_composition_test.py', 'docs/COMPOSITION_TEST_PROTOCOL.md']
STATES = {'discovery': ['dS1-myofibroblast', 'dS2', 'dS3', 'dS4', 'eF1', 'eF2', 'eF3', 'eF4-CXCL14', 'fib C7', 'fib C7 secretory', 'fib C7 secretory_late', 'fib C7-SFRP2'],
          'heca': ['dStromal_early', 'dStromal_mid', 'dStromal_late', 'eStromal', 'eStromal_MMPs', 'eStromal_cycling', 'Fibroblast_basalis']}


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


class CompositionFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'project'
        for rel in FILES:
            (self.root / rel).parent.mkdir(parents=True, exist_ok=True)
            if (ROOT / rel).exists(): shutil.copyfile(ROOT / rel, self.root / rel)
        self.cfg = json.loads((ROOT / 'config/composition_test.json').read_text(encoding='utf-8')); self.cfg['synthetic'] = True

    def build(self, shift_log_odds):
        rng = np.random.default_rng(11); donors_rows, denom_rows = [], []
        for c in ('discovery', 'heca'):
            ids = [f'{c[0]}{i:02d}' for i in range(22)]; disease = np.array([1] * 12 + [0] * 10) if c == 'discovery' else np.array([1] * 15 + [0] * 7)
            age = rng.normal(0, 6, 22); batch = np.array([i % 3 for i in range(22)]); adeno = np.array([1, 1, 0, 0, 1] + [0] * 17)
            cols = ['intercept', 'ENDO_minus_non_ENDO', 'age_centered_years', 'batch:b1_vs_b0', 'batch:b2_vs_b0', 'adenomyosis']
            design = [{'donor_id': d, 'intercept': 1.0, 'ENDO_minus_non_ENDO': float(disease[i]), 'age_centered_years': float(age[i]),
                       'batch:b1_vs_b0': float(batch[i] == 1), 'batch:b2_vs_b0': float(batch[i] == 2), 'adenomyosis': float(adeno[i])} for i, d in enumerate(ids)]
            write_tsv(self.root / f'reports/final_clinical_design_{c}_primary.tsv', design, ['donor_id'] + cols)
            states = STATES[c]; dec = self.cfg['partition'][c]['decidualized']
            for i, d in enumerate(ids):
                n = int(rng.integers(800, 6000)); logit = -1.0 + shift_log_odds * disease[i] + rng.normal(0, 0.25)
                p_dec = 1 / (1 + np.exp(-logit)); k_dec = rng.binomial(n, p_dec)
                w = rng.dirichlet(np.ones(len(dec))) * k_dec; w_rest = rng.dirichlet(np.ones(len(states) - len(dec))) * (n - k_dec)
                alloc = {}; j = k = 0
                for s in states:
                    if s in dec: alloc[s] = int(round(w[j])); j += 1
                    else: alloc[s] = int(round(w_rest[k])); k += 1
                diff = n - sum(alloc.values()); alloc[states[0]] += diff
                for s in states:
                    donors_rows.append({'cohort': c, 'donor_id': d, 'selected_sample_id': d, 'disease': 'ENDO' if disease[i] else 'non-ENDO', 'source_label_field': 'x', 'source_subtype_label': s,
                                        'cells': alloc[s], 'primary_cells': n, 'proportion_primary': alloc[s] / n, 'availability': 'observed'})
                denom_rows.append({'cohort': c, 'donor_id': d, 'selected_sample_id': d, 'disease': 'ENDO' if disease[i] else 'non-ENDO', 'primary_cells': n})
        write_tsv(self.root / 'reports/stromal_composition_donors.tsv', donors_rows, list(donors_rows[0].keys()))
        write_tsv(self.root / 'reports/stromal_composition_denominators.tsv', denom_rows, list(denom_rows[0].keys()))
        sys.path.insert(0, str(ROOT / 'scripts')); from followon_freeze import sha
        for key in self.cfg['inputs']: self.cfg['inputs'][key]['sha256'] = sha(self.root / self.cfg['inputs'][key]['path'])
        (self.root / 'config').mkdir(exist_ok=True)
        (self.root / 'config/composition_test.json').write_text(json.dumps(self.cfg, indent=1), encoding='utf-8')

    def run_py(self, script, *args):
        return subprocess.run([sys.executable, str(self.root / script), '--root', str(self.root), *args], capture_output=True, text=True)

    def freeze(self): return self.run_py('scripts/followon_freeze.py', '--config', 'config/composition_test.json')

    def test_planted_shift_is_detected_in_both_cohorts_and_replicated(self):
        self.build(shift_log_odds=1.2)
        self.assertEqual(self.freeze().returncode, 0)
        a = self.run_py('scripts/composition_test_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        s = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8'))
        for c in ('discovery', 'heca'):
            self.assertEqual(s['cohorts'][c]['label'], 'composition_shift_detected'); self.assertGreater(s['cohorts'][c]['primary']['estimate'], 0)
            self.assertEqual(s['cohorts'][c]['primary']['df'], 22 - 6)
            self.assertEqual(len(s['cohorts'][c]['per_state_clr']), len(STATES[c]))
        self.assertEqual(s['pair_label'], 'replicated_shift'); self.assertTrue(s['synthetic'])
        receipt = json.loads((self.root / self.cfg['outputs']['results_receipt']).read_text(encoding='utf-8'))
        self.assertEqual(receipt['labels']['pair'], 'replicated_shift')
        self.assertNotEqual(self.run_py('scripts/composition_test_analysis.py').returncode, 0, 'runs once')

    def test_null_fixture_gets_no_detectable_shift_and_a_bound(self):
        self.build(shift_log_odds=0.0)
        self.assertEqual(self.freeze().returncode, 0)
        a = self.run_py('scripts/composition_test_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        s = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8'))
        for c in ('discovery', 'heca'):
            self.assertEqual(s['cohorts'][c]['label'], 'no_detectable_shift')
            self.assertGreater(s['cohorts'][c]['primary']['mde_odds_ratio'], 1.0)
        self.assertEqual(s['pair_label'], 'no_replicated_shift')

    def test_unknown_author_label_is_refused(self):
        self.build(shift_log_odds=0.0)
        p = self.root / 'reports/stromal_composition_donors.tsv'; text = p.read_text(encoding='utf-8').replace('\teStromal\t', '\teStromal_new\t', 1); p.write_text(text, encoding='utf-8')
        sys.path.insert(0, str(ROOT / 'scripts')); from followon_freeze import sha
        self.cfg['inputs']['composition_donors']['sha256'] = sha(p); (self.root / 'config/composition_test.json').write_text(json.dumps(self.cfg), encoding='utf-8')
        self.assertEqual(self.freeze().returncode, 0)
        r = self.run_py('scripts/composition_test_analysis.py'); self.assertNotEqual(r.returncode, 0); self.assertIn('not covered by the frozen partition', r.stderr)

    def test_analysis_refuses_without_freeze_and_after_code_change(self):
        self.build(shift_log_odds=0.0)
        self.assertNotEqual(self.run_py('scripts/composition_test_analysis.py').returncode, 0)
        self.assertEqual(self.freeze().returncode, 0)
        p = self.root / 'scripts/composition_test_analysis.py'; p.write_text(p.read_text(encoding='utf-8') + '\n# changed\n', encoding='utf-8')
        r = self.run_py('scripts/composition_test_analysis.py'); self.assertNotEqual(r.returncode, 0); self.assertIn('changed since the freeze', r.stderr)


if __name__ == '__main__':
    unittest.main()
