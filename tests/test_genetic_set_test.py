# Synthetic boundary checks for the genetic set-level test: a planted enriched set is labelled
# enriched, a random set is not, small sets are not labelled, and the freeze is enforced.
import csv, gzip, json, pathlib, shutil, subprocess, sys, tempfile, unittest
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
FILES = ['scripts/genetic_set_test_analysis.py', 'scripts/followon_freeze.py', 'tests/test_genetic_set_test.py', 'docs/GENETIC_SET_TEST_PROTOCOL.md']


class GeneticSetFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'project'
        for rel in FILES:
            (self.root / rel).parent.mkdir(parents=True, exist_ok=True)
            if (ROOT / rel).exists(): shutil.copyfile(ROOT / rel, self.root / rel)
        self.cfg = json.loads((ROOT / 'config/genetic_set_test.json').read_text(encoding='utf-8')); self.cfg['synthetic'] = True
        self.cfg['parameters']['matched_draws'] = 500
        rng = np.random.default_rng(5); n = 2000
        genes = [f'G{i:05d}' for i in range(n)]; bins = rng.integers(1, 6, n)
        t = rng.standard_t(20, n); self.enriched = genes[:60]; t[:60] += np.sign(rng.standard_normal(60)) * 2.5
        self.random_set = [genes[i] for i in rng.choice(np.arange(100, n), 60, replace=False)]
        per_gene = self.root / self.cfg['inputs']['combined_per_gene']['path']; per_gene.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(per_gene, 'wt', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(['gene_id', 't', 'logFC_discovery', 'SE_discovery', 'logFC_heca', 'SE_heca', 'beta_ivw', 'ci_low', 'ci_high', 'q_bh', 'same_sign'])
            for i, g in enumerate(genes):
                b1, b2 = t[i] * 0.2 + rng.normal(0, 0.05), t[i] * 0.2 + rng.normal(0, 0.05)
                w.writerow([g, repr(float(t[i])), repr(float(b1)), '0.2', repr(float(b2)), '0.2', repr(float(t[i] * 0.14)), '0', '0', '0.5', 'TRUE' if np.sign(b1) == np.sign(b2) else 'FALSE'])
        mde = self.root / self.cfg['inputs']['discovery_mde_table']['path']; mde.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(mde, 'wt', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 'expression_bin'])
            for g, b in zip(genes, bins): w.writerow([g, int(b)])
        sets = self.root / self.cfg['inputs']['gene_sets']['path']; sets.parent.mkdir(parents=True, exist_ok=True)
        with sets.open('w', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['set', 'ensembl_gene_id', 'symbol', 'mapping'])
            for g in self.enriched: w.writerow(['rahmioglu_nearest', g, g, 'approved'])
            for g in self.random_set: w.writerow(['opentargets_genetic_any', g, g, 'ensembl_from_source'])
            for g in genes[200:205]: w.writerow(['rahmioglu_smr', g, g, 'approved'])
        sys.path.insert(0, str(ROOT / 'scripts')); from followon_freeze import sha
        for key in self.cfg['inputs']: self.cfg['inputs'][key]['sha256'] = sha(self.root / self.cfg['inputs'][key]['path'])
        (self.root / 'config').mkdir(exist_ok=True)
        (self.root / 'config/genetic_set_test.json').write_text(json.dumps(self.cfg, indent=1), encoding='utf-8')

    def run_py(self, script, *args):
        return subprocess.run([sys.executable, str(self.root / script), '--root', str(self.root), *args], capture_output=True, text=True)

    def test_planted_enrichment_random_set_and_small_set(self):
        self.assertEqual(self.run_py('scripts/followon_freeze.py', '--config', 'config/genetic_set_test.json').returncode, 0)
        a = self.run_py('scripts/genetic_set_test_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        s = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8'))
        self.assertEqual(s['sets']['rahmioglu_nearest']['label'], 'enriched')
        self.assertLess(s['sets']['rahmioglu_nearest']['tests']['pooled']['wilcoxon_p_one_sided'], 1e-6)
        self.assertEqual(s['sets']['opentargets_genetic_any']['label'], 'no_enrichment')
        self.assertEqual(s['sets']['rahmioglu_smr']['label'], 'not_labelled')
        self.assertEqual(s['sets']['opentargets_genetic_0.5']['set_size_in_family'], 0)
        self.assertEqual(s['sets']['rahmioglu_nearest']['set_size_in_family'], 60)
        with (self.root / self.cfg['outputs']['set_members']).open(encoding='utf-8') as f: rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(sum(r['set'] == 'rahmioglu_nearest' for r in rows), 60)
        self.assertNotEqual(self.run_py('scripts/genetic_set_test_analysis.py').returncode, 0, 'runs once')

    def test_freeze_required_and_code_change_refused(self):
        self.assertNotEqual(self.run_py('scripts/genetic_set_test_analysis.py').returncode, 0)
        self.assertEqual(self.run_py('scripts/followon_freeze.py', '--config', 'config/genetic_set_test.json').returncode, 0)
        p = self.root / 'scripts/genetic_set_test_analysis.py'; p.write_text(p.read_text(encoding='utf-8') + '\n# changed\n', encoding='utf-8')
        r = self.run_py('scripts/genetic_set_test_analysis.py'); self.assertNotEqual(r.returncode, 0); self.assertIn('changed since the freeze', r.stderr)


if __name__ == '__main__':
    unittest.main()
