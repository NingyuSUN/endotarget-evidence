# Synthetic boundary checks for the published-signature test: a planted concordant signature is
# recovered by T1, T2 and T3, a random signature is not labelled as anything, a tiny set stays
# unlabelled, and the freeze / once-only refusals hold. Fixtures are invented; no real data.
import csv, gzip, hashlib, importlib.util, json, pathlib, shutil, sys, tempfile, unittest
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
FILES = ['scripts/published_signature_test_analysis.py', 'scripts/published_signature_build.py', 'scripts/followon_freeze.py', 'tests/test_published_signature_test.py']
spec = importlib.util.spec_from_file_location('pst', ROOT / 'scripts/published_signature_test_analysis.py'); pst = importlib.util.module_from_spec(spec); spec.loader.exec_module(pst)
sys.path.insert(0, str(ROOT / 'scripts')); import followon_freeze  # noqa: E402


def wt(path, text): path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text, encoding='utf-8', newline='\n')


def tsv(path, header, rows, gz=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    f = gzip.open(path, 'wt', encoding='utf-8', newline='') if gz else path.open('w', encoding='utf-8', newline='')
    with f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(header); w.writerows(rows)


class SignatureFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'project'
        for rel in FILES:
            (self.root / rel).parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT / rel, self.root / rel)
        self.cfg = json.loads((ROOT / 'config/published_signature_test.json').read_text(encoding='utf-8')); self.cfg['synthetic'] = True
        self.cfg['parameters']['matched_draws'] = 300
        self.cfg['set_order'] = ['planted', 'random', 'tiny']
        self.cfg['freeze']['files'] = FILES
        self.cfg['freeze']['attestations'] = ['synthetic fixture']
        self.cfg['host_inputs']['discovery_counts']['path'] = 'data/discovery_counts.tsv.gz'; self.cfg['host_inputs']['heca_counts']['path'] = 'data/heca_counts.tsv.gz'

    def build(self, planted_shift=1.0, planted_t=4.0, seed=5):
        rng = np.random.default_rng(seed); n = 1500; genes = [f'ENSG{i:011d}' for i in range(n)]
        planted = genes[:60]; random_set = genes[100:160]; tiny = genes[200:205]
        dirs = {g: int(rng.choice([-1, 1])) for g in planted + random_set + tiny}
        # pooled per-gene table: planted genes carry t aligned with their direction
        t = rng.normal(0, 1, n); se = np.full(n, 0.2)
        for i, g in enumerate(planted): t[i] = dirs[g] * (planted_t + rng.normal(0, 0.5))
        beta = t * se; rows = []
        for i, g in enumerate(genes):
            rows.append([g, f'{beta[i] * 1.1:.6f}', f'{se[i] * 1.4:.6f}', f'{beta[i] * 0.9:.6f}', f'{se[i] * 1.4:.6f}', f'{beta[i]:.6f}', f'{se[i]:.6f}', '20', f'{t[i]:.6f}', '0.5', '0.9', f'{beta[i] - 0.4:.6f}', f'{beta[i] + 0.4:.6f}', 'TRUE' if rng.random() < 0.55 else 'FALSE'])
        tsv(self.root / self.cfg['inputs']['combined_per_gene']['path'], ['gene_id', 'logFC_discovery', 'SE_discovery', 'logFC_heca', 'SE_heca', 'beta_ivw', 'se_ivw', 'df_satterthwaite', 't', 'p', 'q_bh', 'ci_low', 'ci_high', 'same_sign'], rows, gz=True)
        tsv(self.root / self.cfg['inputs']['discovery_mde_table']['path'], ['gene_id', 'expression_bin'], [[g, str(1 + i % 5)] for i, g in enumerate(genes)], gz=True)
        srows = [['planted', g, f'S{g[-4:]}', str(dirs[g]), f'{dirs[g] * 1.5:.3f}', '0.001', 'approved'] for g in planted]
        srows += [['random', g, f'S{g[-4:]}', str(dirs[g]), f'{dirs[g] * 1.5:.3f}', '0.001', 'approved'] for g in random_set]
        srows += [['tiny', g, f'S{g[-4:]}', str(dirs[g]), f'{dirs[g] * 1.5:.3f}', '0.001', 'approved'] for g in tiny]
        tsv(self.root / self.cfg['inputs']['signature_sets']['path'], ['set', 'ensembl_gene_id', 'symbol', 'direction', 'published_logFC_disease_minus_control', 'published_adj_p', 'mapping'], srows)
        # two cohorts: 22 donors each, 11 ENDO; counts with the planted shift in ENDO donors
        for cname, prefix in (('discovery', 'D'), ('heca', 'H')):
            donors = [f'{prefix}{k:02d}' for k in range(22)]; endo = np.array([1.0] * 11 + [0.0] * 11); age = rng.normal(0, 5, 22); batch = (np.arange(22) % 2).astype(float); aden = np.zeros(22); aden[3] = 1
            c = self.cfg['cohorts'][cname]
            tsv(self.root / self.cfg['inputs'][c['design']]['path'], ['donor_id', 'intercept', 'ENDO_minus_non_ENDO', 'age_centered_years', 'batch:b_vs_a', 'adenomyosis'], [[d, '1.0', str(endo[k]), f'{age[k]:.3f}', str(batch[k]), str(aden[k])] for k, d in enumerate(donors)])
            lib = rng.uniform(5e6, 2e7, 22)
            tsv(self.root / self.cfg['inputs'][c['normalization']]['path'], ['donor_id', 'raw_library_size', 'filtered_library_size', 'norm_factor', 'effective_library_size'], [[d, f'{lib[k]:.1f}', f'{lib[k]:.1f}', '1.0', f'{lib[k]:.1f}'] for k, d in enumerate(donors)])
            tsv(self.root / self.cfg['inputs'][c['expression_filter']]['path'], ['gene_id', 'expression_pass'], [[g, 'TRUE' if i < 1400 else 'FALSE'] for i, g in enumerate(genes)])
            base = rng.uniform(3, 9, n); M = np.empty((n, 22))
            for k in range(22):
                lc = base + rng.normal(0, 0.3, n)
                for i, g in enumerate(planted): lc[i] += dirs[g] * planted_shift * endo[k]
                M[:, k] = np.round((2 ** lc) * lib[k] / 1e6)
            tsv(self.root / self.cfg['host_inputs'][c['counts']]['path'], ['gene_id'] + donors, [[g] + [str(int(M[i, k])) for k in range(22)] for i, g in enumerate(genes)], gz=True)
            self.cfg['host_inputs'][c['counts']]['sha256'] = hashlib.sha256((self.root / self.cfg['host_inputs'][c['counts']]['path']).read_bytes()).hexdigest()
        for spec_ in self.cfg['inputs'].values(): spec_['sha256'] = followon_freeze.sha(self.root / spec_['path'])
        wt(self.root / 'config/published_signature_test.json', json.dumps(self.cfg, indent=1))
        wt(self.root / 'docs/PUBLISHED_SIGNATURE_PROTOCOL.md', '# synthetic protocol\n')
        followon_freeze.main(self.root, 'config/published_signature_test.json')

    def test_planted_signature_recovered_and_random_not(self):
        self.build()
        s = pst.main(self.root)
        p, r, t = s['sets']['planted'], s['sets']['random'], s['sets']['tiny']
        self.assertEqual(p['label_T1'], 'enriched'); self.assertEqual(p['label_T2'], 'concordant')
        self.assertGreater(p['T2_direction']['sign_agreement_share'], 0.95); self.assertGreaterEqual(p['T2_direction']['matched_percentile'], 99)
        for c in ('discovery', 'heca'):
            f = p['T3_donor_score'][c]['primary_sign_weighted']
            self.assertEqual(f['status'], 'ESTIMATED'); self.assertAlmostEqual(f['estimate'], 1.0, delta=0.25); self.assertLess(f['p'], 0.001)
            self.assertEqual(p['T3_donor_score'][c]['label'], 'difference_detected_toward_published_direction'); self.assertGreater(f['mde'], 0)
            self.assertEqual(p['T3_donor_score'][c]['genes_in_cohort_family'], 60)
        self.assertEqual(p['label_T3_pair'], 'replicated')
        self.assertEqual(r['label_T1'], 'no_enrichment'); self.assertEqual(r['label_T2'], 'no_detectable_concordance'); self.assertEqual(r['label_T3_pair'], 'not_replicated')
        for c in ('discovery', 'heca'): self.assertEqual(r['T3_donor_score'][c]['label'], 'no_detectable_difference')
        self.assertEqual(t['label_T1'], 'not_labelled'); self.assertEqual(t['label_T2'], 'not_labelled'); self.assertEqual(t['label_T3_pair'], 'not_labelled')
        # outputs, receipt and the once-only rule
        out = self.cfg['outputs']
        for rel in (out['set_results'], out['set_members'], out['donor_scores'], out['summary'], out['results_receipt']): self.assertTrue((self.root / rel).exists(), rel)
        rec = json.loads((self.root / out['results_receipt']).read_text(encoding='utf-8'))
        self.assertEqual(rec['labels']['planted']['T2'], 'concordant'); self.assertIn('host_input_sha256', rec)
        with (self.root / out['donor_scores']).open(encoding='utf-8', newline='') as fh: donor_rows = list(csv.DictReader(fh, delimiter='\t'))
        self.assertEqual(len(donor_rows), 2 * 22 * 2)   # planted and random sets, two cohorts, 22 donors
        with self.assertRaises(SystemExit): pst.main(self.root)

    def test_direction_flip_gives_discordant(self):
        self.build(planted_shift=0.0)
        # flip the published directions of the planted set after the build: t stays aligned with the old direction
        p = self.root / self.cfg['inputs']['signature_sets']['path']
        with p.open(encoding='utf-8', newline='') as fh: rows = list(csv.DictReader(fh, delimiter='\t'))
        for r in rows:
            if r['set'] == 'planted': r['direction'] = str(-int(r['direction'])); r['published_logFC_disease_minus_control'] = str(-float(r['published_logFC_disease_minus_control']))
        with p.open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
        with self.assertRaises(SystemExit): pst.main(self.root)   # pin mismatch: the frozen set table changed
        self.cfg['inputs']['signature_sets']['sha256'] = followon_freeze.sha(p)
        wt(self.root / 'config/published_signature_test.json', json.dumps(self.cfg, indent=1))
        (self.root / self.cfg['outputs']['freeze_receipt']).unlink()
        followon_freeze.main(self.root, 'config/published_signature_test.json')
        s = pst.main(self.root)
        self.assertEqual(s['sets']['planted']['label_T2'], 'discordant'); self.assertEqual(s['sets']['planted']['label_T1'], 'enriched')
        self.assertEqual(s['sets']['planted']['T3_donor_score']['discovery']['label'], 'no_detectable_difference')

    def test_mde_and_matched_draws(self):
        self.assertAlmostEqual(pst.mde_exact(1.0, 30, 0.05, 0.8), 2.90, delta=0.05)   # ncp ~ 2.90 for df 30
        bins = np.array([1, 1, 2, 2, 3, 3, 3]); members = np.array([0, 2, 4, 5]); rng = np.random.default_rng(0)
        idx, order = pst.matched_draws(bins, members, 50, rng)
        self.assertEqual(idx.shape, (50, 4)); self.assertTrue((np.sort(bins[idx], axis=1) == np.sort(bins[members])).all()); self.assertEqual(sorted(order.tolist()), members.tolist())


if __name__ == '__main__':
    unittest.main()
