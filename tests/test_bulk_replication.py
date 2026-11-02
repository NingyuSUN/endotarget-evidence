# Synthetic boundary checks for the bulk replication: GEO matrix and annotation parsing, probe
# collapsing, the E1 gate and interaction, E2 transfer counting, E3 moderation, labels, and
# the freeze and data-receipt refusals. Fixtures are invented; no real expression data.
import csv, gzip, hashlib, importlib.util, json, pathlib, shutil, subprocess, sys, tempfile, unittest
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
FILES = ['scripts/bulk_replication_analysis.py', 'scripts/bulk_replication_fetch.py', 'scripts/bulk_replication_pin.py', 'scripts/followon_freeze.py',
         'scripts/bulk_replication_metadata_audit.py', 'tests/test_bulk_replication.py', 'docs/BULK_REPLICATION_PROTOCOL.md']
spec = importlib.util.spec_from_file_location('bulk', ROOT / 'scripts/bulk_replication_analysis.py'); bulk = importlib.util.module_from_spec(spec); spec.loader.exec_module(bulk)
sys.path.insert(0, str(ROOT / 'scripts')); from followon_freeze import sha  # noqa: E402


def wt(path, text): path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text, encoding='utf-8')


class BulkFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'project'
        for rel in FILES:
            (self.root / rel).parent.mkdir(parents=True, exist_ok=True)
            if (ROOT / rel).exists(): shutil.copyfile(ROOT / rel, self.root / rel)
        self.cfg = json.loads((ROOT / 'config/bulk_replication.json').read_text(encoding='utf-8')); self.cfg['synthetic'] = True
        self.cfg['parameters']['matched_draws'] = 300
        self.cfg['cohorts'] = {'GSE51981': self.cfg['cohorts']['GSE51981'], 'GSE141549': self.cfg['cohorts']['GSE141549']}
        self.cfg['inputs'] = {k: v for k, v in self.cfg['inputs'].items() if k not in ('samples_GSE6364',)}
        self.cfg['downloads']['files'] = {'GSE51981_series_matrix.txt.gz': 'file://x', 'GPL570.annot.gz': 'file://y'}
        self.cfg['downloads']['hgnc'] = {'path': 'data/hgnc.tsv', 'sha256_lf': ''}

    def build(self, theta_true, b0_agree_fraction=1.0, seed=3):
        rng = np.random.default_rng(seed); n_genes = 600; genes = [f'ENSG{i:011d}' for i in range(n_genes)]
        # HGNC: symbol SYMi, entrez i+1
        wt(self.root / 'data/hgnc.tsv', 'hgnc_id\tsymbol\tstatus\tentrez_id\tensembl_gene_id\n' + ''.join(f'HGNC:{i}\tSYM{i}\tApproved\t{i + 1}\t{g}\n' for i, g in enumerate(genes)))
        self.cfg['downloads']['hgnc']['sha256_lf'] = hashlib.sha256((self.root / 'data/hgnc.tsv').read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        # program: first 200 genes with beta; secretory shift 1 unit
        prog = {g: float(rng.choice([-1, 1]) * rng.uniform(0.5, 2.0)) for g in genes[:200]}
        for key in ('phase_program_heca', 'phase_program_discovery'):
            p = self.root / self.cfg['inputs'][key]['path']; p.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(p, 'wt', encoding='utf-8', newline='') as f:
                w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 'logFC', 'adj.P.Val', 'fit_estimable'])
                for g in genes: w.writerow([g, repr(prog.get(g, 0.0)), '0.001' if g in prog else '0.9', 'TRUE'])
        # samples: PE 20 non-ENDO (10 no pathology, 10 other), PE 20 ENDO, SEC 15 non-ENDO, SEC 20 ENDO
        design = [('proliferative', 'non-ENDO', 'No Uterine Pelvic Pathology')] * 10 + [('proliferative', 'non-ENDO', 'Uterine Pelvic Pathology')] * 10 + [('proliferative', 'ENDO', '')] * 20 + \
                 [('early_secretory', 'non-ENDO', 'No Uterine Pelvic Pathology')] * 8 + [('mid_secretory', 'non-ENDO', 'Uterine Pelvic Pathology')] * 7 + [('early_secretory', 'ENDO', '')] * 10 + [('mid_secretory', 'ENDO', '')] * 10
        gsms = [f'GSM{i:07d}' for i in range(len(design))]; tissue = {'proliferative': 'Proliferative Endometrial tissue', 'early_secretory': 'Early Secretory Endometrial tissue', 'mid_secretory': 'Mid-Secretory Endometrial tissue'}
        rows = []
        for gsm, (ph, dz, sub) in zip(gsms, design):
            ch = f'tissue: {tissue[ph]} | endometriosis/no endometriosis: {"Endometriosis" if dz == "ENDO" else "Non-Endometriosis"} | ' + (f'presence or absence of uterine/pelvic pathology: {sub}' if sub else 'endometriosis severity: Moderate/Severe')
            rows.append({'gsm': gsm, 'title': gsm, 'source': 'x', 'platform': 'GPL570', 'characteristics': ch, 'phase_guess': '', 'disease_guess': '', 'tissue_guess': '', 'stage_text': '', 'hormone_mention': ''})
        for key in ('samples_GSE51981', 'samples_GSE141549'):
            p = self.root / self.cfg['inputs'][key]['path']; p.parent.mkdir(parents=True, exist_ok=True)
            with p.open('w', encoding='utf-8', newline='') as f: w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
        wt(self.root / self.cfg['inputs']['metadata_audit']['path'], json.dumps({'series': {}}))
        # expression: two probes per gene (second one dimmer), program genes shift by beta in secretory, interaction theta_true
        X = rng.normal(7, 1, (n_genes, len(design)))
        for j, (ph, dz, _) in enumerate(design):
            sec = ph != 'proliferative'
            for i, g in enumerate(genes[:200]): X[i, j] += prog[g] * (1.0 * sec + theta_true * sec * (dz == 'ENDO'))
        # B0 top-20: genes 300-319 with planted direction in proliferative ENDO
        b0 = {genes[300 + i]: int(rng.choice([-1, 1])) for i in range(20)}
        for g, sgn in b0.items():
            i = genes.index(g); agree = rng.random() < b0_agree_fraction
            for j, (ph, dz, _) in enumerate(design):
                if ph == 'proliferative' and dz == 'ENDO': X[i, j] += (sgn if agree else -sgn) * 1.5
        with (self.root / 'reports/replication_v2_dossier_pool.tsv').open('w', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 'B0_rank', 'in_B0_top20', 'discovery_sign'])
            for k, (g, s) in enumerate(b0.items(), 1): w.writerow([g, k, 'TRUE', s])
        (self.root / 'reports/combined_bound').mkdir(parents=True, exist_ok=True)
        with gzip.open(self.root / 'reports/combined_bound/combined_per_gene.tsv.gz', 'wt', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n'); w.writerow(['gene_id', 't']); [w.writerow([g, repr(float(rng.normal()))]) for g in genes]
        ddir = self.root / self.cfg['downloads']['directory']; ddir.mkdir(parents=True, exist_ok=True)
        with gzip.open(ddir / 'GSE51981_series_matrix.txt.gz', 'wt', encoding='utf-8') as f:
            f.write('!Series_title\t"synthetic"\n!series_matrix_table_begin\n"ID_REF"\t' + '\t'.join(f'"{g}"' for g in gsms) + '\n')
            for i, g in enumerate(genes):
                f.write(f'"{i}_at"\t' + '\t'.join(f'{v:.4f}' for v in X[i]) + '\n'); f.write(f'"{i}_dim_at"\t' + '\t'.join(f'{v - 3:.4f}' for v in X[i]) + '\n')
            f.write('!series_matrix_table_end\n')
        with gzip.open(ddir / 'GPL570.annot.gz', 'wt', encoding='utf-8') as f:
            f.write('#ID = probe\n^ANNOTATION\nID\tGene title\tGene symbol\tGene ID\n')
            for i, g in enumerate(genes): f.write(f'{i}_at\tt\tSYM{i}\t{i + 1}\n{i}_dim_at\tt\tSYM{i}\t{i + 1}\n')
        rec = {'files': {n: {'sha256': hashlib.sha256((ddir / n).read_bytes()).hexdigest()} for n in ('GSE51981_series_matrix.txt.gz', 'GPL570.annot.gz')}}
        wt(self.root / self.cfg['downloads']['data_receipt'], json.dumps(rec))
        for name, s in self.cfg['inputs'].items(): s['sha256'] = sha(self.root / s['path'])
        wt(self.root / 'config/bulk_replication.json', json.dumps(self.cfg, indent=1))
        return prog, b0

    def run_py(self, script, *args): return subprocess.run([sys.executable, str(self.root / script), '--root', str(self.root), *args], capture_output=True, text=True)
    def freeze(self): return self.run_py('scripts/followon_freeze.py', '--config', 'config/bulk_replication.json')

    def test_planted_attenuation_and_full_transfer(self):
        self.build(theta_true=-0.4, b0_agree_fraction=1.0)
        self.assertEqual(self.freeze().returncode, 0)
        a = self.run_py('scripts/bulk_replication_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        s = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8'))
        r = s['series']['GSE51981']; self.assertEqual(r['n_genes'], 600); self.assertEqual(r['groups'], {'proliferative_ENDO': 20, 'proliferative_non_ENDO': 20, 'secretory_ENDO': 20, 'secretory_non_ENDO': 15})
        e1 = r['E1']['heca_program']; self.assertTrue(e1['gate']['pass']); self.assertEqual(e1['label'], 'attenuation_detected'); self.assertLess(e1['M1_theta']['estimate'], -0.2)
        e2 = r['E2']['all_non_ENDO']; self.assertEqual(e2['evaluable'], 20); self.assertEqual(e2['same_direction'], 20); self.assertEqual(e2['label'], 'transfers')
        self.assertIn('other_pathology_controls_only', r['E2']); self.assertGreaterEqual(r['E3']['all_non_ENDO']['genes_BH_q_lt_0_05'], 15)
        self.assertEqual(s['series']['GSE141549']['status'], 'not_analyzed'); self.assertEqual(s['E1_across_series']['label'], 'not_replicated')
        self.assertTrue((self.root / self.cfg['outputs']['results_receipt']).exists())
        self.assertNotEqual(self.run_py('scripts/bulk_replication_analysis.py').returncode, 0, 'runs once')

    def test_null_fixture_random_directions(self):
        self.build(theta_true=0.0, b0_agree_fraction=0.5, seed=11)
        self.assertEqual(self.freeze().returncode, 0)
        a = self.run_py('scripts/bulk_replication_analysis.py'); self.assertEqual(a.returncode, 0, a.stderr)
        s = json.loads((self.root / self.cfg['outputs']['summary']).read_text(encoding='utf-8')); r = s['series']['GSE51981']
        self.assertEqual(r['E1']['heca_program']['label'], 'no_detectable_difference'); self.assertGreater(r['E1']['heca_program']['M1_theta']['mde'], 0)
        self.assertEqual(r['E2']['all_non_ENDO']['label'], 'does_not_transfer')

    def test_data_receipt_mismatch_and_freeze_are_enforced(self):
        self.build(theta_true=0.0)
        self.assertNotEqual(self.run_py('scripts/bulk_replication_analysis.py').returncode, 0, 'no freeze, no run')
        self.assertEqual(self.freeze().returncode, 0)
        p = self.root / self.cfg['downloads']['directory'] / 'GPL570.annot.gz'; p.write_bytes(p.read_bytes() + b'\n')
        r = self.run_py('scripts/bulk_replication_analysis.py'); self.assertNotEqual(r.returncode, 0); self.assertIn('differs from the data receipt', r.stderr)

    def test_parsers_and_moderation(self):
        self.build(theta_true=0.0)
        gsms, probes, vals = bulk.read_series_matrix(self.root / 'data/bulk_replication/GSE51981_series_matrix.txt.gz')
        self.assertEqual(len(gsms), 75); self.assertEqual(len(probes), 1200)
        annot = bulk.read_annot(self.root / 'data/bulk_replication/GPL570.annot.gz'); self.assertEqual(len(annot), 1200)
        genes, M = bulk.probes_to_genes(probes, vals, annot, bulk.hgnc_maps(self.root / 'data/hgnc.tsv'))
        self.assertEqual(len(genes), 600); self.assertTrue(np.all(M.mean(axis=1) > vals.mean()))   # brighter probe kept
        # trigamma inverse: trigamma(2) = pi^2/6 - 1
        self.assertAlmostEqual(bulk.trigamma_inverse(np.pi ** 2 / 6 - 1), 2.0, places=5)
        self.assertAlmostEqual(bulk.trigamma_inverse(0.1), 10.42, delta=0.1)   # trigamma(10.42) ~ 0.1
        # homogeneous true variances: chi-square sampling noise only -> large prior df, strong shrinkage
        s2 = np.random.default_rng(0).chisquare(10, 500) / 10 * 0.5; mod, d0, s0 = bulk.squeeze_var(s2, 10)
        self.assertGreater(d0, 20); self.assertLess(np.var(mod), np.var(s2)); self.assertAlmostEqual(s0, 0.5, delta=0.1)
        # strongly heterogeneous true variances (log-normal spread) -> small finite prior df, weak shrinkage
        rng = np.random.default_rng(1); true = np.exp(rng.normal(0, 1.5, 2000)); s2h = true * rng.chisquare(10, 2000) / 10
        modh, d0h, _ = bulk.squeeze_var(s2h, 10)
        self.assertLess(d0h, 5); self.assertGreater(np.corrcoef(modh, s2h)[0, 1], 0.95)


if __name__ == '__main__':
    unittest.main()
