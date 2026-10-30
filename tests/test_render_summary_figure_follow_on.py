# Boundary tests for the follow-on summary figure renderer: it must draw only numbers that
# agree across the published inputs, record them in its receipt, and refuse tampered input.
import hashlib, importlib.util, json, pathlib, shutil, subprocess, sys, tempfile, unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/render_summary_figure_follow_on_v1.py'
spec = importlib.util.spec_from_file_location('follow_on_figure', SCRIPT)
fig = importlib.util.module_from_spec(spec); spec.loader.exec_module(fig)


class FollowOnFigureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'project'
        for rel in fig.INPUTS:
            dst = self.root / rel; dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / rel, dst)

    def run_script(self):
        return subprocess.run([sys.executable, str(SCRIPT), '--root', str(self.root)], capture_output=True, text=True)

    def edit_json(self, rel, change):
        p = self.root / rel; d = json.loads(p.read_text(encoding='utf-8')); change(d)
        p.write_text(json.dumps(d), encoding='utf-8')

    def test_renders_and_receipt_binds_inputs_outputs_and_numbers(self):
        r = self.run_script(); self.assertEqual(r.returncode, 0, r.stderr)
        receipt = json.loads((self.root / fig.RECEIPT).read_text(encoding='utf-8'))
        svg = (self.root / fig.OUTPUT).read_bytes()
        self.assertEqual(receipt['output_sha256'], hashlib.sha256(svg).hexdigest())
        self.assertEqual(receipt['script_sha256'], hashlib.sha256(SCRIPT.read_bytes()).hexdigest())
        self.assertEqual(set(receipt['input_sha256']), set(fig.INPUTS))
        for rel, digest in receipt['input_sha256'].items():
            self.assertEqual(digest, hashlib.sha256((self.root / rel).read_bytes()).hexdigest())
        self.assertIs(receipt['new_scientific_analysis'], False)
        n = receipt['numbers_drawn']
        self.assertEqual(n['disease_genes_q_lt_0_05'], {'discovery': 0, 'heca': 0})
        self.assertEqual(n['v1_B0_top20_transfer']['evaluable'], 20)
        self.assertEqual(n['phase_B0_top20_transfer']['selected'], 20)
        text = svg.decode('utf-8')
        for c in ('discovery', 'heca'):
            self.assertIn(f"{n['phase_genes_q_lt_0_05'][c]:,}", text)
            self.assertIn(f"{n['recall_at_11_per_arm'][c]['median'] * 100:.0f}%", text)
            self.assertEqual(n['theta'][c]['n_donors'], 35 if c == 'discovery' else 42)
        self.assertIn('no_detectable_difference', text)
        # The same inputs render the same bytes.
        again = self.run_script(); self.assertEqual(again.returncode, 0, again.stderr)
        self.assertEqual(hashlib.sha256((self.root / fig.OUTPUT).read_bytes()).hexdigest(), receipt['output_sha256'])

    def assert_rejected(self, message):
        r = self.run_script()
        self.assertNotEqual(r.returncode, 0, message)
        self.assertFalse((self.root / fig.OUTPUT).exists(), 'no figure may be written after a failed check')
        self.assertFalse((self.root / fig.RECEIPT).exists(), 'no receipt may be written after a failed check')

    def test_phase_count_disagreeing_with_the_evaluation_receipt_is_rejected(self):
        self.edit_json('reports/mde_analysis/phase_discovery_summary.json',
                       lambda d: d['observed'].__setitem__('genes_BH_q_lt_0_05', d['observed']['genes_BH_q_lt_0_05'] + 1))
        self.assert_rejected('phase gene count drift accepted')

    def test_synthetic_input_is_rejected(self):
        self.edit_json('reports/secretory_program/heca_summary.json', lambda d: d.__setitem__('synthetic', True))
        self.assert_rejected('synthetic secretory summary accepted')

    def test_changed_label_is_rejected(self):
        self.edit_json('reports/secretory_program/discovery_summary.json', lambda d: d.__setitem__('label', 'attenuation_detected'))
        self.assert_rejected('label not matching the frozen result accepted')

    def test_disease_fit_with_a_discovery_is_rejected(self):
        self.edit_json('reports/mde_analysis/v1_heca_summary.json', lambda d: d['observed'].__setitem__('genes_BH_q_lt_0_05', 1))
        self.assert_rejected('v1 disease discovery accepted')

    def test_missing_frozen_b0_gene_is_rejected(self):
        p = self.root / 'reports/replication_v2_dossier_pool.tsv'
        lines = p.read_text(encoding='utf-8').splitlines(keepends=True)
        header = lines[0].rstrip('\r\n').split('\t'); col = header.index('in_B0_top20')
        drop = next(i for i, l in enumerate(lines[1:], 1) if l.split('\t')[col] == 'TRUE')
        p.write_text(''.join(lines[:drop] + lines[drop + 1:]), encoding='utf-8')
        self.assert_rejected('pool with 19 B0 genes accepted')

    def test_recall_without_valid_draws_is_rejected(self):
        p = self.root / 'reports/mde_analysis/curve_heca_curve_summary.tsv'
        rows = p.read_text(encoding='utf-8').splitlines(keepends=True)
        header = rows[0].rstrip('\r\n').split('\t'); v, m = header.index('valid'), header.index('recall_full_median')
        for i, l in enumerate(rows[1:], 1):
            cells = l.rstrip('\r\n').split('\t')
            if cells[header.index('design')] == 'primary' and cells[v] == '0':
                cells[m] = '0.5'; rows[i] = '\t'.join(cells) + '\n'; break
        p.write_text(''.join(rows), encoding='utf-8')
        self.assert_rejected('recall reported for a level with no valid draw accepted')


if __name__ == '__main__':
    unittest.main()
