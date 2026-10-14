import pathlib
import sys
import tempfile
import unittest

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import check_local_primary_replay as replay


class PrimaryReplayBoundaries(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def table(self, name, content):
        path = self.root / name
        path.write_text(content)
        return path

    def compare(self, old, new):
        return replay.compare_tsv(self.table('old.tsv', old), self.table('new.tsv', new))

    def test_continuous_tolerance_is_absolute_not_relative(self):
        self.compare('gene_id\tlogFC\nENSG1\t1\n', 'gene_id\tlogFC\nENSG1\t1.000000005\n')
        with self.assertRaises(ValueError):
            self.compare('gene_id\tlogFC\nENSG1\t1000000000\n', 'gene_id\tlogFC\nENSG1\t1000000001\n')

    def test_numeric_identifiers_ranks_counts_directions_and_flags_are_exact(self):
        for field in ['gene_id', 'rank', 'selected', 'direction', 'expression_pass', 'draw', 'donor_id', 'raw_library_size', 'filtered_library_size']:
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.compare(f'{field}\tlogFC\n1\t0.5\n', f'{field}\tlogFC\n1.000000001\t0.5\n')

    def test_nonfinite_numeric_changes_and_late_rows_are_checked(self):
        with self.assertRaises(ValueError):
            self.compare('gene_id\tlogFC\nA\t1\nB\t2\n', 'gene_id\tlogFC\nA\t1\nB\tNaN\n')
        with self.assertRaises(ValueError):
            self.compare('gene_id\tlogFC\nA\t1\nB\t2\n', 'gene_id\tlogFC\nB\t2\nA\t1\n')

    def test_unknown_numeric_columns_do_not_get_tolerance(self):
        with self.assertRaises(ValueError):
            self.compare('future_count\n10\n', 'future_count\n10.000000001\n')

    def test_semantic_json_counts_and_booleans_are_exact(self):
        replay.compare_json_values({'effect': 1.0}, {'effect': 1.000000005})
        for old, new in [({'selected':20},{'selected':20.000000001}), ({'flag':True},{'flag':1})]:
            with self.subTest(old=old), self.assertRaises(ValueError):
                replay.compare_json_values(old, new)

    def test_project_relative_paths_cannot_bypass_root(self):
        outside = self.root.parent / (self.root.name + '_outside')
        outside.mkdir()
        self.addCleanup(outside.rmdir)
        (self.root / 'link').symlink_to(outside, target_is_directory=True)
        for name in ['/tmp/source', '../source', 'link/source', 'a/../source']:
            with self.subTest(name=name), self.assertRaises(ValueError):
                replay.project_path(self.root, name)

    def test_provenance_cannot_reference_an_outside_source_with_correct_basename(self):
        with self.assertRaises(ValueError):
            replay.generated_input_path(self.root, '/tmp/de_fitted.tsv')
        with self.assertRaises(ValueError):
            replay.generated_input_path(self.root, '/home/analyst/te/../evil/de_fitted.tsv')

    def test_output_inventory_rejects_missing_extra_or_symlinked_sources(self):
        (self.root / 'de').mkdir()
        (self.root / 'de/de_all.tsv').write_text('gene_id\nA\n')
        with self.assertRaises(ValueError):
            replay.inventory(self.root, {'de/de_all.tsv', 'de/de_fitted.tsv'})
        with self.assertRaises(ValueError):
            replay.inventory(self.root, set())
        outside = self.table('outside.tsv', 'gene_id\nA\n')
        (self.root / 'de/de_fitted.tsv').symlink_to(outside)
        with self.assertRaises(ValueError):
            replay.inventory(self.root / 'de', {'de_all.tsv', 'de_fitted.tsv'})


if __name__ == '__main__':
    unittest.main()
