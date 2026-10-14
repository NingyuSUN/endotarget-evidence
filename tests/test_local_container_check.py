"""Boundaries for local synthetic container verification (no daemon required)."""
import importlib.util
import pathlib
import tempfile
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / 'scripts/check_local_container_nextflow.py'


class LocalContainerCheck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('local_check', SCRIPT)
        cls.check = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.check)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def table(self, name, text):
        p = self.root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        return p

    def test_numeric_tolerance_does_not_hide_labels_or_missing_values(self):
        a = self.table('a.tsv', 'gene\tx\tmissing\nENSG1\t0.1\tNA\n')
        b = self.table('b.tsv', 'gene\tx\tmissing\nENSG1\t0.10000000001\tNA\n')
        self.assertEqual(self.check.compare_tsv(a, b)['rows'], 1)
        b.write_text('gene\tx\tmissing\nENSG2\t0.1\tNA\n')
        with self.assertRaisesRegex(ValueError, 'cell'):
            self.check.compare_tsv(a, b)
        b.write_text('gene\tx\tmissing\nENSG1\t0.1\tnan\n')
        with self.assertRaisesRegex(ValueError, 'cell'):
            self.check.compare_tsv(a, b)

    def test_numeric_difference_row_loss_and_extra_columns_fail(self):
        a = self.table('a.tsv', 'gene\tx\nA\t0.1\nB\t2\n')
        for text in ('gene\tx\nA\t0.11\nB\t2\n', 'gene\tx\nA\t0.1\n',
                     'gene\tx\nA\t0.1\textra\nB\t2\n'):
            with self.subTest(text=text):
                b = self.table('b.tsv', text)
                with self.assertRaises(ValueError):
                    self.check.compare_tsv(a, b)

    def test_comparison_requires_identical_complete_table_inventory(self):
        for name in ('a.tsv', 'sub/b.tsv'):
            self.table('native/' + name, 'id\tx\nA\t1\n')
            self.table('new/' + name, 'id\tx\nA\t1\n')
        result = self.check.compare_tree(self.root / 'native', self.root / 'new', 2)
        self.assertEqual(len(result), 2)
        self.table('new/unexpected.tsv', 'id\tx\nA\t1\n')
        with self.assertRaisesRegex(ValueError, 'inventory'):
            self.check.compare_tree(self.root / 'native', self.root / 'new', 2)

    def test_trace_requires_exact_stage_names_and_cache_status(self):
        stages = {'DONOR_DE', 'RANK_TARGETS', 'MATCHED_NULLS'}
        trace = self.table('trace.tsv', 'name\tstatus\texit\n' + ''.join(
            f'{x} (tag)\tCOMPLETED\t0\n' for x in sorted(stages)))
        self.assertEqual(len(self.check.check_trace(trace, stages, False)), 3)
        with self.assertRaisesRegex(ValueError, 'CACHED'):
            self.check.check_trace(trace, stages, True)
        trace.write_text('name\tstatus\texit\nDONOR_DE (tag)\tCACHED\t0\n' * 3)
        with self.assertRaises(ValueError):
            self.check.check_trace(trace, stages, True)

    def test_environment_drops_native_runtime_and_arbitrary_path(self):
        env = self.check.clean_environment({'PATH': '/opt/conda/bin:/bin', 'CONDA_PREFIX': '/opt/conda',
                                          'R_LIBS_USER': '/native/R', 'PYTHONPATH': '/native/python',
                                          'WSL_DISTRO_NAME': 'EndoTargetRuntime'}, self.root, 'run-123')
        self.assertNotIn('CONDA_PREFIX', env)
        self.assertNotIn('R_LIBS_USER', env)
        self.assertNotIn('PYTHONPATH', env)
        self.assertNotIn('conda', env['PATH'])
        self.assertEqual(env['DOCKER_HOST'], 'unix:///run/endotarget-docker.sock')
        self.assertEqual(env['ENDOTARGET_CONTAINER_UID'], '1003')
        self.assertEqual(env['ENDOTARGET_CONTAINER_GID'], '1004')

    def test_rejects_tag_image_path_escape_and_symlink(self):
        self.check.validate_image('127.0.0.1:5000/endo@sha256:' + 'a' * 64)
        for name in ('endo:latest', 'sha256:' + 'a' * 64, 'endo@sha256:' + 'a' * 63):
            with self.assertRaisesRegex(ValueError, 'digest'):
                self.check.validate_image(name)
        self.table('inside.txt', 'x')
        self.assertEqual(self.check.safe_path(self.root, 'inside.txt'), self.root / 'inside.txt')
        for name in ('../outside', '/etc/passwd'):
            with self.assertRaisesRegex(ValueError, 'project'):
                self.check.safe_path(self.root, name)
        (self.root / 'link').symlink_to('/etc/passwd')
        with self.assertRaisesRegex(ValueError, 'project'):
            self.check.safe_path(self.root, 'link')

    def test_registry_image_is_bound_to_completed_build_and_current_sources(self):
        image = 'registry/endo@sha256:' + 'a' * 64
        source = self.table('Dockerfile', 'FROM frozen\n')
        build = {'build_succeeded': True, 'image_id': 'sha256:' + 'b' * 64,
                 'source_hashes': {'Dockerfile': self.check.sha(source)}}
        inspected = {'Id': build['image_id'], 'RepoDigests': [image]}
        self.check.check_image_binding(self.root, image, inspected, build)
        inspected['Id'] = 'sha256:' + 'c' * 64
        with self.assertRaisesRegex(ValueError, 'build'):
            self.check.check_image_binding(self.root, image, inspected, build)
        inspected['Id'] = build['image_id']
        inspected['RepoDigests'] = []
        with self.assertRaisesRegex(ValueError, 'digest'):
            self.check.check_image_binding(self.root, image, inspected, build)
        inspected['RepoDigests'] = [image]
        source.write_text('FROM changed\n')
        with self.assertRaisesRegex(ValueError, 'source'):
            self.check.check_image_binding(self.root, image, inspected, build)

    def test_actual_inspect_requires_isolation_and_both_mounts(self):
        work = '/home/analyst/te/work/local-container/run/nfwork/aa/bb'
        inspect = {'Config': {'User': '1003:1004', 'Image': 'registry/endo@sha256:' + 'a' * 64},
                   'HostConfig': {'ReadonlyRootfs': True, 'NetworkMode': 'none',
                                  'Memory': 2147483648, 'Privileged': False,
                                  'Tmpfs': {'/tmp': 'rw,nosuid,nodev,size=1g,mode=0700,uid=1003,gid=1004'}},
                   'Mounts': [{'Destination': '/home/analyst/te', 'RW': False},
                              {'Destination': work, 'RW': True}]}
        self.check.check_container_inspect(inspect, pathlib.Path('/home/analyst/te'), pathlib.Path(work))
        for key, value in [('NetworkMode', 'bridge'), ('ReadonlyRootfs', False), ('Memory', 0),
                           ('Privileged', True)]:
            original = inspect['HostConfig'][key]
            inspect['HostConfig'][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.check.check_container_inspect(inspect, pathlib.Path('/home/analyst/te'), pathlib.Path(work))
            inspect['HostConfig'][key] = original
        inspect['Mounts'][0]['RW'] = True
        with self.assertRaises(ValueError):
            self.check.check_container_inspect(inspect, pathlib.Path('/home/analyst/te'), pathlib.Path(work))


if __name__ == '__main__':
    unittest.main()
