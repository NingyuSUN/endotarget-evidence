"""Focused delivery rejection tests: real temporary files, no cohort execution."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'scripts'))
try:
    import verify_local_delivery as delivery
except ModuleNotFoundError as error:
    if error.name != 'verify_local_delivery':
        raise
    delivery = None

SOURCE = '4cc7580851cbfcb776385ac0ae858fe5f81e818b'
OLD_README = 'd73153d84255856c794640814c9c43a4b522c88b64e0aab8099dabea4588273e'
ARCHIVE = 'docs/archive/README_srv01_4cc7580.md'


class LocalDeliveryBoundaries(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(delivery, 'Delivery verifier feature has not been implemented')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / 'project'
        self.root.mkdir()
        for folder in ('docs/archive', 'reports', 'scripts', 'receipts', 'work'):
            (self.root / folder).mkdir(parents=True, exist_ok=True)
        self.original = subprocess.check_output(['git', '-C', str(REPO), 'show', SOURCE + ':README.md'])
        self.assertEqual(hashlib.sha256(self.original).hexdigest(), OLD_README)
        (self.root / ARCHIVE).write_bytes(self.original)
        (self.root / 'README.md').write_text('Current local delivery\n')
        self.checker = delivery.DeliveryChecker(self.root)

    def test_only_exact_old_readme_hash_is_redirected(self):
        self.checker.check_hash('README.md', OLD_README)
        self.assertIn((self.root / ARCHIVE).resolve(), self.checker.cache)
        self.assertNotIn((self.root / 'README.md').resolve(), self.checker.cache)

    def test_unexpected_expected_hash_rejected_even_when_current_file_matches(self):
        current = hashlib.sha256((self.root / 'README.md').read_bytes()).hexdigest()
        with self.assertRaisesRegex(ValueError, 'Unexpected historical README hash'):
            self.checker.check_hash('README.md', current)

    def test_tampered_archive_has_no_fallback_to_matching_current_readme(self):
        (self.root / ARCHIVE).write_bytes(b'tampered archive')
        (self.root / 'README.md').write_bytes(self.original)
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            self.checker.check_hash('README.md', OLD_README)

    def test_explicit_base_and_absolute_path_do_not_use_archive_exception(self):
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            self.checker.check_hash('README.md', OLD_README, base='.')
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            self.checker.check_hash(str(self.root / 'README.md'), OLD_README)

    def test_current_file_inventory_hashes_current_readme_without_redirection(self):
        result = delivery.verify_current_files(self.checker, ['README.md'])
        self.assertEqual(result['sha256']['README.md'], hashlib.sha256(b'Current local delivery\n').hexdigest())
        self.assertNotEqual(result['sha256']['README.md'], OLD_README)

    def test_final_inventory_cannot_omit_method_card_and_current_execution_docs(self):
        self.assertTrue(hasattr(delivery, 'verify_required_current_files'), 'Core delivery file gate is not implemented')
        with self.assertRaisesRegex(ValueError, 'Missing core current delivery files'):
            delivery.verify_required_current_files(self.checker, ['README.md'])

    def test_explicit_paths_reject_traversal_noncanonical_case_scope_and_self_reference(self):
        (self.root / 'reports/ok.md').write_text('ok')
        for name in ('../outside.md', './reports/ok.md', 'reports//ok.md', '/tmp/a',
                     'reports\\ok.md', 'Reports/ok.md', 'work/secret',
                     '.git/config', 'receipts/local_delivery_verification.json'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                delivery.current_path(self.checker, name)

    def test_symlink_outside_project_rejected(self):
        outside = Path(self.tmp.name) / 'outside.md'
        outside.write_text('external')
        (self.root / 'reports/link.md').symlink_to(outside)
        with self.assertRaisesRegex(ValueError, 'symlink|escapes'):
            delivery.current_path(self.checker, 'reports/link.md')

    def test_current_file_list_rejects_duplicates(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            delivery.verify_current_files(self.checker, ['README.md', 'README.md'])

    def test_only_two_byte_archives_are_exempt_from_historical_broken_links(self):
        (self.root / 'docs/archive/other.md').write_text('[bad](missing.tsv)')
        with self.assertRaisesRegex(ValueError, 'Missing local Markdown target'):
            delivery.verify_markdown(self.checker, ['docs/archive/other.md'])
        (self.root / ARCHIVE).write_text('[historical](missing.tsv)')
        self.assertEqual(delivery.verify_markdown(self.checker, [ARCHIVE])['explicit_archives_skipped'], [ARCHIVE])

    def test_active_markdown_checks_references_decoded_paths_and_ignores_code(self):
        (self.root / 'reports/real file.tsv').write_text('data')
        (self.root / 'docs/active.md').write_text('[ok](../reports/real%20file.tsv#row) [web](https://example.com) [mail](mailto:a@example.com) [here](#header)\n[x][ref]\n[ref]: ../reports/real%20file.tsv\n```\n[not a link](missing)\n```\n')
        result = delivery.verify_markdown(self.checker, ['docs/active.md'])
        self.assertEqual(result['local_targets_checked'], 2)
        (self.root / 'docs/active.md').write_text('[missing](../reports/no.tsv)')
        with self.assertRaisesRegex(ValueError, 'Missing local Markdown target'):
            delivery.verify_markdown(self.checker, ['docs/active.md'])

    def test_shortcut_reference_links_are_not_silently_ignored(self):
        (self.root / 'docs/active.md').write_text('[report]\n\n[report]: missing.tsv\n')
        with self.assertRaisesRegex(ValueError, 'Missing local Markdown target'):
            delivery.verify_markdown(self.checker, ['docs/active.md'])

    def test_first_run_markdown_can_link_exact_atomic_receipt_output(self):
        (self.root / 'docs/active.md').write_text('[delivery](../receipts/local_delivery_verification.json)')
        self.assertFalse((self.root / 'receipts/local_delivery_verification.json').exists())
        result = delivery.verify_markdown(self.checker, ['docs/active.md'])
        self.assertEqual(result['planned_atomic_output_links'], [{
            'document': 'docs/active.md',
            'target': 'receipts/local_delivery_verification.json',
        }])

    def test_first_run_output_exception_does_not_allow_another_missing_receipt(self):
        (self.root / 'docs/active.md').write_text('[other](../receipts/another_delivery_verification.json)')
        with self.assertRaisesRegex(ValueError, 'Missing local Markdown target'):
            delivery.verify_markdown(self.checker, ['docs/active.md'])

    def test_markdown_cannot_escape_through_symlink_or_parent(self):
        (Path(self.tmp.name) / 'outside.md').write_text('x')
        (self.root / 'docs/active.md').write_text('[outside](../../outside.md)')
        with self.assertRaisesRegex(ValueError, 'escapes'):
            delivery.verify_markdown(self.checker, ['docs/active.md'])

    def test_visual_review_binds_all_required_artifacts(self):
        import json
        names = ['reports/executive_brief_en.pdf', 'reports/figures/donor_flow_v1.svg',
                 'reports/figures/donor_composition_qc_v1.svg', 'reports/figures/matched_null_sensitivity_v1.svg',
                 'reports/figures/evidence_matrix_v1.svg']
        artifacts = {}
        for name in names:
            path = self.root / name
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(b'%PDF-1.4' if name.endswith('.pdf') else b'<svg/>')
            artifacts[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        obj = {'status': 'PASS', 'reviewer': 'Codex visual inspection', 'pdf_pages': 2, 'artifacts': artifacts}
        path = self.root / 'receipts/visual.json'
        path.write_text(json.dumps(obj))
        result = delivery.verify_visual_review(self.checker, 'receipts/visual.json')
        self.assertFalse(result['automated_visual_proof'])
        (self.root / names[-1]).write_bytes(b'<svg>changed</svg>')
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch|changed during'):
            delivery.verify_visual_review(self.checker, 'receipts/visual.json')

    def test_visual_review_missing_matrix_is_not_a_pass(self):
        import json
        (self.root / 'receipts/visual.json').write_text(json.dumps({
            'status': 'PASS', 'reviewer': 'Codex visual inspection', 'pdf_pages': 2,
            'artifacts': {'reports/executive_brief_en.pdf': '0' * 64}}))
        with self.assertRaisesRegex(ValueError, 'coverage'):
            delivery.verify_visual_review(self.checker, 'receipts/visual.json')

    def test_pending_runtime_receipts_cannot_pass_from_absent_other_fields(self):
        import json
        for name, function in [('local_container_nextflow', delivery.verify_nextflow),
                               ('local_primary_container_replay', delivery.verify_primary)]:
            (self.root / 'receipts' / (name + '.json')).write_text(json.dumps({'status': 'RUNNING'}))
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'has not passed'):
                function(self.checker, {})

    def test_stale_pdf_source_hash_is_rejected(self):
        import json
        md = self.root / 'reports/executive_brief_en.md'
        pdf = self.root / 'reports/executive_brief_en.pdf'
        script = self.root / 'scripts/render_executive_brief_pdf.py'
        md.write_text('original brief'); pdf.write_bytes(b'%PDF-1.4\n/Type /Page\n/Type /Page\n');script.write_text('# renderer')
        obj = {'source': 'reports/executive_brief_en.md', 'output': 'reports/executive_brief_en.pdf',
               'source_sha256': hashlib.sha256(md.read_bytes()).hexdigest(),
               'output_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
               'script_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
               'pages': 2, 'new_scientific_analysis': False}
        (self.root / 'receipts/executive_brief_pdf.json').write_text(json.dumps(obj))
        md.write_text('new claims added')
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            delivery.verify_pdf(self.checker)

    def test_pdf_receipt_two_pages_cannot_hide_actual_third_page(self):
        import json
        md = self.root / 'reports/executive_brief_en.md'
        pdf = self.root / 'reports/executive_brief_en.pdf'
        script = self.root / 'scripts/render_executive_brief_pdf.py'
        md.write_text('brief'); pdf.write_bytes(b'%PDF-1.4\n' + b'/Type /Page\n' * 3);script.write_text('# renderer')
        obj = {'source': 'reports/executive_brief_en.md', 'output': 'reports/executive_brief_en.pdf',
               'source_sha256': hashlib.sha256(md.read_bytes()).hexdigest(),
               'output_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
               'script_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
               'pages': 2, 'new_scientific_analysis': False}
        (self.root / 'receipts/executive_brief_pdf.json').write_text(json.dumps(obj))
        with self.assertRaisesRegex(ValueError, 'page objects'):
            delivery.verify_pdf(self.checker)


class LocalBuildRuntimeBoundaries(unittest.TestCase):
    def setUp(self):
        self.assertTrue(hasattr(delivery, 'verify_os_git_runtime'), 'OS/Git runtime gate is not implemented')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        folder = self.root / 'containers/statistics-local'
        folder.mkdir(parents=True)
        lock = (REPO / 'containers/statistics-local/os-packages.lock.tsv').read_text()
        (folder / 'os-packages.lock.tsv').write_text(lock)
        expected = dict(line.split('\t') for line in lock.splitlines())
        self.assertEqual(len(expected), 37)
        self.checker = delivery.DeliveryChecker(self.root)
        self.build = {
            'required_os_packages': expected,
            'provenance_git_available': True,
            'runtime_inspection': {
                'git_version': 'git version 2.47.3\n', 'git_head': SOURCE,
                'os_packages': ''.join(name + ':amd64\t' + version + '\n' for name, version in expected.items())}}

    def test_all_37_locked_packages_and_actual_git_are_accepted(self):
        result = delivery.verify_os_git_runtime(self.checker, self.build)
        self.assertEqual(result['installed_os_packages'], 37)
        self.assertEqual(result['git_head_at_build_inspection'], SOURCE)

    def test_missing_git_flag_or_actual_version_cannot_pass(self):
        for field in ('provenance_git_available', 'git_version'):
            import copy
            obj = copy.deepcopy(self.build)
            if field == 'provenance_git_available': obj[field] = False
            else: obj['runtime_inspection'].pop(field)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'Git'):
                delivery.verify_os_git_runtime(self.checker, obj)

    def test_actual_git_version_must_match_locked_git_package(self):
        self.build['runtime_inspection']['git_version'] = 'git version 1.0.0\n'
        with self.assertRaisesRegex(ValueError, 'Git version'):
            delivery.verify_os_git_runtime(self.checker, self.build)

    def test_wrong_declared_os_lock_version_is_rejected(self):
        self.build['required_os_packages']['git'] = '9.9.9'
        with self.assertRaisesRegex(ValueError, 'OS package'):
            delivery.verify_os_git_runtime(self.checker, self.build)

    def test_wrong_measured_os_package_version_is_rejected(self):
        self.build['runtime_inspection']['os_packages'] = self.build['runtime_inspection']['os_packages'].replace('git:amd64\t1:2.47.3-0+deb13u1', 'git:amd64\t9.9.9')
        with self.assertRaisesRegex(ValueError, 'OS package'):
            delivery.verify_os_git_runtime(self.checker, self.build)

    def test_duplicate_normalized_dpkg_package_is_rejected(self):
        self.build['runtime_inspection']['os_packages'] += 'git\t1:2.47.3-0+deb13u1\n'
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            delivery.verify_os_git_runtime(self.checker, self.build)

    def test_invalid_build_git_head_is_rejected(self):
        self.build['runtime_inspection']['git_head'] = 'not-a-commit'
        with self.assertRaisesRegex(ValueError, 'Git'):
            delivery.verify_os_git_runtime(self.checker, self.build)


class LocalImageIdentityBoundaries(unittest.TestCase):
    def setUp(self):
        self.assertTrue(hasattr(delivery, 'verify_image_configuration'), 'Separate image default-user gate is not implemented')
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.checker = delivery.DeliveryChecker(Path(self.tmp.name))
        self.config = {'architecture': 'amd64', 'os': 'linux', 'config': {'User': '57439:57439'}}
        self.build = {'architecture': 'amd64', 'os': 'linux', 'image_config': {'User': '57439:57439'}}
        self.recipe = 'FROM locked-base\nUSER root\nRUN install-something\nUSER 57439:57439\n'

    def test_default_user_and_task_user_can_differ_with_both_verified(self):
        delivery.verify_image_configuration(self.config, self.build, self.recipe)
        taskdirs = ['/home/analyst/te/work/local-test/task' + str(i) for i in range(5)]
        items = []
        for work in taskdirs:
            items.append({'inspect': {
                'Config': {'User': '1003:1004', 'Image': 'locked-image'}, 'Image': 'verified-id',
                'HostConfig': {'NetworkMode': 'none', 'ReadonlyRootfs': True, 'Memory': 2147483648,
                               'Privileged': False, 'SecurityOpt': [],
                               'Tmpfs': {'/tmp': 'rw,nosuid,nodev,uid=1003,gid=1004'}},
                'Mounts': [{'Destination': '/home/analyst/te', 'RW': False},
                           {'Destination': work, 'RW': True}]}})
        receipt = {'image': 'locked-image', 'image_id': 'verified-id', 'container_inspections': items,
                   'container_observer_errors': []}
        delivery.verify_containers(self.checker, receipt, taskdirs)

    def test_default_root_is_rejected_even_if_record_and_recipe_agree(self):
        self.config['config']['User'] = 'root'
        self.build['image_config']['User'] = 'root'
        with self.assertRaisesRegex(ValueError, 'default user'):
            delivery.verify_image_configuration(self.config, self.build, 'FROM base\nUSER root\n')


if __name__ == '__main__':
    unittest.main()
