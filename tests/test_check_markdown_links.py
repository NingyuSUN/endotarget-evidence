# Boundary tests for the Markdown link checker: tracked files only, relative links only,
# historical archive skipped, broken links named and failing the run.
import pathlib, subprocess, sys, tempfile, unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/check_markdown_links.py'
GIT_ID = ['-c', 'user.name=test', '-c', 'user.email=test@example.invalid']


class MarkdownLinkTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name) / 'repo'; self.root.mkdir()
        subprocess.run(['git', 'init', '-q'], cwd=self.root, check=True)
        (self.root / 'docs' / 'archive').mkdir(parents=True)
        (self.root / 'docs' / 'guide.md').write_text('# Guide\n\nBack to the [index](../README.md#top).\n', encoding='utf-8')
        (self.root / 'README.md').write_text(
            '# Top\n\nSee the [guide](docs/guide.md), the [site](https://example.invalid/x), '
            '[mail](mailto:a@example.invalid) and [this section](#top).\n', encoding='utf-8')
        (self.root / 'docs' / 'archive' / 'old.md').write_text('Historical [link](../../gone.md).\n', encoding='utf-8')
        self.commit()

    def commit(self):
        subprocess.run(['git', 'add', '-A'], cwd=self.root, check=True)
        subprocess.run(['git', *GIT_ID, 'commit', '-q', '-m', 'test'], cwd=self.root, check=True)

    def run_check(self):
        return subprocess.run([sys.executable, str(SCRIPT), '--root', str(self.root)], capture_output=True, text=True)

    def test_clean_tree_passes_and_archive_external_and_anchor_links_are_skipped(self):
        r = self.run_check()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn('broken: 0', r.stdout)
        self.assertIn('checked 2 relative links', r.stdout)

    def test_broken_tracked_link_fails_and_is_named(self):
        (self.root / 'README.md').write_text('# Top\n\n[missing](docs/missing.md) and [ok](docs/guide.md)\n', encoding='utf-8')
        self.commit()
        r = self.run_check()
        self.assertEqual(r.returncode, 1)
        self.assertIn('broken: 1', r.stdout)
        self.assertIn('README.md -> docs/missing.md', r.stdout)

    def test_percent_encoded_target_is_decoded_before_lookup(self):
        (self.root / 'docs' / 'two words.md').write_text('x\n', encoding='utf-8')
        (self.root / 'README.md').write_text('[spaced](docs/two%20words.md)\n', encoding='utf-8')
        self.commit()
        r = self.run_check(); self.assertEqual(r.returncode, 0, r.stdout)

    def test_untracked_file_with_broken_link_is_not_checked(self):
        (self.root / 'scratch.md').write_text('[gone](nowhere.md)\n', encoding='utf-8')
        r = self.run_check(); self.assertEqual(r.returncode, 0, r.stdout)
        subprocess.run(['git', 'add', 'scratch.md'], cwd=self.root, check=True)
        r = self.run_check(); self.assertEqual(r.returncode, 1, 'staged file must be checked')


if __name__ == '__main__':
    unittest.main()
