"""Bounded migration-verifier rejection tests; no scientific inputs."""
import hashlib
import pathlib
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
from verify_local_migration import IntegrityChecker


class MigrationIntegrityBoundaries(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.tmp.name) / "project"
        self.root.mkdir()
        (self.root / "evidence.txt").write_bytes(b"frozen evidence\n")
        self.digest = hashlib.sha256(b"frozen evidence\n").hexdigest()
        self.checker = IntegrityChecker(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def test_historical_absolute_path_reads_local_copy(self):
        self.checker.check_hash("/home/analyst/te/evidence.txt", self.digest)
        self.assertEqual(self.checker.unique_files_hashed, 1)

    def test_wrong_digest_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            self.checker.check_hash("evidence.txt", "0" * 64)

    def test_relative_and_historical_path_escape_are_rejected(self):
        for name in ("../outside.txt", "/home/analyst/te/../outside.txt"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "escapes"):
                self.checker.resolve(name)

    def test_other_absolute_root_and_symlink_escape_are_rejected(self):
        outside = pathlib.Path(self.tmp.name) / "outside.txt"
        outside.write_bytes(b"outside")
        with self.assertRaisesRegex(ValueError, "absolute path"):
            self.checker.resolve(str(outside))
        (self.root / "outside-link").symlink_to(outside)
        with self.assertRaisesRegex(ValueError, "escapes"):
            self.checker.resolve("outside-link")

    def test_missing_evidence_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Missing"):
            self.checker.check_hash("missing.txt", self.digest)


if __name__ == "__main__":
    unittest.main()
