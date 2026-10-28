#!/usr/bin/env python3
"""Run the synthetic unit tests that a fresh checkout can execute.

Standard library only. Equivalent to ``python -m unittest discover -s tests``
except that the tests listed in SKIPPED are left out by name and reported with
their reason. The delivery tests read a historical README through ``git show``,
so the checkout needs the full published history (fetch-depth 0 in CI). No
check is modified or relaxed, and nothing here reruns a scientific analysis.
"""
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
SKIPPED = {
    "test_rank_targets.RankingTests.test_real_cli_synthetic_only_outputs_complete_same_universe":
        "runs the ranking CLI through the host-specific interpreter envs/m0-h5/bin/python",
}


def flatten(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from flatten(item)
        else:
            yield item


def wanted(test):
    test_id = test.id()
    return not any(test_id == key or test_id.startswith(key + ".") for key in SKIPPED)


def main():
    # Same import environment as `python -m unittest discover -s tests` run from
    # the repository root: the root itself (for `from scripts import ...`) and the
    # tests directory (added by discover) are both importable.
    sys.path.insert(0, str(ROOT))
    discovered = unittest.TestLoader().discover(str(TESTS), pattern="test_*.py")
    all_tests = list(flatten(discovered))
    kept = [test for test in all_tests if wanted(test)]
    for key, reason in SKIPPED.items():
        print(f"SKIPPED {key}: {reason}", flush=True)
    print(f"collected {len(all_tests)} tests, running {len(kept)}", flush=True)
    result = unittest.TextTestRunner(verbosity=1).run(unittest.TestSuite(kept))
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
