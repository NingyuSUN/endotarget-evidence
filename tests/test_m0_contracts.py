import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("m0_contracts", Path(__file__).resolve().parents[1] / "scripts/audit_m0_contracts.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class M0ContractTests(unittest.TestCase):
    def setUp(self):
        Path("tmp").mkdir(exist_ok=True)
        self.tmp = tempfile.TemporaryDirectory(dir="tmp")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def counts(self, text):
        path = self.root / "counts.tsv"
        path.write_text(text)
        return path

    def test_fractional_discovery_preserved_but_not_accepted_as_integer_heca(self):
        path = self.counts("gene_id\ta\tb\nENSG00000000001\t1.25\t2\n")
        summary, _ = m.count_contract(path, ["a", "b"], False)
        self.assertEqual(summary["library_totals"]["a"], "1.25")
        with self.assertRaisesRegex(m.ContractError, "fractional"):
            m.count_contract(path, ["a", "b"], True)

    def test_donor_swap_rejected_even_with_same_donor_set(self):
        path = self.counts("gene_id\tb\ta\nENSG00000000001\t1\t2\n")
        with self.assertRaisesRegex(m.ContractError, "order"):
            m.count_contract(path, ["a", "b"], True)

    def test_nonfinite_negative_and_duplicate_gene_rejected(self):
        for body in ("ENSG00000000001\tNaN\t2\n", "ENSG00000000001\t-1\t2\n",
                     "ENSG00000000001\t1\t2\nENSG00000000001\t3\t4\n"):
            with self.subTest(body=body):
                path = self.counts("gene_id\ta\tb\n" + body)
                with self.assertRaises(m.ContractError):
                    m.count_contract(path, ["a", "b"], True)

    def test_duplicate_person_does_not_satisfy_donor_minimum(self):
        rows = [{"donor_id": "same", "disease": "ENDO", "cells": "100"},
                {"donor_id": "same", "disease": "non-ENDO", "cells": "100"}]
        with self.assertRaisesRegex(m.ContractError, "duplicate biological donor"):
            m.clinical_contract(rows, 1, "cells")

    def test_disease_batch_confounded_design_rejected(self):
        rows = [{"donor_id": str(i), "disease": "ENDO" if i < 5 else "non-ENDO"} for i in range(10)]
        path = self.root / "design.tsv"
        with path.open("w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["donor_id", "intercept", "ENDO_minus_non_ENDO", "batch"])
            for r in rows:
                case = int(r["disease"] == "ENDO")
                w.writerow([r["donor_id"], 1, case, case])
        with self.assertRaisesRegex(m.ContractError, "rank deficiency"):
            m.design_contract(path, rows)

if __name__ == "__main__":
    unittest.main()
