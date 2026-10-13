"""Hand-solvable synthetic ranking contracts; no external effects or biology claims."""
import csv
import gzip
import importlib.util
import json
import pathlib
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIELDS = ["gene_id", "logFC", "t", "P.Value", "adj.P.Val", "AveExpr"]
def gene(n):
    return f"ENSG{n:011d}"
def row(n, t=1, fc=1):
    return dict(zip(FIELDS, [gene(n), str(fc), str(t), "0.5", "0.8", "4"]))
def edge(a, b, w=1):
    return {"gene1": gene(a), "gene2": gene(b), "weight": str(w)}
def impl(test):
    path = ROOT / "scripts/rank_targets.py"
    if not path.exists():
        test.fail("Ranking implementation missing")
    spec = importlib.util.spec_from_file_location("rank_targets", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def write_table(path, rows, fields):
    with (gzip.open(path, "wt") if str(path).endswith(".gz") else path.open("w")) as handle:
        w = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

class RankingTests(unittest.TestCase):
    def test_failed_network_eligible_gene_cannot_shrink_universe(self):
        m = impl(self)
        all_rows = [row(1), row(2), row(3, "nan")]
        with self.assertRaisesRegex(ValueError, "missing.*fit"):
            m.rank([row(1), row(2)], [edge(1, 2), edge(2, 3)],
                   eligibility_rows=all_rows, test_min_universe=2)

    def test_failed_fit_outside_induced_graph_is_reported_excluded(self):
        m = impl(self)
        all_rows = [row(1), row(2), row(3, "nan")]
        result = m.rank([row(1), row(2)], [edge(1, 2), edge(3, 4)],
                        eligibility_rows=all_rows, test_min_universe=2)
        self.assertEqual(result["universe"], [gene(1), gene(2)])
        self.assertEqual(result["counts"]["discovery_expressed_genes"], 3)
        self.assertEqual(result["excluded"][0]["gene_id"], gene(3))

    def test_complete_expression_eligibility_must_be_unique_and_contain_fits(self):
        m = impl(self)
        for all_rows in [[row(1), row(1)], [row(1)], [dict(row(1), gene_id="ENSG1")], []]:
            with self.subTest(all_rows=all_rows), self.assertRaises(ValueError):
                m.rank([row(1), row(2)], [edge(1, 2)],
                       eligibility_rows=all_rows, test_min_universe=2)

    def test_two_node_closed_form_and_mass(self):
        m = impl(self)
        result = m.rank([row(1, 1, -2), row(2, 0, 0)], [edge(1, 2)], test_min_universe=2)
        scores = {r["gene_id"]: r["score"] for r in result["rankings"]["G1"]}
        self.assertAlmostEqual(scores[gene(1)], 2/3, places=9)
        self.assertAlmostEqual(scores[gene(2)], 1/3, places=9)
        self.assertAlmostEqual(sum(scores.values()), 1, places=12)
        self.assertLess(result["propagation"]["final_l1_change"], 1e-10)
        self.assertLess(result["propagation"]["fixed_point_l1_residual"], 1e-10)

    def test_weighted_path_degree_and_column_normalization(self):
        m = impl(self)
        result = m.rank([row(1, 0), row(2, 1), row(3, 0)], [edge(1, 2, .7), edge(2, 3, 1)], test_min_universe=3)
        b1 = result["rankings"]["B1"]
        self.assertEqual([r["gene_id"] for r in b1], [gene(2), gene(3), gene(1)])
        self.assertAlmostEqual(b1[0]["score"], 1.7)
        scores = {r["gene_id"]: r["score"] for r in result["rankings"]["G1"]}
        self.assertAlmostEqual(scores[gene(2)], 2/3, places=9)
        self.assertAlmostEqual(scores[gene(1)], (1/3)*(.7/1.7), places=9)
        self.assertAlmostEqual(scores[gene(3)], (1/3)*(1/1.7), places=9)

    def test_induced_isolation_not_background_isolation(self):
        m = impl(self)
        result = m.rank([row(1), row(2), row(3), row(5)], [edge(1, 2), edge(3, 4)], test_min_universe=2)
        self.assertEqual(result["universe"], [gene(1), gene(2)])
        self.assertEqual([r["gene_id"] for r in result["excluded"]], [gene(3), gene(5)])
        self.assertEqual({r["reason"] for r in result["excluded"]}, {"zero_degree_in_discovery_induced_graph"})
        for records in result["rankings"].values():
            self.assertEqual({r["gene_id"] for r in records}, {gene(1), gene(2)})
        self.assertEqual(result["counts"]["discovery_expressed_genes"], 4)
        self.assertEqual(result["counts"]["excluded_induced_zero_degree"], 2)

    def test_ties_signs_and_permutation_invariance(self):
        m = impl(self)
        rows = [row(3, -1, 0), row(2, 1, -1), row(1, -1, 1)]
        edges = [edge(1, 2), edge(2, 3), edge(1, 3)]
        result = m.rank(rows, edges, test_min_universe=3)
        other = m.rank(rows[::-1], edges[::-1], test_min_universe=3)
        self.assertEqual(result, other)
        for records in result["rankings"].values():
            self.assertEqual([r["gene_id"] for r in records], [gene(1), gene(2), gene(3)])
            self.assertEqual([r["discovery_sign"] for r in records], [1, -1, 0])
            self.assertEqual([r["logFC"] for r in records], [1, -1, 0])

    def test_disconnected_nonisolated_components_allowed(self):
        m = impl(self)
        result = m.rank([row(1), row(2, 0), row(3, 0), row(4, 0)], [edge(1, 2), edge(3, 4)], test_min_universe=4)
        scores = {r["gene_id"]: r["score"] for r in result["rankings"]["G1"]}
        self.assertAlmostEqual(scores[gene(1)], 2/3, places=9)
        self.assertEqual(scores[gene(3)], 0)
        self.assertEqual(scores[gene(4)], 0)
        self.assertEqual(len(result["universe"]), 4)

    def test_minimum_applies_after_induced_isolation(self):
        m = impl(self)
        rows = [row(i) for i in range(1, 51)]
        with self.assertRaisesRegex(ValueError, "universe.*50"):
            m.rank(rows, [edge(1, 2)])
        with self.assertRaises(ValueError):
            m.rank([row(1), row(2)], [edge(1, 2)], test_min_universe=0)

    def test_zero_seeds_and_nonconvergence_fail(self):
        m = impl(self)
        with self.assertRaisesRegex(ValueError, "seed"):
            m.rank([row(1, 0), row(2, 0)], [edge(1, 2)], test_min_universe=2)
        with self.assertRaisesRegex(ValueError, "converge"):
            m.propagate([1., 0.], [0], [1], [1.], [1., 1.], max_iterations=1)

    def test_extreme_finite_t_values_normalize_without_overflow(self):
        m = impl(self)
        result = m.rank([row(1, 1e308), row(2, -1e308)], [edge(1, 2)], test_min_universe=2)
        self.assertEqual([r["score"] for r in result["rankings"]["G1"]], [.5, .5])

    def test_bad_de_values_duplicates_and_schema_fail(self):
        m = impl(self)
        for field, values in {"t": ["nan", "inf", "bad"], "logFC": ["-inf"], "P.Value": ["-0.1", "1.1"], "adj.P.Val": ["nan", "2"], "AveExpr": ["nan"]}.items():
            for value in values:
                bad = row(1)
                bad[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    m.rank([bad, row(2)], [edge(1, 2)], test_min_universe=2)
        for rows in [[row(1), row(1)], [dict(row(1), gene_id="ENSG1"), row(2)], [{k:v for k,v in row(1).items() if k != "t"}, row(2)], []]:
            with self.assertRaises(ValueError):
                m.rank(rows, [edge(1, 2)], test_min_universe=2)

    def test_bad_reference_including_outside_induction_fails(self):
        m = impl(self)
        rows = [row(1), row(2)]
        bad_edges = [[edge(1, 1)], [edge(1, 2), edge(2, 1)], [edge(3, 4, .69)], [edge(3, 4, "nan")], [edge(3, 4, 1.01)], [{"gene1": gene(1), "gene2": gene(2), "combined_score": "1"}], [dict(edge(1, 2), gene2="SYMBOL")]]
        for edges in bad_edges:
            with self.subTest(edges=edges), self.assertRaises(ValueError):
                m.rank(rows, [edge(1, 2)] + edges, test_min_universe=2)

    def test_real_cli_synthetic_only_outputs_complete_same_universe(self):
        m = impl(self)
        with tempfile.TemporaryDirectory(prefix="ranking-test-", dir=ROOT / "work") as tmp:
            folder = pathlib.Path(tmp)
            de, edges, out = folder/"de.tsv", folder/"edges.tsv.gz", folder/"ranks"
            write_table(de, [row(i, i, -1 if i%2 else 0) for i in range(1, 61)], FIELDS)
            write_table(edges, [edge(i, i+1, .8) for i in range(1, 60)], ["gene1", "gene2", "weight"])
            command = [str(ROOT/"envs/m0-h5/bin/python"), str(ROOT/"scripts/rank_targets.py"), "--de", str(de), "--de-all", str(de), "--edges", str(edges), "--output-dir", str(out), "--synthetic-test"]
            p = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            provenance = json.loads((out/"provenance.json").read_text())
            self.assertEqual(provenance["status"], "SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE")
            self.assertEqual(provenance["counts"]["final_universe_genes"], 60)
            self.assertFalse(provenance["external_data_used"])
            for method in ["B0", "B1", "G1"]:
                with (out/f"{method}.tsv").open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len(rows), 60)
                self.assertEqual([int(r["rank"]) for r in rows], list(range(1, 61)))
                with (out/f"{method}_top20.tsv").open() as handle:
                    self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), rows[:20])
                with (out/f"{method}_top50.tsv").open() as handle:
                    self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), rows[:50])
                self.assertEqual(m.digest(out/f"{method}.tsv"), provenance["outputs"][f"{method}.tsv"]["sha256"])
            self.assertEqual(m.digest(de), provenance["inputs"]["discovery_de"]["sha256"])
            self.assertEqual(m.digest(edges), provenance["inputs"]["reference_edges"]["sha256"])
            self.assertEqual(m.digest(out/"universe.tsv"), provenance["outputs"]["universe.tsv"]["sha256"])
            self.assertTrue((out/"induced_edges.tsv").exists())
            # Preserve existing output instead of silently overwriting a freeze.
            again = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(again.returncode, 0)
            # Production must reject an unpinned reference; no external input exists.
            p = subprocess.run(command[:-1] + ["--output-dir", str(folder/"production")], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("hash", p.stderr)
            # A fitted-only production call cannot silently omit expression-pass failures.
            missing_all = command[:-1].copy()
            j = missing_all.index("--de-all")
            del missing_all[j:j+2]
            p = subprocess.run(missing_all, capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("requires --de-all", p.stderr)
            # Exact-U omission fails through the real CLI with >=50 remaining fits.
            de_all = folder/"de_all.tsv"
            write_table(de_all, [row(i, i) for i in range(1, 61)] + [row(61, "nan")], FIELDS)
            write_table(edges, [edge(i, i+1, .8) for i in range(1, 61)], ["gene1", "gene2", "weight"])
            protected = command.copy()
            protected[protected.index("--de-all")+1] = str(de_all)
            protected[protected.index("--output-dir")+1] = str(folder/"protected")
            p = subprocess.run(protected, capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("missing finite discovery fit", p.stderr)
            self.assertFalse((folder/"protected").exists())

    def test_tsv_headers_duplicate_columns_extra_or_missing_fields_fail(self):
        m = impl(self)
        with tempfile.TemporaryDirectory(prefix="ranking-schema-", dir=ROOT/"work") as tmp:
            p = pathlib.Path(tmp)/"bad.tsv"
            for content in ["gene_id\tt\tt\n", "gene_id\tlogFC\tt\tP.Value\tadj.P.Val\tAveExpr\n"+gene(1)+"\t1\t1\t.1\t.2\t4\textra\n", "gene_id\tlogFC\tt\tP.Value\tadj.P.Val\tAveExpr\n"+gene(1)+"\t1\t1\t.1\t.2\n"]:
                p.write_text(content)
                with self.assertRaises(ValueError):
                    m.read_de(p)

if __name__ == "__main__":
    unittest.main()
