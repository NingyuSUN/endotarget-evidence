import importlib.util, json, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("freeze_analysis",ROOT/"scripts/freeze_analysis.py")
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class FreezeTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=ROOT/"tmp")
        self.base=Path(self.temp.name)
        self.de=self.base/"de";self.rank=self.base/"ranking";self.null=self.base/"nulls"
        for p in [self.de,self.rank,self.null]:p.mkdir()
        self.genes=[f"ENSG{i:011d}" for i in range(1,61)]
        table="gene_id\n"+"\n".join(self.genes)+"\n"
        for name in ["de_all.tsv","de_fitted.tsv"]:(self.de/name).write_text(table)
        (self.rank/"universe.tsv").write_text(table)
        (self.rank/"induced_edges.tsv").write_text("gene1\tgene2\tweight\n")
        for method in ["B0","B1","G1"]:
            (self.rank/f"{method}.tsv").write_text(table)
            for k in [20,50]:
                (self.rank/f"{method}_top{k}.tsv").write_text("gene_id\n"+"\n".join(self.genes[:k])+"\n")
            (self.null/f"{method}_null_membership.tsv").write_text("gene_id\n")
        (self.null/"all_U_directions.tsv").write_text(table)
        self.refresh()
    def tearDown(self):self.temp.cleanup()
    def refresh(self):
        dp={"data_kind":"public_cohort","cohort":"discovery",
            "output_sha256":{str(p):mod.sha(p) for p in self.de.iterdir() if p.name!="provenance.json"}}
        (self.de/"provenance.json").write_text(json.dumps(dp))
        rp={"status":"DISCOVERY_RANKING_ARTIFACT_NOT_COMPLETE_EXTERNAL_FREEZE","external_data_used":False,
            "inputs":{"discovery_de":{"sha256":mod.sha(self.de/"de_fitted.tsv")},
                      "discovery_de_all":{"sha256":mod.sha(self.de/"de_all.tsv")}},
            "outputs":{p.name:{"sha256":mod.sha(p)} for p in self.rank.iterdir() if p.name!="provenance.json"}}
        (self.rank/"provenance.json").write_text(json.dumps(rp))
        np={"status":"DISCOVERY_NULL_MEMBERSHIPS_FROZEN_NOT_COMPLETE_EXTERNAL_FREEZE","external_data_used":False,
            "inputs":{p.name:{"sha256":mod.sha(p)} for p in self.rank.iterdir()},
            "outputs":{p.name:{"sha256":mod.sha(p)} for p in self.null.iterdir() if p.name!="status.json"},
            "methods":{m:{"status":"NOT_ESTIMABLE","sampled_sets":0,"reason":"declared_toy_shortage"} for m in ["B0","B1","G1"]}}
        np["inputs"]["discovery_de"]={"sha256":mod.sha(self.de/"de_fitted.tsv")}
        (self.null/"status.json").write_text(json.dumps(np))
    def valid(self):return mod.validate_discovery(self.de,self.rank,self.null)
    def test_explicit_failed_null_is_not_hidden(self):
        result=self.valid();self.assertEqual(result["universe_size"],60)
        self.assertEqual(set(result["null_status"].values()),{"NOT_ESTIMABLE"})
    def test_top_list_replacement_fails_even_with_updated_hash(self):
        p=self.rank/"G1_top20.tsv";p.write_text(p.read_text().replace(self.genes[0],self.genes[-1]))
        self.refresh()
        with self.assertRaisesRegex(ValueError,"prefix"):self.valid()
    def test_synthetic_cannot_be_scientific_freeze(self):
        p=self.de/"provenance.json";j=json.loads(p.read_text());j["data_kind"]="synthetic_test_only";p.write_text(json.dumps(j))
        with self.assertRaisesRegex(ValueError,"real discovery"):self.valid()
    def test_wrong_null_upstream_fails(self):
        p=self.null/"status.json";j=json.loads(p.read_text());j["inputs"]["B0.tsv"]["sha256"]="0"*64;p.write_text(json.dumps(j))
        with self.assertRaisesRegex(ValueError,"different ranking"):self.valid()
    def test_missing_direction_gene_fails_with_updated_hash(self):
        p=self.null/"all_U_directions.tsv";p.write_text("gene_id\n"+"\n".join(self.genes[:-1])+"\n");self.refresh()
        with self.assertRaisesRegex(ValueError,"prediction directions"):self.valid()
    def test_freeze_hash_detects_changed_file(self):
        file=self.base/"thing.txt";file.write_text("original")
        manifest=self.base/"freeze.json";manifest.write_text(json.dumps({"phase":"toy","artifacts":{mod.rel(file):mod.sha(file)}}))
        self.assertEqual(mod.verify_manifest(manifest,"toy")["phase"],"toy")
        file.write_text("changed")
        with self.assertRaisesRegex(ValueError,"Changed frozen artifact"):mod.verify_manifest(manifest,"toy")
if __name__=="__main__":unittest.main()
