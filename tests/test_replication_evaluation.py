"""Synthetic golden cases; never import or read real external effects."""
import importlib.util
import math
import pathlib
import unittest

SPEC = importlib.util.spec_from_file_location("replication", pathlib.Path(__file__).resolve().parents[1] / "scripts/evaluate_replication.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

def row(gene, sign=1, p=.5, q=.5, ci=(-1, 2), valid=True):
    return {"gene_id": gene, "logFC": str(sign), "t": str(2*sign), "P.Value": str(p),
            "adj.P.Val": str(q), "CI_low": str(ci[0]), "CI_high": str(ci[1]),
            "fit_estimable": "TRUE" if valid else "FALSE"}

class GoldenEvaluation(unittest.TestCase):
    def test_exhaustive_reasons_and_denominators(self):
        # 20 selected: 10 agree, 5 opposite, then one of each five exclusion reasons.
        selected = [{"gene_id": f"g{i}", "discovery_sign": 1, "rank": i+1} for i in range(20)]
        selected[18]["discovery_sign"] = 0
        external = {f"g{i}": row(f"g{i}", 1 if i<10 else -1) for i in range(15)}
        external["g17"] = row("g17", valid=False)
        external["g18"] = row("g18")
        external["g19"] = row("g19", 0)
        filters = {f"g{i}": True for i in range(17,20)}
        filters.update({f"g{i}": True for i in range(15)})
        filters["g16"] = False
        details, summary = M.evaluate_set(selected, filters, external, 20)
        self.assertEqual((summary["agree"], summary["evaluable"], summary["selected"]), (10,15,20))
        self.assertEqual(summary["missing_bounds"], [.5,.75])
        self.assertAlmostEqual(summary["concordance"], 2/3)
        self.assertEqual([r["missing_reason"] for r in details[-5:]], [
            "absent_id", "expression_filter_fail", "invalid_fit", "discovery_zero_sign", "external_zero_sign"])
        self.assertEqual(summary["opposite"], 5)
        self.assertFalse(summary["coverage_adequate"])
        self.assertEqual(details[10]["signed_external_t"], -2)
        self.assertEqual(len(details), 20)  # No refill.
    def test_zero_and_tiny_nonzero(self):
        self.assertEqual(M.classify("g",0,{"g":True},{"g":row("g",0)}), "discovery_zero_sign")
        self.assertEqual(M.classify("g",1,{"g":True},{"g":row("g",0)}), "external_zero_sign")
        self.assertEqual(M.classify("g",1,{"g":True},{"g":row("g",1e-300)}), "")
    def test_external_zero_preserves_signed_t_but_not_concordance(self):
        records,summary=M.evaluate_set([{ "gene_id":"g", "discovery_sign":1, "rank":1}],
            {"g":True},{"g":row("g",0)},1)
        self.assertEqual(records[0]["signed_external_t"],0.)
        self.assertEqual(records[0]["missing_reason"],"external_zero_sign")
        self.assertFalse(records[0]["evaluable"])
        self.assertEqual((summary["agree"],summary["opposite"],summary["evaluable"]),(0,0,0))
        self.assertIsNone(summary["concordance"])
    def test_reason_precedence(self):
        bad=row("g",valid=False)
        self.assertEqual(M.classify("g",0,{},{"g":bad}), "absent_id")
        self.assertEqual(M.classify("g",0,{"g":False},{"g":bad}), "expression_filter_fail")
        self.assertEqual(M.classify("g",0,{"g":True},{"g":bad}), "invalid_fit")
    def test_empty_evaluable_bounds(self):
        _, summary=M.evaluate_set([{"gene_id":"gone","rank":1,"discovery_sign":1}],{}, {}, 1)
        self.assertIsNone(summary["concordance"])
        self.assertEqual(summary["missing_bounds"], [0.,1.])
    def test_paired_missing_bounds_prevent_superiority(self):
        b={"missing_bounds":[.5,.7],"concordance":.625,"evaluable":16}
        g={"missing_bounds":[.65,.85],"concordance":.8125,"evaluable":16}
        c=M.primary_contrast(b,g)
        self.assertAlmostEqual(c["missing_bounds"][0],-.05)
        self.assertAlmostEqual(c["missing_bounds"][1],.35)
        self.assertFalse(c["eligible_for_bootstrap_superiority_assessment"])
        self.assertFalse(c["superiority_claim"])
    def test_complete_family_BH_and_legitimate_failed_gene(self):
        # Four family members, only two finite fits: q=(.004,.4,NA,NA), not (.002,.2).
        family=[row("a",p=.001,q=.004),row("b",p=.2,q=.4),
                row("c",p="NA",q="NA",valid=False),row("d",p="NA",q="NA",valid=False)]
        self.assertEqual(M.bh_complete_family(family), {"a":.004,"b":.4,"c":None,"d":None})
        M.validate_bh(family)
        family[0]["adj.P.Val"]=".002"
        with self.assertRaisesRegex(ValueError,"complete testing family"):
            M.validate_bh(family)
    def test_failed_gene_q_zero_rejected(self):
        with self.assertRaisesRegex(ValueError,"failed fit"):
            M.validate_bh([row("bad",p="NA",q=0,valid=False)])
    def test_null_coverage_descriptive_only(self):
        summaries=[{"draw":1,"concordance":.4,"coverage_adequate":True},
                   {"draw":2,"concordance":.8,"coverage_adequate":True},
                   {"draw":3,"concordance":1.,"coverage_adequate":False}]
        ans=M.null_distribution(summaries,.6)
        self.assertEqual((ans["total_sets"],ans["adequate_sets"],ans["inadequate_sets"]),(3,2,1))
        self.assertEqual(ans["observed_descriptive_percentile"],50)
        self.assertIsNone(ans["p_value"])
        self.assertIsNone(M.null_distribution([summaries[2]],.6)["observed_descriptive_percentile"])
    def test_null_own_discovery_sign(self):
        selected=[{"gene_id":"g","rank":1,"discovery_sign":-1}]
        details,ans=M.evaluate_set(selected,{"g":True},{"g":row("g",-1)},1)
        self.assertEqual(ans["agree"],1)
        self.assertEqual(details[0]["signed_external_t"],2)
    def test_dossier_fixed_pool_and_no_forced_positive(self):
        rows=[
          {"gene_id":"A","B0_rank":2,"evaluable":True,"agree":True,"external_q":.03,
           "external_CI_low":.1,"external_CI_high":2,"external_sign":1,"discovery_sign":1},
          {"gene_id":"B","B0_rank":1,"evaluable":True,"agree":True,"external_q":.04,
           "external_CI_low":.1,"external_CI_high":2,"external_sign":1,"discovery_sign":1},
          {"gene_id":"C","B0_rank":3,"evaluable":True,"agree":False,"external_q":.2,
           "external_CI_low":-2,"external_CI_high":-.1,"external_sign":-1,"discovery_sign":1},
          {"gene_id":"D","B0_rank":4,"evaluable":False,"agree":None,"external_q":None,
           "external_CI_low":None,"external_CI_high":None,"external_sign":None,"discovery_sign":1}]
        chosen=M.select_dossiers(rows)
        self.assertEqual([(r["gene_id"],r["dossier_category"]) for r in chosen],
                         [("B","supported_same_direction"),("A","supported_same_direction"),
                          ("C","strongly_conflicting"),("D","inconclusive")])
        self.assertEqual([r["dossier_category"] for r in M.select_dossiers(rows[2:])],
                         ["strongly_conflicting","inconclusive"])
    def test_opposite_CI_touching_zero_is_inconclusive(self):
        r={"gene_id":"Z","B0_rank":1,"evaluable":True,"agree":False,"external_q":.01,
           "external_CI_low":-1,"external_CI_high":0,"external_sign":-1,"discovery_sign":1}
        self.assertEqual(M.select_dossiers([r])[0]["dossier_category"],"inconclusive")
    def test_family_filter_contract_rejects_absent_passed_gene(self):
        with self.assertRaisesRegex(ValueError,"expression-pass family"):
            M.validate_family({"a":True,"b":True},[row("a")],[row("a")])
    def test_fitted_subset_contract(self):
        family=[row("a",q=1.),row("b",p="NA",q="NA",valid=False)]
        M.validate_family({"a":True,"b":True,"c":False},family,[family[0]])
        with self.assertRaisesRegex(ValueError,"de_fitted"):
            M.validate_family({"a":True,"b":True},family,family)



class SyntheticFileContract(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory(dir=M.ROOT / "work/synthetic", prefix="evaluation-contract-")
        self.base = pathlib.Path(self.tmp.name)
        self.discovery = self.base / "discovery"
        self.external = self.base / "external"
        self.discovery.mkdir(); self.external.mkdir()
        genes=[f"g{i:03}" for i in range(90)]
        directions=[{"gene_id":g,"logFC":"1","discovery_sign":"1"} for g in genes]
        products={}
        def write(path,rows,fields):
            M.write_tsv(path,rows,fields)
            return str(path.relative_to(M.ROOT))
        products["all_U_directions"]=write(self.discovery/"directions.tsv",directions,["gene_id","logFC","discovery_sign"])
        for method,start in (("B0",0),("B1",20),("G1",40)):
            ordered=genes[start:]+genes[:start]
            full=[{"rank":str(i+1),"gene_id":g,"method":method,"logFC":"1","discovery_sign":"1"} for i,g in enumerate(ordered)]
            for suffix,k in (("full",90),("top20",20),("top50",50)):
                products[f"{method}_{suffix}"]=write(self.discovery/f"{method}_{suffix}.tsv",full[:k],
                    ["rank","gene_id","method","logFC","discovery_sign"])
        members=[]
        for method in M.METHODS:
            own={r["gene_id"] for r in M.read_tsv(M.ROOT/products[f"{method}_top20"])}
            candidates=[g for g in genes if g not in own][:20]
            for draw in range(1,201):
                for pos,g in enumerate(candidates,1):
                    members.append({"method":method,"draw":draw,"position":pos,"gene_id":g,
                                    "logFC":1,"discovery_sign":1})
        products["null_memberships"]=write(self.discovery/"nulls.tsv",members,
            ["method","draw","position","gene_id","logFC","discovery_sign"])
        null_state=self.discovery/"null_state.json"
        M.write_json(null_state,{"methods":{m:{"status":"ESTIMABLE","sampled_sets":200} for m in M.METHODS}})
        products["null_status"]=str(null_state.relative_to(M.ROOT))
        config=self.discovery/"config.json"
        M.write_json(config,{"synthetic_only":True})
        products["statistical_config"]=str(config.relative_to(M.ROOT))
        counts=self.discovery/"counts_signature.txt"
        counts.write_text("invented fixture, no real counts\n")
        artifacts={str(p.relative_to(M.ROOT)):M.sha(p) for p in self.discovery.iterdir()}
        self.freeze=self.base/"freeze.json"
        M.write_json(self.freeze,{"phase":"SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE",
            "external_disease_effects_computed":False,"artifacts":artifacts,"products":products})
        # g000 absent, g001 filtered, g002 invalid, g003 zero -> 16/20 B0 evaluable.
        filters=[{"gene_id":g,"expression_pass":"FALSE" if g=="g001" else "TRUE"} for g in genes[1:]]
        self.family=[row(g,0 if g=="g003" else 1,p="NA" if g=="g002" else .5,
            q="NA" if g=="g002" else 44/87,valid=g!="g002") for g in genes[2:]]
        write(self.external/"expression_filter.tsv",filters,["gene_id","expression_pass"])
        write(self.external/"de_all.tsv",self.family,list(M.DE_COLUMNS))
        write(self.external/"de_fitted.tsv",self.family[1:],list(M.DE_COLUMNS))
        self.prov={"cohort":"heca","data_kind":"synthetic_test_only","scientific_effect_run":False,
             "tested_family_size":88,"estimable_genes":87,
             "input_sha256":{str(counts.relative_to(M.ROOT)):M.sha(counts),
                             str(self.freeze.relative_to(M.ROOT)):M.sha(self.freeze)},
             "output_sha256":{str(p):M.sha(p) for p in self.external.iterdir()}}
        M.write_json(self.external/"provenance.json",self.prov)
    def tearDown(self):
        self.tmp.cleanup()
    def test_full_synthetic_cli_600_frozen_sets_210_list_rows(self):
        import subprocess,sys
        output=self.base/"evaluation"
        run=subprocess.run([sys.executable,str(M.ROOT/"scripts/evaluate_replication.py"),
            "--discovery-freeze",str(self.freeze),"--external-dir",str(self.external),
            "--output-dir",str(output),"--synthetic-test"],capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stderr)
        result=M.read_json(output/"summary.json")
        self.assertEqual(result["methods"]["B0"]["20"]["evaluable"],16)
        self.assertEqual(result["methods"]["B0"]["20"]["missing_bounds"],[.8,1.])
        self.assertEqual(result["primary"]["missing_bounds"][0],0.)
        self.assertAlmostEqual(result["primary"]["missing_bounds"][1],.2)
        self.assertFalse(result["primary"]["superiority_claim"])
        self.assertEqual(len(M.read_tsv(output/"bootstrap_interface.tsv")),210)
        self.assertEqual(len(M.read_tsv(output/"null_sets.tsv")),600)
        self.assertEqual(result["dossier_pool_size"],40)
        selected=M.read_tsv(output/"dossier_selection.tsv")
        self.assertEqual([r["gene_id"] for r in selected],["g000"])
        self.assertNotIn("supported_same_direction",[r["dossier_category"] for r in selected])
        self.assertNotIn("strongly_conflicting",[r["dossier_category"] for r in selected])
    def test_output_tamper_rejected(self):
        with (self.external/"de_all.tsv").open("a") as f:f.write("tampered\n")
        manifest=M.read_json(self.freeze)
        with self.assertRaisesRegex(ValueError,"output hash"):
            M.validate_external(self.external,self.freeze,manifest,True)
    def test_original_input_signature_mismatch_rejected(self):
        self.prov["input_sha256"][next(iter(self.prov["input_sha256"]))]="0"*64
        (self.external/"provenance.json").write_text(__import__("json").dumps(self.prov))
        with self.assertRaisesRegex(ValueError,"original frozen input"):
            M.validate_external(self.external,self.freeze,M.read_json(self.freeze),True)
    def test_synthetic_cannot_enter_production(self):
        with self.assertRaisesRegex(ValueError,"pre-external discovery freeze"):
            M.verify_discovery(self.freeze,False)
    def test_null_membership_uses_own_sign_rejected(self):
        # Test the predeclared validation itself on artificial rows, after fixture hash gate.
        members=M.read_tsv(M.ROOT/M.read_json(self.freeze)["products"]["null_memberships"])
        members[0]["discovery_sign"]="-1"
        path=self.discovery/"nulls.tsv";path.unlink()
        M.write_tsv(path,members,["method","draw","position","gene_id","logFC","discovery_sign"])
        manifest=M.read_json(self.freeze)
        manifest["artifacts"][str(path.relative_to(M.ROOT))]=M.sha(path)
        self.freeze.write_text(__import__("json").dumps(manifest))
        freeze,lists,ranks,directions=M.verify_discovery(self.freeze,True)
        external={r["gene_id"]:r for r in self.family}
        with self.assertRaisesRegex(ValueError,"own frozen discovery prediction"):
            M.evaluate_nulls(freeze,directions,lists,{r["gene_id"]:True for r in self.family},external,{m:1. for m in M.METHODS})


if __name__=="__main__": unittest.main()
