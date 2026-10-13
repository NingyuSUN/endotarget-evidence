"""Boundary checks for reporting, using invented counts only; never fit a model."""
import importlib.util, json, tempfile, unittest
from unittest.mock import patch
import contextlib, io
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if ROOT.name=="work": ROOT=ROOT.parent
source=ROOT/"scripts/build_replication_report.py"
if not source.exists(): source=ROOT/"work/report_builder_draft/build_replication_report.py"
spec=importlib.util.spec_from_file_location("report_builder",source)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
class ReportingBoundaryTests(unittest.TestCase):
    def setUp(self):
        base=ROOT/"work/synthetic/report_unit"; base.mkdir(parents=True,exist_ok=True)
        self.tmp=tempfile.TemporaryDirectory(dir=base); self.p=Path(self.tmp.name)
    def tearDown(self): self.tmp.cleanup()
    def bootstrap(self,failed):
        rs=[dict(replicate=i,valid=i>failed,G1_minus_B0_top20="NA" if i<=failed else 0,reason="singular" if i<=failed else "") for i in range(1,501)]
        mod.write_tsv(self.p/"paired_differences.tsv",rs,rs[0].keys())
        s=dict(total_required=500,completed=500,valid=500-failed,invalid=failed,interval_suppressed=failed>50,interval=None if failed>50 else [0,0],synthetic=True)
        (self.p/"summary.json").write_text(json.dumps(s)); return s
    def test_50_failures_permits_conditional_interval(self):
        self.bootstrap(50); s,rs=mod.bootstrap_summary(self.p,True)
        self.assertFalse(s["interval_suppressed"]); self.assertEqual(len(rs),500)
    def test_51_failures_suppresses_interval(self):
        self.bootstrap(51); s,_=mod.bootstrap_summary(self.p,True)
        self.assertTrue(s["interval_suppressed"]); self.assertIsNone(s["interval"])
    def test_suppressed_interval_cannot_be_reported(self):
        s=self.bootstrap(51); s["interval"]=[-.1,.1]; (self.p/"summary.json").write_text(json.dumps(s))
        with self.assertRaisesRegex(ValueError,"Suppressed"): mod.bootstrap_summary(self.p,True)
    def test_false_summary_counts_rejected(self):
        s=self.bootstrap(1); s["valid"]=500; (self.p/"summary.json").write_text(json.dumps(s))
        with self.assertRaisesRegex(ValueError,"counts mismatch"): mod.bootstrap_summary(self.p,True)
    def lodo(self,missing=False):
        rs=[]
        for i in range(1,23):
            for m in mod.METHODS:
                rs.append(dict(fold=i,dropped_donor="fixture"+str(i),method=m,valid=i!=3,reason="design fails" if i==3 else "",intersection=15,retention=.75,jaccard=.6))
        if missing: rs.pop()
        mod.write_tsv(self.p/"fold_metrics.tsv",rs,rs[0].keys())
        (self.p/"summary.json").write_text(json.dumps(dict(total_folds=22,valid_folds=21,failed_folds=1,failed_folds_replaced=False,synthetic=True)))
    def test_one_failed_lodo_fold_preserved(self):
        self.lodo(); s,rs=mod.lodo_summary(self.p,True)
        self.assertEqual(s["failed_folds"],1); self.assertEqual(len(rs),66)
        self.assertEqual(sum(not mod.truth(r["valid"]) for r in rs),3)
    def test_missing_fold_cannot_disappear(self):
        self.lodo(True)
        with self.assertRaisesRegex(ValueError,"Incomplete"): mod.lodo_summary(self.p,True)
    def test_retention_has_fixed20_denominator(self):
        self.lodo(); rs=mod.rows(self.p/"fold_metrics.tsv"); rs[0]["retention"]=".8"
        mod.write_tsv(self.p/"fold_metrics.tsv",rs,rs[0].keys())
        with self.assertRaisesRegex(ValueError,"denominator"): mod.lodo_summary(self.p,True)
    def test_no_forced_supported_dossier(self):
        pool=[dict(gene_id="fixtureB",B0_rank=2,dossier_category="inconclusive"),dict(gene_id="fixtureA",B0_rank=1,dossier_category="inconclusive")]
        selected=mod.dossier_choices(pool)
        self.assertEqual([r["gene_id"] for r in selected],["fixtureA"])
    def test_default_main_final_receipt(self):
        with patch("sys.argv", ["build_replication_report.py"]), patch.object(mod,"build",return_value={}) as build, contextlib.redirect_stdout(io.StringIO()):
            mod.main()
        self.assertEqual(build.call_args.args[0].receipt,"receipts/replication_report.json")
    def test_default_main_partial_receipt(self):
        with patch("sys.argv", ["build_replication_report.py","--partial"]), patch.object(mod,"build",return_value={}) as build, contextlib.redirect_stdout(io.StringIO()):
            mod.main()
        self.assertEqual(build.call_args.args[0].receipt,"receipts/replication_draft_report.json")
    def test_explicit_receipt_override_preserved(self):
        with patch("sys.argv", ["build_replication_report.py","--partial","--receipt","work/synthetic/custom_receipt.json"]), patch.object(mod,"build",return_value={}) as build, contextlib.redirect_stdout(io.StringIO()):
            mod.main()
        self.assertEqual(build.call_args.args[0].receipt,"work/synthetic/custom_receipt.json")
    def test_nonfinite_is_unavailable(self):
        self.assertIsNone(mod.num("nan")); self.assertIsNone(mod.num("inf")); self.assertEqual(mod.num("0"),0)
if __name__=="__main__": unittest.main()

