import pathlib, sys, unittest, tempfile, types, csv, shutil, contextlib, io
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"))
class LodoTests(unittest.TestCase):
    def test_new_implementation_exists(self):
        self.assertTrue(pathlib.Path("scripts/run_discovery_lodo.py").exists())
    def test_missing_genes_count_as_losses(self):
        from run_discovery_lodo import top20_metrics
        ref=[str(i) for i in range(20)]
        got=[str(i) for i in range(10)]+[str(i) for i in range(20,30)]
        result=top20_metrics(ref,got)
        self.assertEqual(result["retention"],.5)
        self.assertAlmostEqual(result["jaccard"],10/30)
    def test_reject_incomplete_reference(self):
        from run_discovery_lodo import top20_metrics
        with self.assertRaises(ValueError): top20_metrics(["a"],["a"])
    def test_failed_folds_are_not_replaced_and_changed_inputs_fail_resume(self):
        from run_discovery_lodo import run
        root=pathlib.Path("work/synthetic");root.mkdir(parents=True,exist_ok=True)
        temp=pathlib.Path(tempfile.mkdtemp(prefix="lodo_minima_",dir=root))
        try:
            baseline=temp/"baseline";baseline.mkdir()
            for method in ("B0","B1","G1"):
                (baseline/f"{method}_top20.tsv").write_text("gene_id\n"+"\n".join(f"ENSG{i:011d}" for i in range(20))+"\n")
            clinical=temp/"clinical.tsv"
            with clinical.open("w",newline="") as f:
                w=csv.writer(f,delimiter="\t");w.writerow(["donor_id","disease","primary_lineage_cells"])
                w.writerows((f"d{i}","ENDO" if i<8 else "non-ENDO",100) for i in range(16))
            counts=temp/"counts.tsv";counts.write_text("gene_id\n")
            edges=temp/"edges.tsv";edges.write_text("gene1\tgene2\tweight\n")
            args=types.SimpleNamespace(synthetic_test=True,counts=counts,clinical=clinical,baseline_dir=baseline,edges=edges,output_dir=temp/"out",config=pathlib.Path("config/statistical_analysis.json"),ranking_config=pathlib.Path("config/ranking_algorithm.json"),discovery_freeze=pathlib.Path("unused"),pre_effect_freeze=pathlib.Path("unused"))
            with contextlib.redirect_stdout(io.StringIO()):
                first=run(args);second=run(args)
            self.assertEqual(first["failed_folds"],16)
            self.assertEqual(first,second)
            self.assertFalse(first["failed_folds_replaced"])
            counts.write_text("gene_id\nENSG00000000000\n")
            with self.assertRaisesRegex(ValueError,"Unsafe LODO resume"):run(args)
        finally:shutil.rmtree(temp)
if __name__=="__main__": unittest.main()

