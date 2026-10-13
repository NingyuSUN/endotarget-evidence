import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"))
class SensitivityReportingTests(unittest.TestCase):
    def test_new_module_exists(self):self.assertTrue(pathlib.Path("scripts/summarize_clinical_sensitivities.py").exists())
    def test_unknown_day_stays_missing_and_day28_is_retained(self):
        from summarize_clinical_sensitivities import ascertainment_row
        base={"donor_id":"d1","sample_id":"s1","disease":"non-ENDO","age":"30","phase_original":"Proliferative","histology_phase_original":"Secretory","expression_phase_original":"","day_original":"unknown","histology_conflict":"True","adenomyosis":"0","hormone_h5_original":"nan","hormone_clinical_original":"NA","hormone_eligibility_basis":"study_defined_natural_cycle_group"}
        row=ascertainment_row("heca",base,None);self.assertIsNone(row["cycle_day_numeric"]);self.assertEqual(row["cycle_day_original"],"unknown");self.assertFalse(row["strict_histology_membership"])
        base["day_original"]="28";row=ascertainment_row("heca",base,None);self.assertEqual(row["cycle_day_numeric"],28);self.assertEqual(row["author_phase_original"],"Proliferative")
    def test_discovery_expression_phase_is_not_invented(self):
        from summarize_clinical_sensitivities import ascertainment_row
        cl={"donor_id":"d","sample_id":"s","disease":"ENDO","age":"31","strict_histology_sensitivity":"True","adenomyosis_free_sensitivity":"True","author_phase":"proliferative","histology_phase_original":"NA","hormone_eligibility_basis":"explicit_FALSE"}
        candidate={"Patient.age":"31","Menstrual cycle day (Base 1)":"9","Progesterone in serum (nmol/L)":".3","pathologist HE evaluation":"NA","MenstrualCyclePhase":"proliferative","Hormonal treatment/IUD":"FALSE"}
        row=ascertainment_row("discovery",cl,candidate);self.assertEqual(row["expression_phase_original"],"");self.assertEqual(row["cycle_day_numeric"],9);self.assertEqual(row["histology_original"],"NA")
    def test_missing_gene_is_retained_not_an_opposite_result(self):
        from summarize_clinical_sensitivities import paired_gene_details,summarize_pairs
        p={"g1":{"gene_id":"g1","logFC":"1","t":"2","P.Value":".1","fit_estimable":"TRUE"},"g2":{"gene_id":"g2","logFC":"-1","t":"-2","P.Value":".1","fit_estimable":"TRUE"}}
        e={"g1":{"gene_id":"g1","logFC":"-2","t":"-3","P.Value":".1","fit_estimable":"TRUE"}}
        rows=paired_gene_details(p,{"g1":True,"g2":False},e);summary=summarize_pairs(rows,"all")
        self.assertEqual(len(rows),2);self.assertEqual(summary["direction_evaluable"],1);self.assertEqual(summary["opposite_direction"],1);self.assertEqual(summary["unevaluable"],1)
if __name__=="__main__":unittest.main()
