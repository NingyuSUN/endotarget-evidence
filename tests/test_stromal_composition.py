"""Boundary tests for descriptive proportions; no biological validation."""
import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location("composition",Path(__file__).resolve().parents[1]/"scripts/summarize_stromal_composition.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class CompositionTests(unittest.TestCase):
    def selected(self):
        return [{"donor_id":"d1","sample_id":"lib1","disease":"ENDO","expected":"2"},
                {"donor_id":"d2","sample_id":"lib2","disease":"non-ENDO","expected":"1"}]
    def cells(self):
        return [{"sample":"lib1","donor":"d1","label":"A","id":"c1","QC":"no"},
                {"sample":"lib1","donor":"d1","label":"B","id":"c2","QC":"no"},
                {"sample":"lib2","donor":"d2","label":"A","id":"c3","QC":"no"},
                {"sample":"lib1_extra","donor":"d1","label":"A","id":"c4","QC":"no"},
                {"sample":"lib1","donor":"d1","label":"A","id":"c5","QC":"yes"},
                {"sample":"lib1","donor":"d1","label":"nan","id":"c6","QC":"no"}]
    def scan(self,cells=None):
        return m.scan_cells(self.cells() if cells is None else cells,self.selected(),["A","B"],"label","id",donor_field="donor",qc_field="QC")
    def test_unselected_library_and_lowQC_are_not_counted(self):
        x=self.scan();self.assertEqual(dict(x["counts"]["d1"]),{"A":1,"B":1});self.assertEqual(x["positions"],[0,1,2])
    def test_complete_grid_uses_primary_denominator_and_real_zero(self):
        x=self.scan();m.reconcile(x,self.selected(),"expected",3)
        rr=m.donor_grid("test",self.selected(),["A","B"],x,"label")
        self.assertEqual(len(rr),4);self.assertEqual(rr[0]["proportion_primary"],0.5)
        self.assertEqual(rr[-1]["cells"],0);self.assertEqual(rr[-1]["proportion_primary"],0)
    def test_missing_label_remains_visible_and_not_inferred(self):
        self.assertEqual(self.scan()["missing"]["d1"],1)
    def test_duplicate_selected_cell_is_rejected(self):
        cells=self.cells();cells.append(cells[0])
        with self.assertRaises(ValueError):self.scan(cells)
    def test_wrong_donor_crosswalk_is_rejected(self):
        cells=self.cells();cells[0]["donor"]="d2"
        with self.assertRaises(ValueError):self.scan(cells)
    def test_repeated_donor_in_selected_ledger_is_rejected(self):
        s=self.selected();s[1]["donor_id"]="d1"
        with self.assertRaises(ValueError):m.validate_selection(s)
    def test_exact_donor_denominator_drift_fails(self):
        x=self.scan();s=self.selected();s[0]["expected"]="3"
        with self.assertRaises(ValueError):m.reconcile(x,s,"expected",3)
    def test_equal_donor_weighting_differs_from_pooled_cells(self):
        rr=[{"cohort":"test","disease":"ENDO","source_label_field":"label","source_subtype_label":"A",
             "donor_id":"d1","cells":1,"primary_cells":2,"proportion_primary":0.5},
            {"cohort":"test","disease":"ENDO","source_label_field":"label","source_subtype_label":"A",
             "donor_id":"d2","cells":9,"primary_cells":9,"proportion_primary":1.0}]
        out=m.group_summary(rr)[0]
        self.assertEqual(out["mean_proportion_equal_donor"],0.75)
        self.assertAlmostEqual(out["pooled_cell_proportion"],10/11)

if __name__=="__main__":unittest.main()
