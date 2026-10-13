import io,sys,pathlib,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"))
import prepare_heca_counts as p
class BoundaryTests(unittest.TestCase):
 def test_genuine_selected_sum(self):
  src=io.BytesIO(b"%%MatrixMarket matrix coordinate integer general\n3 2 4\n1 1 2\n3 1 5\n1 2 7\n2 2 11\n")
  got=p.aggregate_mtx(src,3,[True,False])
  self.assertEqual(got["gene_counts"].tolist(),[2,0,5])
  self.assertEqual(got["cell_counts"].tolist(),[7,18])
  self.assertEqual(got["selected_total"],7)
 def test_missing_selected_barcode_rejected(self):
  with self.assertRaises(ValueError):p.join_barcodes(["AAAAAAAAAAAAAAAA-1"],["LIB_CCCCCCCCCCCCCCCC"],"LIB")
 def test_collapsed_gem_suffix_ambiguity_rejected(self):
  with self.assertRaises(ValueError):p.join_barcodes(["AAAAAAAAAAAAAAAA-1","AAAAAAAAAAAAAAAA-2"],["LIB_AAAAAAAAAAAAAAAA"],"LIB")
 def test_duplicate_selected_cell_rejected(self):
  with self.assertRaises(ValueError):p.join_barcodes(["AAAAAAAAAAAAAAAA-1"],["LIB_AAAAAAAAAAAAAAAA"]*2,"LIB")
 def test_log_or_real_layer_rejected(self):
  with self.assertRaises(ValueError):p.aggregate_mtx(io.BytesIO(b"%%MatrixMarket matrix coordinate real general\n1 1 1\n1 1 0.5\n"),1,[True])
 def test_negative_raw_count_rejected(self):
  with self.assertRaises(ValueError):p.aggregate_mtx(io.BytesIO(b"%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 -1\n"),1,[True])
 def test_duplicate_coordinate_rejected(self):
  with self.assertRaises(ValueError):p.aggregate_mtx(io.BytesIO(b"%%MatrixMarket matrix coordinate integer general\n1 1 2\n1 1 2\n1 1 3\n"),1,[True])
 def test_wrong_feature_identity_rejected(self):
  with self.assertRaises(ValueError):p.check_features([["ENSG1","A","Gene Expression"]],[["ENSG2","A","Gene Expression"]])

class HardenedCheckpointTests(unittest.TestCase):
 def setUp(self):
  import tempfile,numpy as np
  base=p.ROOT/"data/processed/heca/test_scratch";base.mkdir(parents=True,exist_ok=True)
  self.tmp=tempfile.TemporaryDirectory(dir=base);self.path=pathlib.Path(self.tmp.name)/"counts.npz"
  np.savez(self.path,gene_counts=np.array([2.,0.,5.]),cell_counts=np.array([7.,18.]))
  self.pins={"archive_sha256":"archive","selected_barcodes_sha256":"barcodes","source_features_sha256":"features","parser_sha256":"parser"}
  self.meta={**self.pins,"npz_sha256":p.sha(self.path),"observed_nnz":4,"declared_nnz":4,"source_shape":[3,2],"selected_total":7,"all_total":25}
 def tearDown(self):self.tmp.cleanup()
 def test_correct_checkpoint_accepted(self):
  import numpy as np
  with np.load(self.path) as a:self.assertTrue(p.validate_checkpoint(self.meta,a,self.pins,self.path,3,np.array([True,False])))
 def test_mutated_npz_rejected(self):
  import numpy as np
  np.savez(self.path,gene_counts=np.array([3.,0.,4.]),cell_counts=np.array([7.,18.]))
  with np.load(self.path) as a:
   with self.assertRaises(ValueError):p.validate_checkpoint(self.meta,a,self.pins,self.path,3,np.array([True,False]))
 def test_changed_selected_barcode_pin_rejected(self):
  import numpy as np
  changed={**self.pins,"selected_barcodes_sha256":"different"}
  with np.load(self.path) as a:
   with self.assertRaises(ValueError):p.validate_checkpoint(self.meta,a,changed,self.path,3,np.array([True,False]))
 def test_fractional_checkpoint_rejected_even_with_current_hash(self):
  import numpy as np
  np.savez(self.path,gene_counts=np.array([2.5,0.,4.5]),cell_counts=np.array([7.,18.]))
  self.meta["npz_sha256"]=p.sha(self.path)
  with np.load(self.path) as a:
   with self.assertRaises(ValueError):p.validate_checkpoint(self.meta,a,self.pins,self.path,3,np.array([True,False]))
 def test_deterministic_gzip(self):
  one=pathlib.Path(self.tmp.name)/"one.tsv.gz";two=pathlib.Path(self.tmp.name)/"two.tsv.gz"
  p.write_tsv(one,[{"gene_id":"ENSG1","donor":"7"}],["gene_id","donor"]);p.write_tsv(two,[{"gene_id":"ENSG1","donor":"7"}],["gene_id","donor"])
  self.assertEqual(p.sha(one),p.sha(two))

if __name__=="__main__":unittest.main()

