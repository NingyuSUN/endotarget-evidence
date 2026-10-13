"""Ingestion tests use hand-calculated synthetic inputs, never scientific evidence."""
import importlib.util, pathlib, unittest, hashlib, base64
ROOT=pathlib.Path(__file__).resolve().parents[1]
def load_script(name,test):
    p=ROOT/'scripts'/f'{name}.py'
    if not p.exists():test.fail(f'{name} not implemented')
    spec=importlib.util.spec_from_file_location(name,p);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
class CountsTests(unittest.TestCase):
    def test_restore_drops_entire_many_source_group_and_preserves_fraction(self):
        m=load_script('prepare_discovery_counts',self)
        mapping=[{'source_row':'0','source_gene':'ENSG00000000001.4','canonical_gene':'ONE'}, {'source_row':'1','source_gene':'GeneB','canonical_gene':'TWO'}, {'source_row':'2','source_gene':'ENSG00000000003.1','canonical_gene':'MERGED'}, {'source_row':'3','source_gene':'ENSG00000000004.2','canonical_gene':'MERGED'}]
        rows=[{'gene':'MERGED','sample1':'9.875'}, {'gene':'TWO','sample1':'2.0625'}, {'gene':'ONE','sample1':'0.125'}]
        restored,dropped=m.restore_count_rows(mapping,rows,['sample1'])
        self.assertEqual(restored,[{'source_feature_id':'ENSG00000000001.4','sample1':'0.125'},{'source_feature_id':'GeneB','sample1':'2.0625'}])
        self.assertEqual([r['gene'] for r in dropped],['MERGED'])
    def test_duplicate_stable_target_discards_both_rows(self):
        m=load_script('prepare_discovery_counts',self)
        annotation=[{'hgnc_id':'HGNC:1','symbol':'GENEA','status':'Approved','ensembl_gene_id':'ENSG00000000001','prev_symbol':'OLDA'}, {'hgnc_id':'HGNC:2','symbol':'GENEB','status':'Approved','ensembl_gene_id':'ENSG00000000002','prev_symbol':''}]
        mapped=m.map_features(['ENSG00000000001.4','GENEA','GENEB'],annotation)
        self.assertEqual([r['decision'] for r in mapped],['exclude_duplicate_target','exclude_duplicate_target','include'])
        self.assertEqual(mapped[2]['stable_gene_id'],'ENSG00000000002')
        self.assertEqual(mapped[0]['source_version'],'4')
    def test_ambiguous_previous_symbol_and_unmapped_never_guess(self):
        m=load_script('prepare_discovery_counts',self)
        annotation=[{'hgnc_id':'HGNC:1','symbol':'GENEA','status':'Approved','ensembl_gene_id':'ENSG00000000001','prev_symbol':'OLD'}, {'hgnc_id':'HGNC:2','symbol':'GENEB','status':'Approved','ensembl_gene_id':'ENSG00000000002','prev_symbol':'OLD'}]
        mapped=m.map_features(['OLD','MISSING','GENEB'],annotation)
        self.assertEqual([r['decision'] for r in mapped],['exclude_ambiguous_symbol','exclude_unmapped','include'])
    def test_exact_approved_symbol_precedes_other_records_previous_symbol(self):
        m=load_script('prepare_discovery_counts',self)
        annotation=[{'hgnc_id':'HGNC:1','symbol':'GENEA','status':'Approved','ensembl_gene_id':'ENSG00000000001','prev_symbol':''}, {'hgnc_id':'HGNC:2','symbol':'GENEB','status':'Approved','ensembl_gene_id':'ENSG00000000002','prev_symbol':'GENEA'}]
        mapped=m.map_features(['GENEA'],annotation)
        self.assertEqual(mapped[0]['decision'],'include');self.assertEqual(mapped[0]['stable_gene_id'],'ENSG00000000001')
    def test_count_validation_rejects_negative_and_nonfinite(self):
        m=load_script('prepare_discovery_counts',self)
        for value in ['-0.125','nan','inf']:
            with self.assertRaises(ValueError):m.count_value(value)
    def test_actual_pinned_C18orf21_is_ambiguous_and_never_kept(self):
        m=load_script('prepare_discovery_counts',self)
        annotation=ROOT/'data/reference/gene_annotation/hgnc_complete_set_2026-10-06.tsv'
        if not annotation.exists():self.skipTest('Fetch pinned HGNC annotation before real-data ambiguity audit')
        rows=m.table(annotation)
        result=m.map_features(['C18orf21'],rows)[0]
        self.assertEqual(result['decision'],'exclude_ambiguous_symbol')
        self.assertEqual(set(result['hgnc_id'].split('|')),{'HGNC:28802','HGNC:57042'})
class AnnotationTests(unittest.TestCase):
    def test_corrupt_reference_rejected_even_if_schema_valid(self):
        m=load_script('fetch_gene_annotation',self)
        blob=b'hgnc_id\tsymbol\tstatus\tensembl_gene_id\tprev_symbol\nHGNC:1\tGENEA\tApproved\tENSG00000000001\t\n'
        pin={'bytes':len(blob),'sha256':hashlib.sha256(blob).hexdigest(),'official_GCS_md5':base64.b64encode(hashlib.md5(blob).digest()).decode()}
        self.assertEqual(m.validate_annotation(blob,pin),1)
        with self.assertRaises(ValueError):m.validate_annotation(blob.replace(b'GENEA',b'GENEX'),pin)
if __name__=='__main__':unittest.main()
