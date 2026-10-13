"""Literal hand-calculated graph examples; no biological validation claims."""
import importlib.util,pathlib,unittest
from decimal import Decimal
ROOT=pathlib.Path(__file__).resolve().parents[1]
def implementation(test):
    p=ROOT/'scripts/map_string_genes.py'
    if not p.exists():test.fail('STRING gene mapper not implemented')
    spec=importlib.util.spec_from_file_location('map_string_genes',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
ANNOTATION=[
 {'hgnc_id':'HGNC:1','symbol':'GENEA','status':'Approved','ensembl_gene_id':'ENSG00000000001','prev_symbol':'OLDA'},
 {'hgnc_id':'HGNC:2','symbol':'GENEB','status':'Approved','ensembl_gene_id':'ENSG00000000002','prev_symbol':'OLD'},
 {'hgnc_id':'HGNC:3','symbol':'GENEC','status':'Approved','ensembl_gene_id':'ENSG00000000003','prev_symbol':'OLD'},
]
class StringGeneTests(unittest.TestCase):
    def test_many_proteins_same_gene_allowed_approved_and_unique_previous(self):
        m=implementation(self)
        info=[{'#string_protein_id':'9606.a','preferred_name':'GENEA'},{'#string_protein_id':'9606.b','preferred_name':'OLDA'}]
        x=m.map_proteins(info,ANNOTATION)
        self.assertEqual([r['decision'] for r in x],['include','include'])
        self.assertEqual([r['stable_gene_id'] for r in x],['ENSG00000000001','ENSG00000000001'])
        self.assertEqual(x[0]['source_preferred_name'],'GENEA')
        self.assertEqual(x[1]['mapping_method'],'exact_unique_HGNC_previous_symbol')
    def test_ambiguous_unknown_and_case_guess_are_excluded(self):
        m=implementation(self)
        names=['OLD','MISSING','genea']
        x=m.map_proteins([{'#string_protein_id':f'9606.p{j}','preferred_name':n} for j,n in enumerate(names)],ANNOTATION)
        self.assertEqual([r['decision'] for r in x],['exclude_ambiguous','exclude_unmapped','exclude_unmapped'])
    def test_current_approved_symbol_precedes_other_previous_symbol(self):
        m=implementation(self)
        annotation=ANNOTATION+[{'hgnc_id':'HGNC:4','symbol':'GENED','status':'Approved','ensembl_gene_id':'ENSG00000000004','prev_symbol':'GENEA'}]
        x=m.map_proteins([{'#string_protein_id':'9606.a','preferred_name':'GENEA'}],annotation)
        self.assertEqual(x[0]['stable_gene_id'],'ENSG00000000001')
    def test_gene_edge_max_collapse_loops_endpoints_and_mass(self):
        m=implementation(self)
        names={'a':'GENEA','b':'OLDA','c':'GENEB','d':'MISSING','e':'OLD'}
        x=m.map_proteins([{'#string_protein_id':'9606.'+k,'preferred_name':v} for k,v in names.items()],ANNOTATION)
        edges=[{'protein1':'9606.a','protein2':'9606.c','weight':'0.700'}, {'protein1':'9606.c','protein2':'9606.b','weight':'0.900'}, {'protein1':'9606.a','protein2':'9606.b','weight':'0.800'}, {'protein1':'9606.a','protein2':'9606.d','weight':'0.750'}, {'protein1':'9606.a','protein2':'9606.e','weight':'0.850'}]
        graph,stats=m.collapse_gene_edges(edges,x)
        self.assertEqual(graph,{('ENSG00000000001','ENSG00000000002'):Decimal('0.900')})
        self.assertEqual(stats['source_protein_edges'],5)
        self.assertEqual(stats['endpoint_excluded_edges'],2)
        self.assertEqual(stats['gene_self_loops_excluded'],1)
        self.assertEqual(stats['duplicate_gene_edges_collapsed'],1)
        self.assertEqual(stats['source_weight_sum'],Decimal('4.000'))
        self.assertEqual(stats['endpoint_excluded_weight_sum'],Decimal('1.600'))
        self.assertEqual(stats['gene_self_loop_weight_sum'],Decimal('0.800'))
        self.assertEqual(stats['duplicate_collapse_weight_loss'],Decimal('0.700'))
        self.assertEqual(stats['output_gene_weight_sum'],Decimal('0.900'))
    def test_unknown_graph_endpoint_and_invalid_weight_fail(self):
        m=implementation(self);x=m.map_proteins([{'#string_protein_id':'9606.a','preferred_name':'GENEA'}],ANNOTATION)
        with self.assertRaises(ValueError):m.collapse_gene_edges([{'protein1':'9606.a','protein2':'9606.unknown','weight':'0.800'}],x)
        for w in ['nan','-0.1','1.001','0.699']:
            with self.subTest(weight=w),self.assertRaises(ValueError):m.collapse_gene_edges([{'protein1':'9606.a','protein2':'9606.a','weight':w}],x)
    def test_duplicate_protein_ID_never_overwrites_mapping(self):
        m=implementation(self)
        with self.assertRaises(ValueError):m.map_proteins([{'#string_protein_id':'9606.a','preferred_name':'GENEA'},{'#string_protein_id':'9606.a','preferred_name':'GENEB'}],ANNOTATION)
if __name__=='__main__':unittest.main()
