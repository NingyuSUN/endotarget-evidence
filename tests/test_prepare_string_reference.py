"""Hand-calculated STRING channel boundaries; synthetic edges are not biology."""
import importlib.util,json,pathlib,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('prepare_string_reference',ROOT/'scripts/prepare_string_reference.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)
HEADER='protein1 protein2 experiments database experiments_transferred database_transferred textmining combined_score'
class StringFilterTests(unittest.TestCase):
    def setUp(self):
        self.cfg=json.loads((ROOT/'config/network_reference.json').read_text())
        self.proteins={'9606.a','9606.b','9606.c','9606.d'}
        self.assertTrue(hasattr(MODULE,'filter_edges'),'pure streaming edge filter not implemented')
    def run_filter(self,rows,cfg=None):
        return MODULE.filter_edges(iter([HEADER]+rows),self.proteins,cfg or self.cfg)
    def test_combined_transferred_and_text_only_cannot_rescue_low_direct(self):
        edges,stats=self.run_filter(['9606.a 9606.b 699 0 1000 1000 1000 1000','9606.a 9606.c 0 0 1000 1000 0 1000','9606.a 9606.d 0 0 0 0 1000 1000'])
        self.assertEqual(edges,{})
        self.assertEqual(stats['below_direct_score_threshold'],3)
    def test_700_exact_boundary_uses_either_direct_channel(self):
        edges,stats=self.run_filter(['9606.a 9606.b 699 0 0 0 0 1000','9606.a 9606.c 700 0 0 0 0 700','9606.a 9606.d 0 700 0 0 0 700'])
        self.assertEqual(edges,{('9606.a','9606.c'):700,('9606.a','9606.d'):700})
        self.assertEqual(stats['qualifying_directed_rows'],2)
    def test_reverse_and_same_order_duplicate_edges_take_max_permitted_score(self):
        edges,stats=self.run_filter(['9606.a 9606.b 700 0 1000 0 0 1000','9606.b 9606.a 800 900 0 0 0 900','9606.a 9606.b 750 800 0 0 0 850'])
        self.assertEqual(edges,{('9606.a','9606.b'):900})
        self.assertEqual(stats['qualifying_directed_rows'],3)
    def test_qualifying_self_loop_is_excluded(self):
        edges,stats=self.run_filter(['9606.a 9606.a 1000 0 0 0 0 1000'])
        self.assertEqual(edges,{})
        self.assertEqual(stats['self_loops_excluded'],1)
    def test_unknown_or_nonhuman_qualifying_protein_fails(self):
        for protein in ['9606.unknown','10090.a']:
            with self.subTest(protein=protein),self.assertRaises(ValueError):
                self.run_filter([f'9606.a {protein} 700 0 0 0 0 700'])
        edges,_=self.run_filter(['9606.a 9606.unknown 0 0 1000 1000 1000 1000'])
        self.assertEqual(edges,{})
    def test_unsupported_configuration_semantics_fail_instead_of_silent_ignoring(self):
        changes=[('permitted_direct_channels',['experiments','database','textmining']),('source_score_scale',100),('minimum_max_direct_score',699),('edge_weight','combined_score/1000'),('undirected_duplicate_policy','first'),('self_loop_policy','keep'),('species_taxid',10090)]
        for key,value in changes:
            with self.subTest(key=key),self.assertRaises(ValueError):
                self.run_filter(['9606.a 9606.b 700 0 0 0 0 700'],dict(self.cfg,**{key:value}))
    def test_invalid_direct_score_and_malformed_row_fail(self):
        for line in ['9606.a 9606.b 1001 0 0 0 0 1000','9606.a 9606.b -1 800 0 0 0 800','9606.a 9606.b 700']:
            with self.subTest(line=line),self.assertRaises(ValueError):self.run_filter([line])
if __name__=='__main__':unittest.main()
