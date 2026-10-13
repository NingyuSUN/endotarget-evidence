#!/usr/bin/env python3
"""Add explicit source/tissue context to audited rows without rescanning QTL XML."""
import csv,datetime,hashlib,json,re
from pathlib import Path
from genetic_evidence_audit_helpers import digest,write_tsv

def main():
 audit=json.load(open('receipts/genetic_evidence_audit.json'))
 for p,h in audit['output_hashes'].items():assert digest(p)==h
 source=audit['source'];paper='10.1038/s41588-023-01323-z'
 items=[]
 for name,path,url,status,note in [
  ('2023 primary article','data/genetics/pmc.txt','https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/','accepted','Retrieval UTC from saved-file mtime; publication is author manuscript'),
  ('2023 supplementary tables',source['path'],source['url'],'accepted','Publisher XLSX; full table 1-35 object; only defined relevant tables parsed'),
  ('PMC supplement failed response','data/genetics/rahmioglu2023_supplement.xlsx','https://pmc.ncbi.nlm.nih.gov/articles/instance/10042257/bin/NIHMS1873427-supplement-Supplementary_Materials.xlsx','rejected_html','Preparing to download HTML, not accepted as XLSX'),
  ('Nature article failed response','data/genetics/article.txt','https://www.nature.com/articles/s41588-023-01323-z','rejected_html','Challenge HTML; not article evidence')]:
  p=Path(path);items.append(dict(name=name,publication=paper,url=url,path=path,bytes=p.stat().st_size,sha256=digest(p),retrieved_utc=source['retrieved_utc'] if path==source['path'] else datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat(),status=status,note=note))
 write_tsv('manifests/genetic_sources.tsv',items)
 hits=list(csv.DictReader(open('reports/genetic_evidence_rows.tsv'),delimiter='\t'));context=[]
 for h in hits:
  c=json.loads(h['source_cells_json']);v={re.sub(r'\d+$','',k):val for k,val in c.items()};n=int(re.match(r'Supp(\d+)',h['source_sheet'])[1])
  if n==11:
   extra=dict(locus='MLLT10/10p12.31',variant=v['C'],index_variant=v['A'],genome_coordinates=v['D'],genome_build='hg19 as published',ancestry='European GWAS fine-mapping (58,961 cases/700,352 controls)',tissue='not applicable: genomic annotation',mapping_method='author Nearest gene column; VEP consequence '+v['G'],source_analysis='99% credible set at single-trait endometriosis locus',effect='',p_value='',interpretation='positional support only; 99% refers to set coverage, not this variant posterior or gene causal probability')
  elif n==31:
   extra=dict(locus='MLLT10/10p12.31',variant=v['A'],index_variant='rs7098100 (ENDO MTAG lead in main Table 2; proxy relationship per Supp31 title)',genome_coordinates=v['B']+':'+v['C'],genome_build='hg19 as published',ancestry='UK/Australian endometrial QTL; source European GWAS analyses; QTL per-person ancestry not verified',tissue='bulk endometrium, Australian+UK meta-analysis',mapping_method='published eQTL lookup at MTAG shared endometriosis/pain loci',source_analysis='MTAG-shared-locus QTL lookup; not independent single-trait replication',effect='A1='+v['D']+';A2='+v['E']+';beta='+v['H']+';SE='+v['I']+';published orientation='+v['G'],p_value=v['J'],interpretation='limited regulatory support; nominal table p-value; no pool-gene SMR/HEIDI or formal colocalization result; treatment direction unknown')
  else:raise AssertionError('Unexpected hit needs explicit evidence context: '+h['source_sheet'])
  context.append({**h,**extra,'source_url':source['url'],'cohort_reuse':'same GWAS/QTL publication; no expression-cohort reuse documented; individual nonoverlap unverified','independent_publications':1,'causal_claim':'none'})
 write_tsv('reports/genetic_evidence_context.tsv',context)
 assert len(hits)==2 and {h['gene_id'] for h in hits}=={'ENSG00000180592'}
 pool=list(csv.DictReader(open('reports/genetic_evidence_frozen_pool.tsv'),delimiter='\t'));actual={r['gene_id'] for r in pool};fr=json.load(open('receipts/discovery_rank_freeze.json'));expected=set()
 for m in ['B0','G1']:
  p=fr['products'][m+'_top20'];assert digest(p)==fr['artifacts'][p];expected.update(r['gene_id'] for r in csv.DictReader(open(p),delimiter='\t'))
 assert actual==expected and len(actual)==40
 rejected=list(csv.DictReader(open('reports/genetic_evidence_rejected_prose_hits.tsv'),delimiter='\t'));assert {r['gene_id'] for r in rejected}=={'ENSG00000075624','ENSG00000146648'}
 summary=list(csv.DictReader(open('reports/genetic_evidence_summary.tsv'),delimiter='\t'));assert len(summary)==40 and {r['gene_id'] for r in summary}==actual
 q=list(csv.DictReader(open('reports/genetic_evidence_uk_qtl_context.tsv'),delimiter='\t'))
 results=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pool_membership_reconstructed_and_equal=True,all_40_rows_present=True,true_mapping_gene='SKIDA1',actual_source_cells=['Supp11.CrediblesetSNPannot!J76','Supp31.EndometriumeQTLs-shared!K11'],rejected_prose_mentions=['ACTB','EGFR'],n_pool_genes_in_UK_cis_QTL_context=sum(int(r['tested_variant_rows'])>0 for r in q),UK_QTL_context_variant_rows=sum(int(r['tested_variant_rows']) for r in q),source_and_output_hashes={p:digest(p) for p in ['manifests/genetic_sources.tsv','reports/genetic_evidence_context.tsv']},rankings_unchanged=True,scientific_causal_validation=False)
 Path('receipts/genetic_evidence_validation.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
if __name__=='__main__':main()
