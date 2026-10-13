#!/usr/bin/env python3
"""Audit the fixed 40-gene pool against exact published 2023 XLSX columns.
Run from the SRV01 project root; never reads external effects or writes rankings.
"""
import csv,datetime,hashlib,json,re,time
from pathlib import Path
from collections import defaultdict
from genetic_evidence_xlsx import Workbook

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def write_tsv(p,rows,fields=None):
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def vals(c):return {re.sub(r'\d+$','',k):v for k,v in c.items()}
start=time.monotonic();pool_path='reports/genetic_evidence_frozen_pool.tsv';pool=list(csv.DictReader(open(pool_path),delimiter='\t'))
freeze=json.load(open('receipts/genetic_evidence_pool_freeze.json'));assert sha(pool_path)==freeze['pool_sha256'];assert len(pool)==40
fr=json.load(open('receipts/discovery_rank_freeze.json'))
assert sha('receipts/discovery_rank_freeze.json')==freeze['freeze_receipt_sha256']
for m in ['B0','G1']:
 p=fr['products'][m+'_top20'];assert sha(p)==fr['artifacts'][p]
ref=json.load(open('config/gene_mapping_reference.json'));assert sha(ref['path'])==ref['sha256']
hg=list(csv.DictReader(open(ref['path']),delimiter='\t'));approved={r['symbol']:r for r in hg if r['status']=='Approved'}
lookup=defaultdict(set)
for r in hg:
 if r['status']!='Approved':continue
 for name in [r['symbol']]+r['prev_symbol'].split('|'):
  if name:lookup[name].add(r['ensembl_gene_id'])
ids={r['gene_id'] for r in pool}; names={k:next(iter(v)) for k,v in lookup.items() if len(v)==1 and next(iter(v)) in ids}
for r in pool:assert approved[r['approved_symbol']]['ensembl_gene_id']==r['gene_id']
source=json.load(open('data/genetics/supplement_download.json'));assert sha(source['path'])==source['sha256']
w=Workbook(source['path']);hits=[];invent=[];excluded=[]
# Exact column assignments transcribed from actual worksheet headers.
# (column, literal header, category); generic prose never creates gene evidence.
cols={
 11:[('J','Nearest gene','nearest_gene_with_consequence')],
 12:[('L','Gene','coding_variant_annotation')],
 13:[('L','Genes within Chromatin Interactiond','chromatin_contact')],
 14:[('C','Gene','published_SMR_eQTL')],
 15:[('A','Gene','published_SMR_eQTL'),('K','Genes with Chromatin Interaction','chromatin_contact')],
 16:[('A','Gene','published_SMR_eQTL'),('J','Genes with Chromatin Interaction','chromatin_contact')],
 17:[('D','Closest Gene','nearest_gene_to_mQTL'),('K','Genes with Chromatin Interaction','chromatin_contact')],
 18:[('D','Closest Gene','nearest_gene_to_mQTL'),('K','Genes with Chromatin Interaction','chromatin_contact')],
 31:[('K','Gene ID','eQTL_MTAG_shared_locus')],32:[('L','Gene ID','eQTL_MTAG_shared_locus')]}
qstats={gid:dict(rows=0,min_p=None,min_row=None,top_rows=[]) for gid in ids}
for sheet in w.sheets:
 num=int(re.match(r'Supp(\d+)',sheet)[1]);nrows=0;nonempty=0;headers={};count_before=len(hits)
 if num not in set(cols)|{3,10,34,35}:continue
 for rn,cells in w.rows(sheet):
  nrows+=1;nonempty+=bool(cells);v=vals(cells)
  if rn<=3:headers.update(cells)
  if num in (34,35):
   gid=v.get('A','').split('.')[0]
   if rn<=3 or gid not in ids:continue
   st=qstats[gid];st['rows']+=1;pval=float(v['L']);entry=dict(sheet=sheet,row=rn,cells=cells)
   if st['min_p'] is None or pval<st['min_p']:st['min_p']=pval;st['min_row']=entry
   if v.get('N')=='1':st['top_rows'].append(entry)
   continue
  for col,header,cat in cols.get(num,[]):
   hr=3 if num==13 else 2
   if rn==hr:assert v.get(col)==header,(sheet,col,v.get(col),header)
   if rn<=hr:continue
   for token in re.split(r'[,;\s]+',v.get(col,'')):
    token=token.strip('*%^')
    if token in names:
     hits.append(dict(gene_id=names[token],source_symbol=token,source_sheet=sheet,source_row=rn,source_cell=col+str(rn),source_column=header,evidence_category=cat,source_cells_json=json.dumps(cells,ensure_ascii=False,sort_keys=True),source_sha256=source['sha256'],publication='10.1038/s41588-023-01323-z'))
  for col in ['L'] if num==31 else ['M'] if num==32 else []:
   if rn<=2:continue
   for name,gid in names.items():
    if re.search(r'(?<![A-Za-z0-9_])'+re.escape(name)+r'(?![A-Za-z0-9_])',v.get(col,'')):
     excluded.append(dict(gene_id=gid,source_sheet=sheet,source_cell=col+str(rn),reason='Gene Function prose mentions this gene but gene-ID column names a different target',actual_mapped_gene=v.get('K' if num==31 else 'L','')))
 invent.append(dict(sheet=sheet,rows_including_empty=nrows,nonempty_rows=nonempty,pool_evidence_rows=len(hits)-count_before,headers_json=json.dumps(headers,ensure_ascii=False,sort_keys=True)))
 print('AUDITED',sheet,nrows,flush=True)
write_tsv('reports/genetic_evidence_rows.tsv',hits)
write_tsv('reports/genetic_evidence_table_inventory.tsv',invent)
write_tsv('reports/genetic_evidence_rejected_prose_hits.tsv',excluded,['gene_id','source_sheet','source_cell','reason','actual_mapped_gene'])
qout=[]
for gid in sorted(ids):
 st=qstats[gid]
 qout.append(dict(gene_id=gid,tested_variant_rows=st['rows'],minimum_nominal_p=st['min_p'] if st['min_p'] is not None else '',minimum_row_json=json.dumps(st['min_row'],sort_keys=True) if st['min_row'] else '',author_top_variant_rows=len(st['top_rows']),author_top_rows_json=json.dumps(st['top_rows'],sort_keys=True),interpretation='cis_QTL_context_only_not_disease_colocalization' if st['rows'] else 'not_present_in_restricted_UK_QTL_table_not_negative'))
write_tsv('reports/genetic_evidence_uk_qtl_context.tsv',qout)
summary=[]
for r in pool:
 hh=[h for h in hits if h['gene_id']==r['gene_id']];st=qstats[r['gene_id']]
 summary.append({**r,'genetic_status':'supporting_mapping_limited' if hh else 'no_available_mapping_in_audited_source','mapping_categories':';'.join(sorted({h['evidence_category'] for h in hh})),'evidence_rows':len(hh),'publication_count':1 if hh else 0,'uk_qtl_context_rows':st['rows'],'conflicting_evidence':'not_established','causal_direction':'unknown','stromal_specificity':'not_established','participant_independence':'no_expression_cohort_reuse_documented_individual_nonoverlap_not_verified','limitation':'same_publication_correlated_annotations;QTL_context_not_disease_support;absence_not_negative'})
write_tsv('reports/genetic_evidence_summary.tsv',summary)
receipt=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_pool=40,n_supporting_mapping=sum(r['genetic_status']=='supporting_mapping_limited' for r in summary),n_no_available_mapping=sum(r['genetic_status']=='no_available_mapping_in_audited_source' for r in summary),n_conflicting=0,source=source,pool_sha256=sha(pool_path),hgnc_sha256=ref['sha256'],elapsed_seconds=time.monotonic()-start,external_effects_read=False,rankings_modified=False,de_novo_inference=False,validation='exact source headers and frozen input hashes asserted; 40/40 summary coverage',output_hashes={str(p):sha(p) for p in Path('reports').glob('genetic_evidence*.tsv')})
Path('receipts/genetic_evidence_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
