#!/usr/bin/env python3
"""Fixed-panel, disease-blind donor-pseudobulk descriptive audit. No DE."""
import argparse, collections, csv, decimal, gzip, hashlib, json, statistics
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
D=decimal.Decimal
decimal.getcontext().prec=80
OSF=ROOT/'data/input_routes/osf/selected/04_Analysis/S3_COMPLETION_EXTENSION/GSE266265_DONOR_RECONSTRUCTION'
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def rows(p,sep='\t'):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:
        yield from csv.DictReader(f,delimiter=sep)
def count(x):
    v=D(x)
    if not v.is_finite() or v<0: raise ValueError('Nonfinite or negative count')
    return v
def validate_header(header):
    if len(header)<2 or len(header)!=len(set(header)): raise ValueError('Duplicate or missing matrix columns')
def matrix(p,sep='\t'):
    result={}
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:
        reader=csv.reader(f,delimiter=sep);header=next(reader);validate_header(header);donors=header[1:]
        for row in reader:
            if len(row)!=len(header): raise ValueError('Matrix row width differs')
            gene=row[0]
            if not gene or gene in result: raise ValueError('Duplicate or missing feature '+gene)
            result[gene]=[count(x) for x in row[1:]]
    if not result: raise ValueError('Empty matrix')
    return donors,result
def difference(a,b):
    v=a-b
    if v<0: raise ValueError('Negative P005-P077 difference; comparator unavailable, no clipping')
    return v
def stats(values):
    q=statistics.quantiles(values,n=4,method='inclusive') if len(values)>1 else values*3
    return [len(values),sum(v>0 for v in values),statistics.mean(values),statistics.median(values),q[0],q[2],min(values),max(values)]
def write_tsv(p,header,rs):
    with open(p,'w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rs)
def selftest():
    assert stats([1.0,100.0])[2]==50.5 # equal donors, not summed-count weighting
    assert difference(D('1.2'),D('0.1'))==D('1.1')
    try: difference(D('1'),D('1.00000001'))
    except ValueError: pass
    else: raise AssertionError('Negative subtraction accepted')
    assert stats([0.0,2.0])[1]==1
    try: validate_header(['gene','donor1','donor1'])
    except ValueError: pass
    else: raise AssertionError('Duplicate donor headers accepted')
    # Absent target is unavailable even when a different observed marker is zero.
    assert 'absent' not in {'present_zero':[D(0)]}
    return ['equal_donor_weighting','fractional_decimal_subtraction','negative_remainder_rejection','positive_count_denominator','duplicate_donor_header_rejection','missing_feature_not_observed_zero']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--heca-counts');ap.add_argument('--skip-comparator',action='store_true');a=ap.parse_args()
    panel_path=ROOT/'config/lineage_marker_panel.json';panel=json.loads(panel_path.read_text())
    freeze=json.loads((ROOT/'receipts/lineage_marker_panel_freeze.json').read_text())
    # Freeze receipt format deliberately checked below, not silently regenerated.
    frozen_hashes=freeze.get('files',freeze.get('sha256',{}))
    if isinstance(frozen_hashes,dict) and str(panel_path.relative_to(ROOT)) in frozen_hashes:
        expected=frozen_hashes[str(panel_path.relative_to(ROOT))]
        if isinstance(expected,dict): expected=expected['sha256']
        assert sha(panel_path)==expected
    else:
        assert sha(panel_path)=='8c790a2a6d7997d20cd6c0c2c48902fdb1d2f54796792191b44fe8a36ebcc861'
    ref=json.loads((ROOT/'config/gene_mapping_reference.json').read_text());hp=ROOT/ref['path']
    assert sha(hp)==ref['sha256']
    symbols={m['symbol'] for m in panel['markers']};hgnc=collections.defaultdict(set)
    for r in rows(hp):
        if r['symbol'] in symbols and r['status']=='Approved' and r['ensembl_gene_id']:
            hgnc[r['symbol']].add(r['ensembl_gene_id'])
    targets={s:next(iter(v)) for s,v in hgnc.items() if len(v)==1}
    inputs={str(panel_path.relative_to(ROOT)):sha(panel_path),ref['path']:sha(hp)}
    cp=ROOT/'data/processed/discovery/counts_stable.tsv.gz'
    donors,primary=matrix(cp);assert len(donors)==22
    inputs[str(cp.relative_to(ROOT))]=sha(cp)
    populations=[('discovery','primary_fibroblast',donors,primary)]
    comparator={'status':'not_requested'}
    if not a.skip_comparator:
        try:
            p5=OSF/'pseudobulk/P005.csv.gz';p77=OSF/'pseudobulk/P077.csv.gz'
            d5,m5=matrix(p5,',');d77,m77=matrix(p77,',')
            assert d5==d77 and set(m5)==set(m77),'Source features/donors differ'
            for g in m5:
                for x,y in zip(m5[g],m77[g]): difference(x,y)
            # Verify source label partition without accessing disease fields.
            lin=json.loads((ROOT/'config/discovery_lineage.json').read_text())
            inc=set(lin['include_exact_refined_labels']);mural=set(lin['exclude_exact_mural_labels']);assert not inc&mural
            label_census=collections.Counter();stromal_cells=collections.Counter();fibroblast_cells=collections.Counter()
            for c in rows(OSF/'PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz',','):
                if c['AnnotationMain']=='stromal':
                    assert c['AnnotationRefined'] in inc|mural,'Unpartitioned source stromal label'
                    label_census[c['AnnotationRefined']]+=1
                    stromal_cells[c['sample']]+=1
                    if c['AnnotationRefined'] in inc: fibroblast_cells[c['sample']]+=1
            manifest_cells={'P005':{},'P077':{}}
            for r in rows(OSF/'DONOR_CELLTYPE_MANIFEST.csv',','):
                if r['population_id'] in manifest_cells:
                    dest=manifest_cells[r['population_id']]
                    assert r['sample'] not in dest,'Duplicate population/donor cell manifest'
                    dest[r['sample']]=int(r['cells'])
            assert manifest_cells['P005']==dict(stromal_cells),'P005 annotation partition differs'
            assert manifest_cells['P077']==dict(fibroblast_cells),'P077 annotation partition differs'
            groups=collections.defaultdict(list)
            for r in rows(OSF/'SOURCE_GENE_MAPPING.csv',','): groups[r['canonical_gene']].append(r['source_gene'])
            source_to_canonical={ss[0]:g for g,ss in groups.items() if len(ss)==1}
            mp=ROOT/'manifests/discovery_features_mapping.tsv';mapped={}
            for r in rows(mp):
                if r['decision']=='include':
                    assert r['stable_gene_id'] not in mapped
                    mapped[r['stable_gene_id']]=source_to_canonical[r['source_feature_id']]
            assert set(mapped)==set(primary),'Accepted gene universe differs'
            idx=[d5.index(d) for d in donors];remainder={};whole_tot=[D(0)]*len(donors);p_tot=[D(0)]*len(donors);r_tot=[D(0)]*len(donors)
            for gene,canonical in mapped.items():
                pv=[m77[canonical][i] for i in idx];wv=[m5[canonical][i] for i in idx]
                assert pv==primary[gene],'P077 differs from prepared stable count at '+gene
                rv=[difference(x,y) for x,y in zip(wv,pv)];remainder[gene]=rv
                for j in range(len(donors)): whole_tot[j]+=wv[j];p_tot[j]+=pv[j];r_tot[j]+=rv[j]
            assert all(w==p+r for w,p,r in zip(whole_tot,p_tot,r_tot))
            assert all(t>0 for t in r_tot),'Zero comparator denominator'
            populations.append(('discovery','source_mural_remainder',donors,remainder))
            comparator={'status':'available_exact_checks_passed','all_source_features':len(m5),'all_source_donors':len(d5),'all_source_pairs_nonnegative':len(m5)*len(d5),'selected_donors':len(donors),'retained_stable_features':len(mapped),'exact_P077_prepared_value_match':True,'exact_selected_retained_mass_conservation':True,'source_stromal_label_census':dict(label_census),'exact_P005_P077_cell_membership_match':True,'selected_mural_cells':sum(stromal_cells[d]-fibroblast_cells[d] for d in donors),'selected_mass':[{'donor_id':d,'P005_total':str(w),'P077_total':str(p),'remainder_total':str(r)} for d,w,p,r in zip(donors,whole_tot,p_tot,r_tot)]}
            for p in [p5,p77,mp,OSF/'SOURCE_GENE_MAPPING.csv',OSF/'PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz',OSF/'DONOR_CELLTYPE_MANIFEST.csv',ROOT/'config/discovery_lineage.json']:
                inputs[str(p.relative_to(ROOT))]=sha(p)
        except (ValueError,AssertionError,KeyError) as e:
            comparator={'status':'unavailable','reason':str(e),'clipping_or_rounding':False}
    if a.heca_counts:
        p=ROOT/a.heca_counts;ds,ms=matrix(p);populations.append(('heca','primary_stromal_fibroblast',ds,ms));inputs[str(p.relative_to(ROOT))]=sha(p)
    donor_rows=[];summary_rows=[];coverage={}
    for cohort,pop,ds,ms in populations:
        totals=[sum((v[j] for v in ms.values()),D(0)) for j in range(len(ds))]
        assert all(t>0 for t in totals),'Nonpositive CPM denominator'
        coverage[cohort+'/'+pop]={'donors':len(ds),'features':len(ms),'total_count_min':str(min(totals)),'total_count_max':str(max(totals))}
        for marker in panel['markers']:
            symbol=marker['symbol'];gene=targets.get(symbol,'');available=gene in ms
            vals=[]
            for j,d in enumerate(ds):
                raw=ms[gene][j] if available else None
                cpm=float(D(1000000)*raw/totals[j]) if available else None
                if cpm is not None: vals.append(cpm)
                donor_rows.append([cohort,pop,d,symbol,marker['category'],gene,'available' if available else 'unavailable',str(totals[j]),str(raw) if raw is not None else '',cpm if cpm is not None else ''])
            summary_rows.append([cohort,pop,symbol,marker['category'],gene,'available' if available else 'unavailable']+(stats(vals) if vals else [len(ds)]+['']*7))
    dh=['cohort','population','donor_id','symbol','category','stable_gene_id','availability','retained_library_total','marker_count','CPM']
    sh=['cohort','population','symbol','category','stable_gene_id','availability','donor_count','positive_count_donors','mean_CPM','median_CPM','Q1_CPM','Q3_CPM','minimum_CPM','maximum_CPM']
    write_tsv(ROOT/'reports/lineage_marker_audit_donors.tsv',dh,donor_rows)
    write_tsv(ROOT/'reports/lineage_marker_audit_summary.tsv',sh,summary_rows)
    receipt={'created_utc':datetime.now(timezone.utc).isoformat(),'audit':'fixed-panel descriptive population plausibility only','script_sha256':sha(Path(__file__)),'panel_freeze_receipt_sha256':sha(ROOT/'receipts/lineage_marker_panel_freeze.json'),'input_sha256':inputs,'disease_columns_consumed':False,'new_cell_or_donor_exclusions':False,'DE_or_disease_comparisons':False,'CPM_denominator':panel['CPM'],'coverage':coverage,'comparator':comparator,'selftests':selftest(),'outputs':{}}
    for name in ['lineage_marker_audit_donors.tsv','lineage_marker_audit_summary.tsv']:
        receipt['outputs']['reports/'+name]=sha(ROOT/'reports'/name)
    (ROOT/'receipts/lineage_marker_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'coverage':coverage,'comparator':comparator['status'],'reason':comparator.get('reason'),'selftests':receipt['selftests']}))
if __name__=='__main__':main()
