"""Pre-effect count preparation only; original lexical fractional counts retained."""
import argparse,collections,csv,decimal,gzip,hashlib,io,json,pathlib,re
ROOT=pathlib.Path(__file__).resolve().parents[1]
OSF=ROOT/'data/input_routes/osf/selected/04_Analysis/S3_COMPLETION_EXTENSION/GSE266265_DONOR_RECONSTRUCTION'
decimal.getcontext().prec=50
ENSG=re.compile(r'^(ENSG[0-9]+)(?:\.([0-9]+))?$')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def table(path,delimiter='\t'):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))
def count_value(text):
    x=decimal.Decimal(text)
    if not x.is_finite() or x<0:raise ValueError('Counts must be finite nonnegative original estimates')
    return x
def restore_count_rows(mapping,counts,donors):
    groups=collections.defaultdict(list)
    for r in mapping:groups[r['canonical_gene']].append(r)
    if len({r['source_gene'] for r in mapping})!=len(mapping):raise ValueError('Duplicate original source feature ID')
    by={r['gene']:r for r in counts}
    if len(by)!=len(counts) or set(by)!=set(groups):raise ValueError('Count genes differ from declared canonical mapping')
    restored=[]
    for r in sorted(mapping,key=lambda r:int(r['source_row'])):
        if len(groups[r['canonical_gene']])!=1:continue
        observed=by[r['canonical_gene']]
        for donor in donors:count_value(observed[donor])
        restored.append({'source_feature_id':r['source_gene'],**{d:observed[d] for d in donors}})
    dropped=[r for r in counts if len(groups[r['gene']])>1]
    return restored,dropped

def map_features(features,annotation):
    approved=collections.defaultdict(list);previous=collections.defaultdict(list);stable=collections.defaultdict(list)
    for r in annotation:
        if r['status']!='Approved':continue
        approved[r['symbol']].append(r)
        for old in filter(None,r.get('prev_symbol','').split('|')):previous[old].append(r)
        if ENSG.fullmatch(r['ensembl_gene_id']):stable[r['ensembl_gene_id']].append(r)
    results=[]
    for feature in features:
        hit=ENSG.fullmatch(feature);version=hit.group(2) or '' if hit else ''
        if hit:records=stable.get(hit.group(1),[]);method='source_Ensembl_version_stripped_HGNC_verified'
        elif feature in approved:records=approved[feature];method='exact_HGNC_approved_symbol'
        else:records=previous.get(feature,[]);method='exact_unique_HGNC_previous_symbol'
        hgnc_ids={r['hgnc_id'] for r in records};targets={r['ensembl_gene_id'] for r in records if ENSG.fullmatch(r['ensembl_gene_id'])}
        decision='include';reason='one official approved gene and one stable Ensembl target'
        if not records:decision='exclude_unmapped';reason='not found in pinned HGNC named-gene universe'
        elif len(hgnc_ids)!=1 or len(targets)>1:decision='exclude_ambiguous_symbol';reason='multiple official gene records or stable targets'
        elif len(targets)!=1:decision='exclude_no_stable_target';reason='HGNC approved gene lacks one valid stable Ensembl target'
        target=next(iter(targets)) if len(targets)==1 else ''
        results.append({'source_feature_id':feature,'source_version':version,'stable_gene_id':target,'approved_symbol':'|'.join(sorted({r['symbol'] for r in records})),'hgnc_id':'|'.join(sorted(hgnc_ids)),'mapping_method':method,'decision':decision,'reason':reason})
    groups=collections.defaultdict(list)
    for r in results:
        if r['decision']=='include':groups[r['stable_gene_id']].append(r)
    for group in groups.values():
        if len(group)>1:
            for r in group:r['decision']='exclude_duplicate_target';r['reason']='entire conflicting source-to-target group excluded; no sum or maximum selection'
    return results

def write_table(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.partial')
    if str(path).endswith('.gz'):
        raw=tmp.open('wb');gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0);f=io.TextIOWrapper(gz,encoding='utf-8',newline='')
    else:raw=None;f=tmp.open('w',newline='')
    try:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    finally:
        f.close()
        if raw:raw.close()
    tmp.replace(path)
def totals(rows,donors):return {d:sum((count_value(r[d]) for r in rows),decimal.Decimal(0)) for d in donors}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='data/processed/discovery');a=ap.parse_args();out=ROOT/a.output
    reference=json.loads((ROOT/'config/gene_mapping_reference.json').read_text());annotation_path=ROOT/reference['path']
    if sha(annotation_path)!=reference['sha256']:raise ValueError('Frozen HGNC annotation changed')
    mapping=table(OSF/'SOURCE_GENE_MAPPING.csv',',');counts=table(OSF/'pseudobulk/P077.csv.gz',',')
    receipt=json.loads((ROOT/'receipts/osf_route_ranges.json').read_text());expected=next(r['sha256'] for r in receipt['requests'] if r.get('member','').endswith('/pseudobulk/P077.csv.gz'))
    if sha(OSF/'pseudobulk/P077.csv.gz')!=expected:raise ValueError('Original P077 bytes changed')
    clinical=table(ROOT/'manifests/discovery_clinical_candidates.tsv');byclinical={r['sample_id']:r for r in clinical}
    flagged={'sample81','sample108','sample115','sample119','sample121'}
    primary=[r['sample_id'] for r in clinical if r['author_proliferative_candidate']=='True' and r['sample_id'] not in flagged|{'sample73'}]
    strict=[d for d in primary if byclinical[d]['phase_histology_conflict']!='True']
    if len(primary)!=22 or len(strict)!=20:raise ValueError('Clinical freeze cohort size changed')
    final_clinical_path=ROOT/'manifests/final_discovery_clinical_selected.tsv'
    final_clinical={r['donor_id']:r for r in table(final_clinical_path)}
    if len(final_clinical)!=22 or set(final_clinical)!=set(primary):raise ValueError('Final clinical ledger differs from primary selection')
    for d in primary:
        c=byclinical[d];final=final_clinical[d]
        if final['disease']!=c['endometriosis status'] or decimal.Decimal(final['age'])!=decimal.Decimal(c['age_years']) or final['capture_batch']!=c['10x capturing batch'] or final['adenomyosis']!=c['Adenomyosis'] or (final['strict_histology_sensitivity']=='True')!=(d in strict):
            raise ValueError('Final clinical ledger metadata conflict '+d)
    source_donors=[k for k in counts[0] if k!='gene']
    if set(source_donors)!=set(primary)|{'sample73'}:raise ValueError('Compact original donor pool changed')
    if collections.Counter(byclinical[d]['endometriosis status'] for d in primary)!={'ENDO':12,'non-ENDO':10}:raise ValueError('Primary disease groups changed')
    for fname in ['SOURCE_GENE_MAPPING.csv','PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz']:
        pinned=next(r['sha256'] for r in receipt['requests'] if r.get('member','').endswith('/'+fname))
        if sha(OSF/fname)!=pinned:raise ValueError('Original OSF member hash changed '+fname)
    lineage=json.loads((ROOT/'config/discovery_lineage.json').read_text());included=set(lineage['include_exact_refined_labels']);mural=set(lineage['exclude_exact_mural_labels']);cellcounts=collections.Counter();nn=0
    with gzip.open(OSF/'PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz','rt') as f:
        for cell in csv.DictReader(f):
            nn+=1
            if cell['AnnotationRefined'] in included:
                if cell['AnnotationMain']!='stromal':raise ValueError('Fibroblast refined label conflicts with main stroma')
                if cell['AnnotationRefined'] in mural:raise ValueError('Mural label enters lineage')
                cellcounts[cell['sample']]+=1
    author={r['sample']:int(r['cells']) for r in table(OSF/'DONOR_CELLTYPE_MANIFEST.csv',',') if r['population_id']=='P077'}
    if dict(cellcounts)!=author or any(cellcounts[d]<50 for d in primary):raise ValueError('Fixed own lineage disagrees with P077 membership/coverage')
    restored,many=restore_count_rows(mapping,counts,primary)
    if len(restored)!=34550 or len(many)!=50:raise ValueError('Expected singleton conservation changed')
    features=map_features([r['source_feature_id'] for r in restored],table(annotation_path));fby={r['source_feature_id']:r for r in features}
    accepted=[r for r in restored if fby[r['source_feature_id']]['decision']=='include'];rejected=[r for r in restored if fby[r['source_feature_id']]['decision']!='include']
    stable_rows=[{'gene_id':fby[r['source_feature_id']]['stable_gene_id'],**{d:r[d] for d in primary}} for r in accepted];stable_rows.sort(key=lambda r:r['gene_id'])
    if len({r['gene_id'] for r in stable_rows})!=len(stable_rows):raise ValueError('Duplicate final gene ID')
    allsum=totals(counts,primary);singletonsum=totals(restored,primary);manysum=totals(many,primary);mapped_sum=totals(accepted,primary);rejected_sum=totals(rejected,primary)
    for d in primary:
        if allsum[d]!=singletonsum[d]+manysum[d] or singletonsum[d]!=mapped_sum[d]+rejected_sum[d]:raise ValueError('Exact decimal count mass not conserved')
    write_table(out/'source_counts.tsv.gz',restored,['source_feature_id']+primary);write_table(out/'counts_stable.tsv.gz',stable_rows,['gene_id']+primary)
    for name,rows,first in [('source_counts.tsv.gz',restored,'source_feature_id'),('counts_stable.tsv.gz',stable_rows,'gene_id')]:
        replay=table(out/name)
        if replay!=rows:raise ValueError(f'Lexical count/ID roundtrip failed {name}')
    write_table(ROOT/'manifests/discovery_features_mapping.tsv',features,list(features[0]))
    original_groups=collections.defaultdict(list)
    for r in mapping:original_groups[r['canonical_gene']].append(r)
    excluded=[{'canonical_gene':g,'original_source_features':'|'.join(r['source_gene'] for r in rows),'source_row_count':len(rows),'reason':'irreversible_multi_source_group'} for g,rows in original_groups.items() if len(rows)>1]
    write_table(ROOT/'manifests/discovery_features_irreversible_exclusions.tsv',excluded,list(excluded[0]))
    ledger=[]
    for d in primary:
        c=byclinical[d];ledger.append({'donor_id':d,'disease':c['endometriosis status'],'age':final_clinical[d]['age'],'batch':final_clinical[d]['capture_batch'],'adenomyosis':final_clinical[d]['adenomyosis'],'GSM':c['GSM'],'primary_included':True,'strict_phase_sensitivity_included':d in strict,'phase_histology_conflict':c['phase_histology_conflict'],'phase_basis':'author_adjudicated_proliferative','hormone_eligibility_basis':c['hormone_interpretation'],'qualifying_fibroblast_cells':cellcounts[d],'original_P077_estimated_UMIs':str(allsum[d]),'irreversible_group_excluded_UMIs':str(manysum[d]),'restored_singleton_UMIs':str(singletonsum[d]),'mapping_excluded_UMIs':str(rejected_sum[d]),'mapped_retained_UMIs':str(mapped_sum[d]),'mapping_excluded_fraction_of_singletons':str(rejected_sum[d]/singletonsum[d]),'overall_excluded_fraction':str((allsum[d]-mapped_sum[d])/allsum[d]),'raw_count_rounding':'none','normalization':'none'})
    write_table(out/'selected_donor_count_ledger.tsv',ledger,list(ledger[0]));write_table(ROOT/'manifests/discovery_features_donor_count_ledger.tsv',ledger,list(ledger[0]))
    outputs=[out/'source_counts.tsv.gz',out/'counts_stable.tsv.gz',out/'selected_donor_count_ledger.tsv',ROOT/'manifests/discovery_features_mapping.tsv',ROOT/'manifests/discovery_features_irreversible_exclusions.tsv']
    total_source_features=len(mapping);decisions=dict(collections.Counter(r['decision'] for r in features));methods=dict(collections.Counter(r['mapping_method'] for r in features if r['decision']=='include'))
    summary={'status':'COUNTS_PREPARED_NO_DE_M0_STILL_OPEN','input_source':'public OSF P077 unrounded fractional RNA Alevin donor sums, inherited author QC and annotation','original_count_sha256':expected,'original_RDS_reconstructed_this_round':False,'source_features':total_source_features,'original_canonical_rows':len(counts),'irreversible_groups_excluded':len(many),'source_rows_in_irreversible_groups':sum(len(v) for v in original_groups.values() if len(v)>1),'singleton_source_rows':len(restored),'versioned_Ensembl_singletons':sum(bool(ENSG.fullmatch(r['source_feature_id'])) for r in restored),'primary_donors':primary,'primary_groups':{'ENDO':12,'non-ENDO':10},'strict_sensitivity_donors':strict,'strict_groups':{'ENDO':11,'non-ENDO':9},'annotation_rows':reference['rows'],'annotation_sha256':reference['sha256'],'annotation_generation':reference['object_generation'],'mapped_stable_rows':len(stable_rows),'mapping_decisions':decisions,'accepted_mapping_methods':methods,'original_source_IDs_preserved':True,'exact_decimal_mass_conservation':True,'final_clinical_ledger_exact_match':True,'original_lexical_count_roundtrip':True,'fractional_estimates_not_rounded':True,'normalized':False,'source_primary_lineage_cells':sum(cellcounts[d] for d in primary),'minimum_primary_lineage_cells':min(cellcounts[d] for d in primary),'annotation_limit':'Pinned HGNC approved named-gene universe, not full Ensembl feature annotation; unknown and ambiguous features explicitly excluded','lineage_limit':lineage['independent_marker_expression_audit'],'clinical_source_limit':'archived primary clinical supplement publisher byte equivalence not independently verified','disease_effects_computed':False,'network_gene_universe_frozen':False,'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [OSF/'SOURCE_GENE_MAPPING.csv',OSF/'PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz',ROOT/'manifests/discovery_clinical_candidates.tsv',ROOT/'docs/PRE_EFFECT_METHOD_DECISIONS.md',ROOT/'config/discovery_lineage.json',ROOT/'config/gene_mapping_reference.json',ROOT/'config/gene_mapping_policy.json',final_clinical_path]},'outputs':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in outputs]}
    (out/'provenance.json').write_text(json.dumps(summary,indent=2)+'\n');(ROOT/'receipts/discovery_counts_preparation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
