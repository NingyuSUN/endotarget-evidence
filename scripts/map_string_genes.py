#!/usr/bin/env python3
"""Pinned name-based STRING gene background; no expression or disease inputs."""
import collections,csv,datetime,decimal,gzip,hashlib,io,json,pathlib,re,resource,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
ENSG=re.compile(r'^ENSG[0-9]{11}$')
decimal.getcontext().prec=50
D=decimal.Decimal

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def map_proteins(info,annotation):
    approved=collections.defaultdict(list);previous=collections.defaultdict(list)
    for r in annotation:
        if r['status']!='Approved':continue
        approved[r['symbol']].append(r)
        for name in filter(None,r.get('prev_symbol','').split('|')):previous[name].append(r)
    result=[];seen=set()
    for protein in info:
        pid=protein['#string_protein_id'];name=protein['preferred_name']
        if pid in seen:raise ValueError('Duplicate protein ID; cannot overwrite crosswalk')
        if not pid.startswith('9606.'):raise ValueError('Nonhuman protein info ID')
        seen.add(pid)
        if name in approved:records=approved[name];method='exact_HGNC_approved_symbol'
        else:records=previous.get(name,[]);method='exact_unique_HGNC_previous_symbol'
        ids={r['hgnc_id'] for r in records};targets={r['ensembl_gene_id'] for r in records if ENSG.fullmatch(r['ensembl_gene_id'])}
        decision='include';reason='not_excluded'
        if not records:decision='exclude_unmapped';reason='source preferred_name absent from official approved/previous symbols; no guess'
        elif len(ids)!=1 or len(targets)>1:decision='exclude_ambiguous';reason='multiple approved HGNC records or stable targets for source name'
        elif len(targets)!=1:decision='exclude_no_stable_target';reason='approved HGNC record lacks one valid unversioned stable Ensembl gene ID'
        result.append({'protein_id':pid,'source_preferred_name':name,'stable_gene_id':next(iter(targets)) if len(targets)==1 else '', 'approved_symbol':'|'.join(sorted({r['symbol'] for r in records})),'hgnc_id':'|'.join(sorted(ids)),'mapping_method':method,'decision':decision,'excluded_reason':reason})
    return result

def collapse_gene_edges(edges,crosswalk):
    by={r['protein_id']:r for r in crosswalk}
    if len(by)!=len(crosswalk):raise ValueError('Duplicate crosswalk protein IDs')
    graph={};reasons=collections.Counter();reason_mass=collections.defaultdict(lambda:D(0))
    stats={'source_protein_edges':0,'endpoint_excluded_edges':0,'gene_self_loops_excluded':0,'eligible_nonloop_protein_edges':0,'source_weight_sum':D(0),'endpoint_excluded_weight_sum':D(0),'gene_self_loop_weight_sum':D(0),'eligible_nonloop_weight_sum':D(0)}
    for edge in edges:
        a,b=edge['protein1'],edge['protein2'];weight=D(edge['weight'])
        if not weight.is_finite() or weight<D('0.700') or weight>D('1.000'):raise ValueError('Invalid pinned direct protein edge weight')
        if a not in by or b not in by:raise ValueError('Graph endpoint absent from source protein info crosswalk')
        stats['source_protein_edges']+=1;stats['source_weight_sum']+=weight
        excluded={by[p]['decision'] for p in [a,b] if by[p]['decision']!='include'}
        if excluded:
            stats['endpoint_excluded_edges']+=1;stats['endpoint_excluded_weight_sum']+=weight
            reason='|'.join(sorted(excluded));reasons[reason]+=1;reason_mass[reason]+=weight;continue
        ga,gb=by[a]['stable_gene_id'],by[b]['stable_gene_id']
        if not ENSG.fullmatch(ga) or not ENSG.fullmatch(gb):raise ValueError('Included endpoint lacks valid stable gene ID')
        if ga==gb:
            stats['gene_self_loops_excluded']+=1;stats['gene_self_loop_weight_sum']+=weight;continue
        key=(ga,gb) if ga<gb else (gb,ga)
        graph[key]=max(weight,graph.get(key,D(0)))
        stats['eligible_nonloop_protein_edges']+=1;stats['eligible_nonloop_weight_sum']+=weight
    stats['output_gene_edges']=len(graph);stats['output_gene_nodes']=len({g for pair in graph for g in pair})
    stats['output_gene_weight_sum']=sum(graph.values(),D(0))
    stats['duplicate_gene_edges_collapsed']=stats['eligible_nonloop_protein_edges']-len(graph)
    stats['duplicate_collapse_weight_loss']=stats['eligible_nonloop_weight_sum']-stats['output_gene_weight_sum']
    stats['endpoint_exclusion_reasons']=dict(reasons);stats['endpoint_exclusion_reason_weight_sums']=dict(reason_mass)
    if stats['source_weight_sum']!=stats['endpoint_excluded_weight_sum']+stats['gene_self_loop_weight_sum']+stats['duplicate_collapse_weight_loss']+stats['output_gene_weight_sum']:
        raise ValueError('Exact edge-weight conservation failed')
    if stats['source_protein_edges']!=stats['endpoint_excluded_edges']+stats['gene_self_loops_excluded']+stats['duplicate_gene_edges_collapsed']+stats['output_gene_edges']:
        raise ValueError('Exact edge-count conservation failed')
    return graph,stats

def write_tsv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.partial')
    if tmp.exists():raise ValueError('Existing partial output preserved; inspect before replay')
    if str(path).endswith('.gz'):
        raw=tmp.open('xb');gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0);f=io.TextIOWrapper(gz,encoding='utf-8',newline='')
    else:raw=None;f=tmp.open('x',newline='')
    try:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    finally:
        f.close()
        if raw:raw.close()
    if path.exists():
        if digest(path)!=digest(tmp):raise ValueError('Previously generated output differs; original and partial preserved')
        tmp.unlink()
    else:tmp.rename(path)

def validate_policy(cfg):
    expected={'STRING_version':'12.5','species_taxid':9606,'HGNC_snapshot':'2026-10-06','source_name_policy':'exact_HGNC_approved_symbol_then_unique_official_previous_symbol','symbol_case_policy':'exact_case_sensitive_no_case_guess','generic_alias_lookup':False,'sequence_equivalence_claim':False,'unique_mapping_requirement':'one_Approved_HGNC_record_and_one_stable_ENSG_target','stable_target_pattern':'^ENSG[0-9]{11}$','multiple_proteins_per_gene':'allowed_name_based_many_to_one','unmapped_or_ambiguous_endpoint_policy':'exclude_edge_and_record_reason','unknown_protein_endpoint_policy':'fail_source_integrity','gene_edge_duplicate_policy':'maximum_original_protein_edge_weight','gene_self_loop_policy':'exclude_after_mapping','weight_scale':'unchanged_original_max_direct_score_divided_by_1000','graph_scope':'all_cached_filtered_protein_reference_no_expression_filter_or_induced_U','disease_input_used':False,'network_downloads':0}
    for k,v in expected.items():
        if cfg.get(k)!=v:raise ValueError('Unsupported frozen gene-mapping policy '+k)

def main():
    start=time.time();config_path=ROOT/'config/network_gene_mapping.json';cfg=json.loads(config_path.read_text());validate_policy(cfg)
    source={}
    for k,pin in cfg['inputs'].items():
        p=ROOT/pin['path']
        if digest(p)!=pin['sha256']:raise ValueError('Pinned input changed '+k)
        source[k]=p
    with gzip.open(source['protein_info'],'rt') as f:info=list(csv.DictReader(f,delimiter='\t'))
    with source['HGNC'].open() as f:annotation=list(csv.DictReader(f,delimiter='\t'))
    crosswalk=map_proteins(info,annotation)
    with gzip.open(source['protein_edges'],'rt') as f:graph,stats=collapse_gene_edges(csv.DictReader(f,delimiter='\t'),crosswalk)
    out=ROOT/'data/reference/string/direct_evidence_gene_edges.tsv.gz'
    edge_rows=[{'gene1':a,'gene2':b,'weight':format(w,'.3f')} for (a,b),w in sorted(graph.items())]
    write_tsv(out,edge_rows,['gene1','gene2','weight'])
    crosswalk_path=ROOT/'manifests/string_protein_gene_crosswalk.tsv';write_tsv(crosswalk_path,crosswalk,list(crosswalk[0]))
    included=[r for r in crosswalk if r['decision']=='include'];gene_proteins=collections.Counter(r['stable_gene_id'] for r in included)
    mapping_stats={'source_protein_info_rows':len(info),'mapping_decisions':dict(collections.Counter(r['decision'] for r in crosswalk)),'included_proteins':len(included),'unique_mapped_genes':len(gene_proteins),'genes_with_multiple_proteins':sum(v>1 for v in gene_proteins.values()),'maximum_proteins_per_gene':max(gene_proteins.values(),default=0),'included_mapping_methods':dict(collections.Counter(r['mapping_method'] for r in included))}
    stats_path=ROOT/'manifests/string_gene_mapping_statistics.json';statistics={'mapping':mapping_stats,'edges':stats,'weight_sum_units':'sum of operational direct evidence weights, not count, probability, independent evidence or biological effect','exact_edge_count_and_weight_conservation':True}
    stats_path.write_text(json.dumps(statistics,indent=2,default=str)+'\n')
    for k,pin in cfg['inputs'].items():
        if digest(source[k])!=pin['sha256']:raise ValueError('Original input changed during mapping '+k)
    with gzip.open(out,'rt') as f:
        replay=list(csv.DictReader(f,delimiter='\t'))
    if replay!=edge_rows:raise ValueError('Gene-edge row/weight roundtrip mismatch')
    report={'status':'NAME_BASED_GENE_REFERENCE_PREPARED_NOT_EXPRESSION_UNIVERSE_OR_RANKING','requested_model_config':cfg['requested_model_config'],'config_frozen_at_utc':cfg['frozen_at_utc'],'run_finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'configuration_path':str(config_path.relative_to(ROOT)),'configuration_sha256':digest(config_path),'code_sha256':digest(pathlib.Path(__file__)),'inputs':cfg['inputs'],'original_protein_graph_preserved':True,'outputs':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in [out,crosswalk_path,stats_path]],'statistics':statistics,'runtime_wall_seconds':time.time()-start,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),'network_downloads':0,'disease_inputs_used':False,'expression_inputs_used':False,'induced_gene_universe_frozen':False,'sequence_or_genome_release_equivalence_verified':False,'mapping_limit':cfg['mapping_limitation']}
    (ROOT/'receipts/string_gene_mapping.json').write_text(json.dumps(report,indent=2,default=str)+'\n');print(json.dumps(report,indent=2,default=str))
if __name__=='__main__':main()
