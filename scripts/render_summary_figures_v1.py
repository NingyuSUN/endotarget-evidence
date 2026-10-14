#!/usr/bin/env python3
"""Render existing donor/evidence tables only; never fit models or change inputs."""
from pathlib import Path
import csv, hashlib, html, json, math
from collections import Counter
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parents[1]
INPUTS = ['reports/stromal_composition_donors.tsv','reports/stromal_composition_denominators.tsv','reports/replication_v2_dossier_pool.tsv','manifests/final_discovery_clinical_selected.tsv','manifests/final_heca_clinical_selected.tsv','manifests/discovery_clinical_candidates.tsv','manifests/final_heca_clinical_donors.tsv','work/statistical/heca_primary/evaluation/null_sets.tsv','work/statistical/heca_primary/evaluation/null_summary.json']
def read(rel):
    with (ROOT / rel).open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def esc(s): return html.escape(str(s), quote=True)
def text(x,y,s,size=14,fill='#243247',extra=''):
    return f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" {extra}>{esc(s)}</text>'
def rect(x,y,w,h,fill,extra=''):
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" {extra}/>'
def line(x1,y1,x2,y2,color='#dce2e9',width=1):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>'
def start(w,h,title,desc):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-labelledby="title desc">',f'<title id="title">{esc(title)}</title><desc id="desc">{esc(desc)}</desc>','<rect width="100%" height="100%" fill="white"/><g font-family="DejaVu Sans, Arial, sans-serif">']
def save(parts,rel):
    p=ROOT/rel; p.parent.mkdir(parents=True,exist_ok=True); p.write_text('\n'.join(parts+['</g></svg>'])+'\n')
# Table validation is deliberately separate from any disease hypothesis test.
comp=read(INPUTS[0]); den=read(INPUTS[1]); pool=read(INPUTS[2]); checks={}
assert len({r['gene_id'] for r in pool}) == len(pool) == 40
assert sum(r['genetic_genetic_status']=='supporting_mapping_limited' for r in pool)==1
assert all(r['genetic_causal_direction']=='unknown' for r in pool)
assert all(r['genetic_conflicting_evidence']=='not_established' for r in pool)
assert all(float(r['external_q'])>=.05 for r in pool if r['external_q'] not in ('','NA'))
for cohort,group_n,total in [('discovery',{'ENDO':12,'non-ENDO':10},150788),('heca',{'ENDO':15,'non-ENDO':7},52744)]:
    clin=read(f'manifests/final_{cohort}_clinical_selected.tsv'); d=[r for r in den if r['cohort']==cohort]; cr=[r for r in comp if r['cohort']==cohort]
    assert len(d)==len(clin)==22 and len({r['donor_id'] for r in d})==22
    assert Counter(r['disease'] for r in clin)==group_n
    assert {r['donor_id'] for r in d}=={r['donor_id'] for r in clin}
    assert sum(int(r['primary_cells']) for r in d)==total
    assert len({(r['donor_id'],r['source_subtype_label']) for r in cr})==len(cr)
    for r in d:
        rs=[z for z in cr if z['donor_id']==r['donor_id']]
        assert sum(int(z['cells']) for z in rs)==int(r['primary_cells'])
        assert all(int(z['primary_cells'])==int(r['primary_cells']) for z in rs)
        assert math.isclose(sum(float(z['proportion_primary']) for z in rs),1,abs_tol=1e-12)
        assert all(math.isclose(float(z['proportion_primary']),int(z['cells'])/int(z['primary_cells']),abs_tol=1e-12) for z in rs)
    checks[cohort]={'donors':22,'disease_counts':dict(Counter(r['disease'] for r in clin)),'primary_cells_or_nuclei':total,'author_subtypes':len({r['source_subtype_label'] for r in cr})}
# Figure 2: entire author-label composition, denominator size, age and literal batch.
parts=start(1640,1200,'Selected donor composition, measurement depth, age and batch','Two separately labeled cohort panels; one column per donor. Author state labels are not harmonized across cohorts. Counts are primary cells or nuclei, not read depth or purity.')
parts += [text(35,42,'Selected donors: composition, primary-lineage counts, age and batch',24),text(35,70,'Descriptive audit only. Each column is one person; source states are retained without cross-cohort harmonization.',16),text(35,95,'Author QC is inherited. Cell/nucleus count is measurement depth, not purity or independent patient replication.',15)]
palette=['#4477AA','#EE7733','#228833','#CC3311','#66CCEE','#AA3377','#BBBBBB','#000000','#DDCC77','#009988','#EE99AA','#332288']
assert len(set(palette)) == 12
for ci,cohort in enumerate(['discovery','heca']):
    ox=35+ci*810; x0=ox+100; step=27; top=192; base=452; cbase=625; abase=771
    clin=sorted(read(f'manifests/final_{cohort}_clinical_selected.tsv'),key=lambda r:(r['disease']=='ENDO',r['donor_id']))
    cr=[r for r in comp if r['cohort']==cohort]; labels=list(dict.fromkeys(r['source_subtype_label'] for r in cr)); cmap={lab:palette[i] for i,lab in enumerate(labels)}
    unit='cells' if cohort=='discovery' else 'nuclei'; title='Discovery: 12 cases / 10 controls' if ci==0 else 'HECA: 15 cases / 7 controls'
    parts += [text(ox,142,title,21),text(ox,167,f"22 donors; {checks[cohort]['primary_cells_or_nuclei']:,} primary {unit}; {len(labels)} author labels",14)]
    for pct in [0,25,50,75,100]:
        y=base-pct/100*260;parts += [line(x0,y,x0+22*step,y),text(x0-12,y+5,f'{pct}%',12,extra='text-anchor="end"')]
    parts += [text(ox,185,'Within-donor composition (%)',12)]
    maxcount=16000 if ci==0 else 10000
    for value in [0,maxcount//2,maxcount]:
        y=cbase-value/maxcount*95;parts += [line(x0,y,x0+22*step,y),text(x0-12,y+4,f'{value:,}',11,extra='text-anchor="end"')]
    parts += [text(ox,503,f'Primary {unit} per donor',12)]
    for age in [15,30,45,50]:
        y=abase-(age-15)/35*95;parts += [line(x0,y,x0+22*step,y),text(x0-12,y+4,str(age),11,extra='text-anchor="end"')]
    parts += [text(ox,675,'Age (years)',12),text(ox,815,'Batch',12),text(ox,905,'Donor ID',12)]
    for i,r in enumerate(clin):
        x=x0+i*step+3; cx=x+10; rs={z['source_subtype_label']:z for z in cr if z['donor_id']==r['donor_id']}; y=base; color='#c45e3b' if r['disease']=='ENDO' else '#325f87'
        for label in labels:
            rr=rs[label]; h=float(rr['proportion_primary'])*260; y-=h
            if h:
                parts.append(f'<rect x="{x}" y="{y}" width="21" height="{h}" fill="{cmap[label]}"><title>{esc(r["donor_id"])}; {esc(label)}; {rr["cells"]}/{rr["primary_cells"]} ({100*float(rr["proportion_primary"]):.2f}%)</title></rect>')
        parts.append(rect(x,base+7,21,8,color))
        count=int(rs[labels[0]]['primary_cells']); h=count/maxcount*95;parts.append(rect(x,cbase-h,21,h,color))
        ay=abase-(float(r['age'])-15)/35*95;parts.append(f'<circle cx="{cx}" cy="{ay}" r="4.5" fill="{color}"><title>{esc(r["donor_id"])}: {r["age"]} years; {count:,} primary {unit}</title></circle>')
        batch=r['capture_batch' if ci==0 else 'library_group_Lane']
        parts += [text(cx,800,batch,11,extra=f'transform="rotate(65 {cx} 800)"'),text(cx,877,r['donor_id'],10,extra=f'transform="rotate(65 {cx} 877)"')]
    parts += [rect(ox+100,965,13,13,'#325f87'),text(ox+120,977,'non-ENDO surgical controls',12),rect(ox+395,965,13,13,'#c45e3b'),text(ox+415,977,'ENDO cases',12)]
    for i,label in enumerate(labels):
        x=ox+100+(i%2)*325;y=1002+(i//2)*22;parts += [rect(x,y,12,12,cmap[label]),text(x+18,y+11,label,11)]
parts += [text(35,1161,'Denominators: all retained primary cells/nuclei in that donor. Observed zero labels remain zero; annotation absence is not biological absence.',13),text(35,1183,'Complete selected cohorts shown; 44/44 ages available. Batch strings are source-specific. No disease test, subtype fit or causal adjustment is added.',13)]
save(parts,'reports/figures/donor_composition_qc_v1.svg')
# Figure 5: retain every frozen candidate and distinguish unknown from tested negative.
parts=start(1500,1515,'Frozen 40-gene evidence matrix','All 40 frozen B0/G1 pool genes; external expression direction, full-family q, source genetics, conflict status and intervention uncertainty. Unavailable mapping is not a negative genetic test.')
parts += [text(30,38,'Frozen pool: 40 genes, separate evidence layers',25),text(30,66,'No supported target shortlist: 0 same-direction external q < 0.05; 1 limited genetic mapping / 40; intervention direction unknown / 40.',15),text(30,91,'Genes are ordered by B0 rank; selection badges identify the frozen expression (B0) or network-propagation (G1) top 20.',14)]
cols=[(35,'Gene'),(155,'List'),(225,'B0 rank'),(315,'G1 rank'),(415,'External direction'),(625,'External q'),(755,'Genetic mapping (one source)'),(1085,'Conflict'),(1275,'Intervention')]
for x,label in cols:parts.append(text(x,130,label,13))
summary=Counter()
for i,r in enumerate(sorted(pool,key=lambda r:int(r['B0_rank']))):
    y=151+i*29
    parts.append(rect(25,y-16,1450,28,'#f4f7fa' if i%2==0 else 'white'))
    ev=r['evaluable']=='TRUE'; ag=r['agree']=='TRUE'
    direction='Same sign' if ev and ag else 'Opposite sign' if ev else 'Not evaluable'
    color='#087f8c' if ev and ag else '#ad5d20' if ev else '#697386';summary[direction]+=1
    mapping='Limited positional / QTL' if r['genetic_genetic_status']=='supporting_mapping_limited' else 'Unavailable in audited source'
    badge='B0' if r['in_B0_top20']=='TRUE' else 'G1'
    vals=[r['gene_symbol'],badge,r['B0_rank'],r['G1_rank'],direction,f"{float(r['external_q']):.3f}" if r['external_q'] not in ('','NA') else 'NA',mapping,'Not established','Unknown']
    for (x,_),value in zip(cols,vals):parts.append(text(x,y,value,13,fill=color if x==415 else '#243247'))
    if not ev:parts.append(text(415,y+11,'expression filter failed',9))
parts += [text(30,1351,'Direction = sign agreement with the frozen discovery coefficient; opposite signs alone do not establish a precise biological conflict.',14),text(30,1375,'q = BH adjustment across all 17,643 externally tested genes. A large q means insufficient evidence here; it does not prove no effect.',14),text(30,1399,'Unavailable genetics = no mapping in the prespecified 2023 source, not 39 negative tests. SKIDA1 annotations reuse one publication.',14),text(30,1423,'Not established conflict / unknown intervention = evidence not established, not absence of conflict or proof that activation/inhibition helps.',14),text(30,1455,'The deterministic dossier rule selected only inconclusive PARP9 (hold); no precise opposite-effect or supported dossier emerged.',14),text(30,1480,'Rows include all frozen pool members; full IDs, effect intervals, missing reasons and provenance remain in replication_v2_dossier_pool.tsv.',13)]
save(parts,'reports/figures/evidence_matrix_v1.svg')
# Figure 1: audited mutually exclusive donor eligibility stages, not sample counts.
disc=read('manifests/discovery_clinical_candidates.tsv'); heca=read('manifests/final_heca_clinical_donors.tsv')
def stage(rows,name,disease): return {'stage':name,'total':len(rows),'cases':sum(r[disease]=='ENDO' for r in rows),'controls':sum(r[disease]=='non-ENDO' for r in rows)}
ds=[stage(disc,'Source clinical donors','endometriosis status')]
disc=[r for r in disc if r['author_proliferative_candidate']=='True']; ds.append(stage(disc,'Author proliferative phase','endometriosis status'))
disc=[r for r in disc if r['Minor exclusion criteria logical']=='0'];ds.append(stage(disc,'No source minor exclusion flag','endometriosis status'))
disc=[r for r in disc if not r['clinical_disease_inconsistency']];ds.append(stage(disc,'Disease-label consistency','endometriosis status'))
assert {r['sample_id'] for r in disc}=={r['donor_id'] for r in read('manifests/final_discovery_clinical_selected.tsv')}
ds.append(dict(ds[-1],stage='At least 50 primary cells'))
hs=[stage(heca,'Source donor ledger','disease')]
heca=[r for r in heca if r['phase_original']=='Proliferative'];hs.append(stage(heca,'Author proliferative phase','disease'))
heca=[r for r in heca if 'mixed_proliferative_menstrual_phase_not_primary' not in r['reasons']];hs.append(stage(heca,'Exclude mixed menstrual phase','disease'))
heca=[r for r in heca if r['primary_selected']=='True'];hs.append(stage(heca,'At least 50 primary nuclei','disease'))
assert {r['donor_id'] for r in heca}=={r['donor_id'] for r in read('manifests/final_heca_clinical_selected.tsv')}
assert [(r['total'],r['cases'],r['controls']) for r in ds]==[(60,35,25),(28,15,13),(23,12,11),(22,12,10),(22,12,10)]
assert [(r['total'],r['cases'],r['controls']) for r in hs]==[(63,44,19),(24,16,8),(23,16,7),(22,15,7)]
parts=start(1500,990,'Source donors to primary analyzed cohorts','Sequential donor selection in discovery and HECA; exclusions assigned to first failed stage, no double counting. Case/control counts shown at each stage.')
parts += [text(35,42,'Donor flow: source eligibility to fixed primary cohorts',25),text(35,72,'All counts are people, not libraries or cell-level observations. Sequential losses are mutually exclusive.',16)]
for ci,(stages,title,final) in enumerate([(ds,'Discovery: GSE266265','22 donors / 150,788 primary cells'),(hs,'Replication: HECA E-MTAB-14039','22 donors / 52,744 primary nuclei')]):
    x=35+ci*745;parts.append(text(x,126,title,22))
    for i,r in enumerate(stages):
        y=158+i*125;parts.append(rect(x,y,680,86,'#eef3f8', 'rx="8"'));parts += [text(x+18,y+28,r['stage'],18),text(x+18,y+61,f"{r['total']} donors = {r['cases']} ENDO + {r['controls']} non-ENDO",19)]
        if i<len(stages)-1:
            nr=stages[i+1];parts += [line(x+330,y+87,x+330,y+119,'#56718b',2),text(x+350,y+110,f"Excluded {r['total']-nr['total']} ({r['cases']-nr['cases']} / {r['controls']-nr['controls']})",13)]
    parts.append(text(x,817,final,20))
parts += [text(35,860,'Exclusion parentheses = cases / controls. HECA mixed phase is removed before the count threshold; FX1141 is counted once.',14),text(35,886,'HECA: 68 source samples represent 63 donors; 26 proliferative samples represent 24 donors. Two retained donors each use one selected library.',14),text(35,912,'The final selected-library rule excludes 2,741 nuclei without removing those two people; FX1215 has zero primary nuclei and is excluded.',14),text(35,938,'Author phase and hormone-eligibility evidence remain inherited. Source separation is documented; genotype-confirmed individual non-overlap is unavailable.',14),text(35,964,'Donor minimums and full-rank model gates apply after these selections. This diagram does not establish cell purity or causal disease effects.',14)]
save(parts,'reports/figures/donor_flow_v1.svg')
# Figure 4 sensitivity: display frozen matched sets; do not generate new sets.
nulls=read('work/statistical/heca_primary/evaluation/null_sets.tsv');ns=json.loads((ROOT/'work/statistical/heca_primary/evaluation/null_summary.json').read_text())
assert len(nulls)==600 and all(r['coverage_adequate']=='TRUE' for r in nulls)
parts=start(1500,650,'Expression and degree matched gene-set sensitivity','Each panel contains 200 frozen matched sets; bars show number of sets by external direction agreement. The observed line uses each fixed list evaluable denominator. Percentiles are descriptive, not p values.')
parts += [text(30,40,'Frozen expression / degree matched gene sets',25),text(30,70,'200 sets per method; 20 genes per set; shared genes and donors mean these are not independent biological replications.',15)]
for i,(method,rate,frac) in enumerate([('B0',14/20,'14 / 20 evaluable'),('B1',3/13,'3 / 13 evaluable'),('G1',6/19,'6 / 19 evaluable')]):
    rows=[r for r in nulls if r['method']==method];assert len(rows)==200
    bins=[0]*20
    for r in rows:bins[min(19,int(float(r['concordance'])*20))]+=1
    x=65+i*490;y0=445;w=400;h=230;limit=60
    assert max(bins)<=limit
    parts += [text(x,125,f'{method}: observed {frac}',19),text(x,152,f"Descriptive percentile: {ns[method]['observed_descriptive_percentile']} (not a p value)",13)]
    if method=='B1':parts.append(text(x,178,'Coverage limited: 13 / 20 < fixed 16 / 20 gate',13,fill='#ad5d20'))
    for val in [0,20,40,60]:
        y=y0-val/limit*h;parts += [line(x,y,x+w,y),text(x-10,y+4,str(val),12,extra='text-anchor="end"')]
    for j,n in enumerate(bins):parts.append(rect(x+j*w/20+1,y0-n/limit*h,w/20-2,n/limit*h,'#91b8ce'))
    xx=x+rate*w;parts += [line(xx,y0-h,xx,y0,'#ad5d20',3),text(x,202,'Number of matched sets',12)]
    for val in [0,.25,.5,.75,1]:parts.append(text(x+val*w,y0+25,f'{val:.0%}',12,extra='text-anchor="middle"'))
    parts.append(text(x+45,y0+53,'External same-sign / evaluable genes',13))
parts += [text(30,555,'Matching uses discovery expression and network degree only; external disease outcomes do not select the sets.',14),text(30,581,'All 600 frozen sets meet coverage requirements. Observed lists retain their own missingness; no replacement genes are added.',14),text(30,607,'Orange line = observed fixed list. The descriptive position does not establish ranking accuracy, therapeutic validity or a calibrated significance test.',14),text(30,633,'See the paired external-donor bootstrap and discovery leave-one-donor-out panels for their distinct uncertainty/stability questions.',14)]
save(parts,'reports/figures/matched_null_sensitivity_v1.svg')

outputs=['reports/figures/donor_composition_qc_v1.svg','reports/figures/evidence_matrix_v1.svg','reports/figures/donor_flow_v1.svg','reports/figures/matched_null_sensitivity_v1.svg']
receipt={'schema_version':1,'created_utc':datetime.now(timezone.utc).isoformat(),'purpose':'Render existing descriptive donor and frozen evidence tables; no scientific model rerun','input_sha256':{p:sha(ROOT/p) for p in INPUTS},'renderer_sha256':sha(Path(__file__)),'output_sha256':{p:sha(ROOT/p) for p in outputs},'checks':checks,'donor_flow':{'discovery':ds,'heca':hs},'matched_sets':{'per_method':200,'total':600,'descriptive_not_p_value':True},'pool_genes':40,'pool_external_direction_counts':dict(summary),'new_disease_tests':0,'scientific_fits_executed':0,'causal_or_clinical_validation':False,'scope_note':'Source-specific author labels, all selected donor denominators and unavailable/unknown evidence retained.'}
(ROOT/'receipts/summary_figures_v1.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'outputs':outputs,'checks':checks,'pool_genes':40,'direction_counts':dict(summary)},indent=2))
