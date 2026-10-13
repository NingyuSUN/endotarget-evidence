#!/usr/bin/env python3
"""Final clinical/sample and metadata-design gate. No expression arrays/effects."""
import collections,csv,datetime,gzip,hashlib,json,pathlib,re
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
LABELS=("eStromal_MMPs","eStromal","eStromal_cycling","dStromal_early","dStromal_mid","dStromal_late","Fibroblast_basalis")
UNKNOWN={"","NA","nan"}
def read(path,delimiter="\t"):
 with (ROOT/path).open() as f:return list(csv.DictReader(f,delimiter=delimiter))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dump(path,obj):(ROOT/path).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n")
def table(path,rows):
 if not rows:raise ValueError("empty ledger")
 with (ROOT/path).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def unique(rows,key):
 out={}
 for r in rows:
  k=r[key]
  if k in out and out[k]!=r:raise ValueError("conflicting duplicate "+key+":"+k)
  out[k]=r
 return out
def model(rows,cat_fields,minimum,adeno=True):
 names=["intercept","ENDO_minus_non_ENDO","age_centered_years"]
 cols=[np.ones(len(rows)),np.array([r["disease"]=="ENDO" for r in rows],float)]
 age=np.array([r["age"] for r in rows],float);cols.append(age-age.mean())
 levels={}
 for field in cat_fields:
  values=sorted({str(r[field]) for r in rows});levels[field]=values
  if "" in values:raise ValueError("missing categorical covariate "+field)
  for v in values[1:]:
   names.append(field+":"+v+"_vs_"+values[0]);cols.append(np.array([str(r[field])==v for r in rows],float))
 adeno_omitted=False
 if adeno and len({int(r["adenomyosis"]) for r in rows})>1:
  names.append("adenomyosis");cols.append(np.array([r["adenomyosis"] for r in rows],float))
 elif adeno:adeno_omitted=True
 x=np.column_stack(cols)
 assert np.isfinite(x).all()
 rank=int(np.linalg.matrix_rank(x));nuis=np.delete(x,1,axis=1)
 identifiable=rank>np.linalg.matrix_rank(nuis)
 disease_resid=x[:,1]-nuis@np.linalg.lstsq(nuis,x[:,1],rcond=None)[0]
 groups=dict(collections.Counter(r["disease"] for r in rows))
 countpass=groups.get("ENDO",0)>=minimum and groups.get("non-ENDO",0)>=minimum
 result={"n":len(rows),"groups":groups,"columns":names,"rank":rank,"n_columns":len(names),"residual_df":len(rows)-rank,"full_rank":rank==len(names),"disease_identifiable":bool(identifiable),"minimum_per_group":minimum,"donor_minimum_pass":countpass,"structural_gate_pass":bool(countpass and rank==len(names) and len(rows)-rank>=5 and identifiable),"disease_information_unweighted":float(disease_resid@disease_resid),"categorical_levels":levels,"adenomyosis_constant_column_omitted":adeno_omitted,"donor_order":[r["donor_id"] for r in rows]}
 return x,result
def main():
 samples=read("reports/heca_h5_sample_crosswalk.tsv")
 donor_original=read("data/metadata/heca_author/metadata-donors.csv",",")
 donor_groups=collections.defaultdict(list)
 for r in donor_original:donor_groups[r["Donor"]].append(r)
 donor={k:v[0] for k,v in donor_groups.items()}
 donor_conflicts={k for k,v in donor_groups.items() if any(r!=v[0] for r in v)}
 library=unique(read("data/metadata/heca_author/metadata-libraries.csv",","),"FolderName")
 deconv=unique([r for r in read("data/metadata/heca_author/metadata-libraries_deconvoluted.csv",",") if r["Library_genotype"] in {x["h5_sample"] for x in samples}],"Library_genotype")
 by_sample=unique(samples,"h5_sample")
 obs_path=ROOT/"data/metadata/heca_h5/obs_selected.tsv"
 # Original row order is original H5 obs order; row0 is first data row.
 obs_counts=collections.Counter();total_counts=collections.Counter();index_by_sample=collections.defaultdict(list)
 with obs_path.open() as f:
  for index,r in enumerate(csv.DictReader(f,delimiter="\t")):
   sample=r["sample"];assert sample in by_sample
   total_counts[sample]+=1
   if r["celltype"] in LABELS and r["discard_lowQC"]=="no":
    obs_counts[sample]+=1;index_by_sample[sample].append(index)
 assert index+1==312246
 prepared=[]
 for r in samples:
  d=donor[r["genotype"]];l=library[r["library"]];dc=deconv[r["h5_sample"]]
  assert dc["Donor"]==r["genotype"] and d["Binary Stage"]==r["Binary Stage"]
  assert float(d["Age"])==float(r["Age"])
  assert obs_counts[r["h5_sample"]]==int(r["primary_stroma"])
  assert total_counts[r["h5_sample"]]==int(r["nuclei"])
  pathology=d["Endometrial pathology_A-F-H-P"]
  tokens=set(x.strip() for x in pathology.split("+"))
  assert tokens.issubset({"A","C","E","F","H","P"}),tokens
  hormone_conflict=r["Hormonal treatment"] not in UNKNOWN or d["Hormonal treatment"] not in UNKNOWN
  eligibility_base=r["genotype"] not in donor_conflicts and r["Binary Stage"]=="Proliferative" and not hormone_conflict and r["Biopsy_type"]=="Superficial_Eutopic" and r["Tissue_sampled"]=="eutopic_endometrium" and r["Stage"]!="Proliferative Menstrual"
  prepared.append({
   "donor_id":r["genotype"],"sample_id":r["h5_sample"],"library_id":r["library"],"disease":"ENDO" if r["Group"]=="Endo_Superficial" else "non-ENDO",
   "Group_original":r["Group"],"age":float(r["Age"]),"phase_original":r["Binary Stage"],"phase_fine_original":r["Stage"],
   "histology_phase_original":d["MensPhase3_pathologist"],"expression_phase_original":d["Binary Stage from scRNAseq data"],"day_original":d["Day (self-reported)"],
   "hormone_h5_original":r["Hormonal treatment"],"hormone_clinical_original":d["Hormonal treatment"],
   "hormone_eligibility_basis":"study_defined_natural_cycle_group" if eligibility_base else "not_eligible_natural_cycle_primary",
   "hormone_explicit_negative_observed":False,"hormone_explicit_contradiction":hormone_conflict,
   "clinical_duplicate_conflict":r["genotype"] in donor_conflicts,"clinical_duplicate_rows_json":json.dumps(donor_groups[r["genotype"]],ensure_ascii=False) if r["genotype"] in donor_conflicts else "",
   "pathology_original":"|".join(sorted({x["Endometrial pathology_A-F-H-P"] for x in donor_groups[r["genotype"]]})),"adenomyosis":(-1 if r["genotype"] in donor_conflicts else int("A" in tokens)),"adenomyosis_coding_basis":"author_pathology_A_token_with_paper_benign_pathology_terms",
   "Other_drug_treatment_original":d["Other drug treatment"],"clinical_flag_original":d["Flag"],"Endometrium_histology_original":d["Endometrium histology"],
   "biopsy_type":r["Biopsy_type"],"tissue":r["Tissue_sampled"],"study_original":l["StudyName"],"library_group_Lane":l["Lane"],
   "kit_original":l["10x kit"],"processing_original":l["Processing"],"preparation_Treatment_original":l["Treatment"],
   "library_QC_original":l["QC_single sample analysis"],"library_CellRanger_QC_original":l["QC_CellRanger"],"author_library_analysed":l["analysed"],
   "sdrf_declared_individual_original":r["sdrf_declared_individual"],"author_deconvolution_Donor":dc["Donor"],
   "primary_nuclei":obs_counts[r["h5_sample"]],"total_nuclei":total_counts[r["h5_sample"]],
   "source_clinical_candidate":eligibility_base,"histology_conflict":"secretory" in d["MensPhase3_pathologist"].lower(),
   "repeat_relation":"NOT_VERIFIED","selected_primary":False,"exclusion_reason":"",
  })
 by_donor=collections.defaultdict(list)
 for r in prepared:by_donor[r["donor_id"]].append(r)
 selected=[];donor_ledger=[];repeat_loss=[]
 for did,rr in sorted(by_donor.items()):
  eligible=[r for r in rr if r["source_clinical_candidate"] and r["primary_nuclei"]>=50]
  best=sorted(eligible,key=lambda r:(-r["primary_nuclei"],r["library_id"],r["sample_id"]))[0] if eligible else None
  for r in rr:
   if r["phase_fine_original"]=="Proliferative Menstrual":r["exclusion_reason"]="mixed_proliferative_menstrual_phase_not_primary"
   elif not r["source_clinical_candidate"]:r["exclusion_reason"]="not_author_proliferative_or_other_clinical_gate"
   elif r["primary_nuclei"]<50:r["exclusion_reason"]="primary_nuclei_below_50"
   elif r is not best:r["exclusion_reason"]="other_library_selected_same_donor_relation_unknown"
   else:r["selected_primary"]=True;selected.append(r)
  representative=best or sorted(rr,key=lambda r:(-r["primary_nuclei"],r["library_id"]))[0]
  donor_ledger.append({"donor_id":did,"disease":representative["disease"],"age":representative["age"],"phase_original":representative["phase_original"],"samples_observed":"|".join(r["sample_id"] for r in rr),"sample_count":len(rr),"selected_sample_id":best["sample_id"] if best else "","primary_selected":best is not None,"selected_primary_nuclei":best["primary_nuclei"] if best else 0,"excluded_primary_nuclei":sum(r["primary_nuclei"] for r in rr if r is not best),"reasons":"|".join(sorted({r["exclusion_reason"] for r in rr if r["exclusion_reason"]}))})
  if best and len(rr)>1:repeat_loss.append({"donor_id":did,"selected_sample":best["sample_id"],"selected_primary_nuclei":best["primary_nuclei"],"excluded_samples":[r["sample_id"] for r in rr if r is not best],"excluded_primary_nuclei":sum(r["primary_nuclei"] for r in rr if r is not best),"excluded_total_nuclei":sum(r["total_nuclei"] for r in rr if r is not best),"relation":"NOT_VERIFIED"})
 selected=sorted(selected,key=lambda r:r["donor_id"])
 assert len(selected)==22 and dict(collections.Counter(r["disease"] for r in selected))=={"ENDO":15,"non-ENDO":7}
 table("manifests/final_heca_clinical_samples.tsv",prepared)
 table("manifests/final_heca_clinical_donors.tsv",donor_ledger)
 table("manifests/final_heca_clinical_selected.tsv",selected)
 indices=np.array(sorted(i for r in selected for i in index_by_sample[r["sample_id"]]),dtype=np.int64)
 assert len(np.unique(indices))==len(indices)==sum(r["primary_nuclei"] for r in selected)
 outdir=ROOT/"work/final_clinical";outdir.mkdir(exist_ok=True)
 np.save(outdir/"heca_primary_h5_obs_indices.npy",indices,allow_pickle=False)
 # Discovery already source-adjudicated phase; retain both conflict sensitivities.
 disc=[];coverage={r["sample"]:r for r in read("data/input_routes/osf/independent_donor_coverage.csv",",") if r["population_id"]=="P077"}
 for r in read("manifests/discovery_clinical_candidates.tsv"):
  if r["author_proliferative_candidate"]!="True" or r["Minor exclusion criteria logical"]!="0" or r["clinical_disease_inconsistency"]:continue
  c=coverage[r["sample_id"]]
  assert c["GSM"]==r["GSM"] and c["disease"]==r["endometriosis status"] and int(c["cells_independent"])>=50
  disc.append({"donor_id":r["sample_id"],"sample_id":r["sample_id"],"GSM":r["GSM"],"disease":r["endometriosis status"],"age":float(r["age_years"]),"capture_batch":r["10x capturing batch"],"sequencing_batch":r["sequencing batch"],"adenomyosis":int(r["Adenomyosis"]),"primary":True,"strict_histology_sensitivity":r["phase_histology_conflict"]!="True","adenomyosis_free_sensitivity":r["Adenomyosis"]=="0","author_phase":r["MenstrualCyclePhase"],"histology_phase_original":r["pathologist HE evaluation"],"hormone_original":r["Hormonal treatment/IUD"],"hormone_eligibility_basis":"explicit_FALSE_and_paper_prior_3_month_exclusion","clinical_disease_inconsistency":r["clinical_disease_inconsistency"],"minor_clinical_flag_original":r["Minor exclusion criteria logical"],"minor_clinical_detail_original":r["Minor exclusion criteria detail"],"primary_lineage_cells":int(c["cells_independent"]),"source_xlsx_row":r["xlsx_row"],"source_workbook_sha256":r["source_sha256"]})
 disc=sorted(disc,key=lambda r:int(r["donor_id"].replace("sample","")))
 assert len(disc)==22
 table("manifests/final_discovery_clinical_selected.tsv",disc)
 scenarios={"heca_primary":(selected,["library_group_Lane"],5),"heca_strict_histology_conflict_exclusion":([r for r in selected if not r["histology_conflict"]],["library_group_Lane"],5),"heca_adenomyosis_free":([r for r in selected if r["adenomyosis"]==0],["library_group_Lane"],5),"discovery_primary":(disc,["capture_batch"],8),"discovery_strict_histology_conflict_exclusion":([r for r in disc if r["strict_histology_sensitivity"]],["capture_batch"],8),"discovery_adenomyosis_free":([r for r in disc if r["adenomyosis"]==0],["capture_batch"],8),"discovery_strict_and_adenomyosis_free":([r for r in disc if r["strict_histology_sensitivity"] and r["adenomyosis"]==0],["capture_batch"],8)}
 designs={}
 for name,(rows,fields,minimum) in scenarios.items():
  x,result=model(rows,fields,minimum);designs[name]=result
  table("reports/final_clinical_design_"+name+".tsv",[{"donor_id":r["donor_id"],**{k:float(v) for k,v in zip(result["columns"],x[i])}} for i,r in enumerate(rows)])
 for threshold in [30,100]:
  subset=[]
  for did,rr in sorted(by_donor.items()):
   er=[r for r in rr if r["source_clinical_candidate"] and r["primary_nuclei"]>=threshold]
   if er:subset.append(sorted(er,key=lambda r:(-r["primary_nuclei"],r["library_id"],r["sample_id"]))[0])
  _,designs["heca_primary_threshold_"+str(threshold)]=model(subset,["library_group_Lane"],5)
 # Granular library dummies leave too few residual degrees of freedom; reject.
 _,bad=model(selected,["library_id"],5)
 assert not bad["structural_gate_pass"]
 # Already selected rows have one biological unit each and immutable original source labels.
 assert len({r["donor_id"] for r in selected})==len(selected)
 nested={}
 for field in ["kit_original","preparation_Treatment_original"]:
  mapping={b:sorted({r[field] for r in selected if r["library_group_Lane"]==b}) for b in sorted({r["library_group_Lane"] for r in selected})}
  nested[field]={"levels_per_Lane":mapping,"nested_in_Lane":all(len(v)==1 for v in mapping.values())}
 summary={"status":"CLINICAL_SAMPLE_AND_DESIGN_GATE_EXECUTED_NOT_FULL_M0","selection_rule":"among clinically eligible libraries with >=50 primary nuclei choose largest primary nucleus count then lexical library_id then sample_id; no summing across unknown biopsy relation","primary_labels":list(LABELS),"original_h5_obs_rows":index+1,"selected_heca_samples":[r["sample_id"] for r in selected],"selected_heca_donors":[r["donor_id"] for r in selected],"selected_heca_primary_nuclei":len(indices),"selected_primary_min_nuclei":min(r["primary_nuclei"] for r in selected),"repeat_losses":repeat_loss,"excluded_low_count_candidates":[{"donor_id":r["donor_id"],"count":r["primary_nuclei"],"fine_phase":r["phase_fine_original"]} for r in prepared if r["source_clinical_candidate"] and r["primary_nuclei"]<50],"mixed_phase_exclusions":[{"donor_id":r["donor_id"],"phase_fine":r["phase_fine_original"],"primary_nuclei":r["primary_nuclei"]} for r in prepared if r["phase_fine_original"]=="Proliferative Menstrual"],"clinical_source_duplicate_conflicts":sorted(donor_conflicts),"hormone_coding":{"proliferative_samples":sum(r["phase_original"]=="Proliferative" for r in prepared),"proliferative_explicit_drug_conflicts":sum(r["phase_original"]=="Proliferative" and r["hormone_explicit_contradiction"] for r in prepared),"missing_not_relabelled_negative":True,"basis":"adopted_study_defined_natural_cycle_group"},"heca_technical_field_interpretation":{"Lane":"literal author library metadata grouping used categorically; physical sequencing lane/date semantics not asserted","library_id":"individual capture/library identifiers, not automatically a multi-donor technical batch","processing":"Frozen constant","kit_and_preparation_nested":nested},"heca_batch_disease_crosstab":{b:dict(collections.Counter(r["disease"] for r in selected if r["library_group_Lane"]==b)) for b in sorted({r["library_group_Lane"] for r in selected})},"heca_pathology_crosstab":{p:dict(collections.Counter(r["disease"] for r in selected if r["pathology_original"]==p)) for p in sorted({r["pathology_original"] for r in selected})},"designs":designs,"negative_guard_library_dummy_design":bad,"no_expression_read":True,"scientific_validation":False,"pending_other_M0_gates":["external_raw_counts","common_gene_IDs","lineage_crosswalk_final_integration","cohort_independence_integration","network_and_ranking_freeze"]}
 dump("reports/final_clinical_gate_summary.json",summary)
 inputs=["data/metadata/heca_h5/obs_selected.tsv","reports/heca_h5_sample_crosswalk.tsv","data/metadata/heca_author/metadata-donors.csv","data/metadata/heca_author/metadata-libraries.csv","data/metadata/heca_author/metadata-libraries_deconvoluted.csv","manifests/discovery_clinical_candidates.tsv","data/input_routes/osf/independent_donor_coverage.csv","docs/PRE_EFFECT_METHOD_DECISIONS.md"]
 outputs=["manifests/final_heca_clinical_samples.tsv","manifests/final_heca_clinical_donors.tsv","manifests/final_heca_clinical_selected.tsv","manifests/final_discovery_clinical_selected.tsv","reports/final_clinical_gate_summary.json","work/final_clinical/heca_primary_h5_obs_indices.npy"]
 receipt={"utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"command":"envs/m0-h5/bin/python scripts/final_clinical_gate.py","numpy_version":np.__version__,"script_sha256":sha(pathlib.Path(__file__)),"input_hashes":{p:sha(ROOT/p) for p in inputs},"output_hashes":{p:sha(ROOT/p) for p in outputs},"h5_indices_zero_based":True,"index_order":"ascending original H5 obs order, no header","index_count":len(indices),"clinical_and_design_checks_pass":designs["heca_primary"]["structural_gate_pass"] and designs["discovery_primary"]["structural_gate_pass"],"M0_pass":False}
 dump("receipts/final_clinical_gate.json",receipt)
 print(json.dumps({"designs":designs,"repeat_losses":repeat_loss,"selected_nuclei":len(indices)},indent=2))
if __name__=="__main__":main()
