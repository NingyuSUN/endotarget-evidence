#!/usr/bin/env python3
"""Read clinical columns only; no publication outcome columns or expression."""
import csv, collections, datetime, hashlib, json, pathlib, re, zipfile
import xml.etree.ElementTree as ET
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE=ROOT/"data/input_routes/osf/selected/04_Analysis/S3_COMPLETION_EXTENSION"
XLSX=BASE/"GSE266265_PUBLISHED_SCREEN/sources/41467_2026_73020_MOESM3_ESM.xlsx"
SECONDARY=BASE/"GSE266265_DONOR_RECONSTRUCTION/PRIMARY_CLINICAL_MANIFEST.csv"
EXPECTED_SHA="5749f05690ffc1c7b6ff7492eba047b5870610c3f6387a47d7342e67538bce69"
NS={"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
COLS={"A":"Patient sample ID","B":"endometriosis status","C":"endometriosis grade",
"D":"Deep Infiltrating Endometriosis","E":"Endometrioma","F":"Peritoneal","G":"Adenomyosis",
"H":"Patient.age","I":"Hormonal treatment/IUD","O":"Menstrual cycle day (Base 1)",
"P":"Progesterone in serum (nmol/L)","R":"pathologist HE evaluation","T":"MenstrualCyclePhase",
"U":"10x capturing batch","V":"sequencing batch","W":"Minor exclusion criteria logical","X":"Minor exclusion criteria detail"}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_tsv(p,rs):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter="\t");w.writeheader();w.writerows(rs)
def read_workbook():
 if sha(XLSX)!=EXPECTED_SHA:raise ValueError("Archived workbook differs from audited OSF member")
 with zipfile.ZipFile(XLSX) as z:
  ss=["".join(x.itertext()) for x in ET.fromstring(z.read("xl/sharedStrings.xml"))]
  rows=[]
  for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).findall(".//s:sheetData/s:row",NS):
   d={}
   for c in r:
    col=re.sub("[0-9]","",c.get("r"))
    if col not in COLS:continue
    if c.find("s:f",NS) is not None:raise ValueError("Unexpected clinical formula: "+c.get("r"))
    v=c.find("s:v",NS);v=v.text if v is not None else ""
    if c.get("t")=="s" and v:v=ss[int(v)]
    if c.get("t")=="inlineStr":v="".join(c.itertext())
    d[COLS[col]]=v
   if r.get("r")=="1":
    assert all(d.get(v)==v for v in COLS.values()),"Workbook schema mismatch"
   elif d.get("Patient sample ID"):
    d["xlsx_row"]=int(r.get("r"));rows.append(d)
 return rows
def counts(rs):return dict(sorted(collections.Counter(r["endometriosis status"] for r in rs).items()))
def main():
 raw=read_workbook()
 clinical={r["Patient sample ID"]:r for r in raw if re.fullmatch("sample[0-9]+",r["Patient sample ID"])}
 assert len(clinical)==60
 assert sum(bool(re.fullmatch("sample[0-9]+",r["Patient sample ID"])) for r in raw)==len(clinical),"Duplicate atlas patient identifier"
 geo=list(csv.DictReader((ROOT/"reports/clinical_gate_discovery_gsm.tsv").open(),delimiter="\t"))
 secondary=list(csv.DictReader(SECONDARY.open()))
 sec={r["Patient sample ID"]:r for r in secondary}
 assert len(sec)==len(secondary)==23
 output=[]
 for g in geo:
  match=re.fullmatch("(CTL|ENDO)_(sample[0-9]+)",g["sample_title"])
  assert match,g["sample_title"]
  sid=match.group(2);c=clinical[sid]
  assert (c["endometriosis status"]=="ENDO")== (g["disease_original"]=="Endometriosis")
  assert c["endometriosis grade"]==g["grade_original"]
  row={"sample_id":sid,"GSM":g["gsm"],"GEO_title":g["sample_title"],
       "GEO_disease_original":g["disease_original"],**c}
  row["hormone_interpretation"]="explicit_FALSE_plus_paper_prior_3_month_exclusion"
  row["age_years"]=int(c["Patient.age"])
  row["clinical_disease_inconsistency"]="non_ENDO_with_endometriosis_lesion_flag" if c["endometriosis status"]=="non-ENDO" and any(c[k]=="1" for k in ["Deep Infiltrating Endometriosis","Endometrioma","Peritoneal"]) else ""
  hist=c.get("pathologist HE evaluation","")
  row["phase_histology_conflict"]=c["MenstrualCyclePhase"]=="proliferative" and "secretory" in hist.lower()
  row["clinical_phase_rule_disagreement"]=c["MenstrualCyclePhase"]=="proliferative" and not(float(c["Progesterone in serum (nmol/L)"])<2 and int(c["Menstrual cycle day (Base 1)"])<=15)
  row["author_proliferative_candidate"]=c["MenstrualCyclePhase"]=="proliferative"
  row["proposal_exclude_minor_clinical_flag"]=c["Minor exclusion criteria logical"]=="1"
  row["proposal_conservative_clinical_candidate"]=bool(row["author_proliferative_candidate"] and c["Hormonal treatment/IUD"]=="FALSE" and not row["proposal_exclude_minor_clinical_flag"] and not row["clinical_disease_inconsistency"] and not row["phase_histology_conflict"])
  row["secondary23_member"]=sid in sec
  row["secondary_clinical_match"]="not_in_secondary23"
  if sid in sec:
   s=sec[sid]
   assert s["GSM"]==g["gsm"] and s["disease"]==c["endometriosis status"] and s["phase"]==c["MenstrualCyclePhase"] and int(s["age"])==row["age_years"]
   assert s["batch"]==c["10x capturing batch"] and s["Hormonal treatment/IUD"].lower()=="false"
   row["secondary_clinical_match"]="GSM_disease_phase_age_batch_hormone_match"
  row["source_sheet"]="Sheet1"
  row["source_sha256"]=EXPECTED_SHA
  row["provenance"]="OSF_archived_publication_supplement_copy;publisher_byte_equivalence_NOT_VERIFIED"
  row["final_eligibility"]="NOT_VERIFIED"
  output.append(row)
 assert len(output)==60 and len({r["sample_id"] for r in output})==60
 phase=[r for r in output if r["author_proliferative_candidate"]]
 minor=[r for r in phase if not r["proposal_exclude_minor_clinical_flag"]]
 disease=[r for r in minor if not r["clinical_disease_inconsistency"]]
 conservative=[r for r in disease if not r["phase_histology_conflict"]]
 noadeno=[r for r in conservative if r["Adenomyosis"]!="1"]
 assert {r["sample_id"] for r in minor}==set(sec),"Secondary23 selection does not match observed simple clinical rule"
 def describe(rs):
  ct=counts(rs)
  return {"n":len(rs),"groups":ct,"minimum_8_per_group_numerically_met":ct.get("ENDO",0)>=8 and ct.get("non-ENDO",0)>=8,"IDs":[r["sample_id"] for r in rs],"age_range_by_group":{g:[min(r["age_years"] for r in rs if r["endometriosis status"]==g),max(r["age_years"] for r in rs if r["endometriosis status"]==g)] for g in ct}}
 summary={
 "status":"METADATA_SCENARIOS_ONLY_M0_NOT_PASSED",
 "source_rows_nonempty_id":len(raw),"atlas_geo_matched":len(output),"IHC_rows_excluded":len(raw)-len(output),
 "GEO_clinical_disease_and_grade_matches":60,"exact_age_available":60,"explicit_hormone_FALSE":sum(r["Hormonal treatment/IUD"]=="FALSE" for r in output),
 "scenarios":{"author_proliferative":describe(phase),"exclude_minor_clinical_flags_PROPOSED":describe(minor),"also_quarantine_disease_inconsistency_PROPOSED":describe(disease),"also_quarantine_phase_histology_conflict_PROPOSED":describe(conservative),"additionally_exclude_adenomyosis_SENSITIVITY_NOT_PRIMARY":describe(noadeno)},
 "minor_flagged_proliferative":[{"id":r["sample_id"],"condition":r["Minor exclusion criteria detail"]} for r in phase if r["proposal_exclude_minor_clinical_flag"]],
 "disease_inconsistencies":[r["sample_id"] for r in output if r["clinical_disease_inconsistency"]],
 "phase_histology_conflicts":[r["sample_id"] for r in phase if r["phase_histology_conflict"]],
 "clinical_phase_rule_disagreement":[r["sample_id"] for r in phase if r["clinical_phase_rule_disagreement"]],
 "secondary23_ids_match_minor_filtered_scenario":True,
 "batch_cross_tab_conservative":{b:counts([r for r in conservative if r["10x capturing batch"]==b]) for b in sorted({r["10x capturing batch"] for r in conservative})},
 "sequencing_cross_tab_conservative":{b:counts([r for r in conservative if r["sequencing batch"]==b]) for b in sorted({r["sequencing batch"] for r in conservative})},
 "cell_lineage_counts_checked":False,"design_matrix_checked":False,"no_expression_effects":True,
 }
 out=ROOT/"manifests/discovery_clinical_candidates.tsv";write_tsv(out,output)
 summ=ROOT/"reports/discovery_clinical_gate_summary.json";summ.write_text(json.dumps(summary,indent=2)+"\n")
 receipt={"utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"command":"python3 scripts/discovery_clinical_gate.py","source_workbook":str(XLSX.relative_to(ROOT)),"source_sha256":sha(XLSX),"source_columns_read":COLS,"sheet":"Sheet1","excluded_columns":"symptoms_and_pseudotime_and_published_analysis_inclusion_flags_not_used","publisher_byte_equivalence":"NOT_VERIFIED","archived_source_status":"OSF_manifest_SHA256_and_member_CRC_verified_by_source_retrieval","inputs":[{"path":str(p.relative_to(ROOT)),"sha256":sha(p)} for p in [XLSX,SECONDARY,ROOT/"reports/clinical_gate_discovery_gsm.tsv"]],"outputs":[{"path":str(p.relative_to(ROOT)),"sha256":sha(p)} for p in [out,summ]],"source_workbook_unchanged":sha(XLSX)==EXPECTED_SHA,"script_sha256":sha(pathlib.Path(__file__))}
 (ROOT/"receipts/discovery_clinical_gate_join.json").write_text(json.dumps(receipt,indent=2)+"\n")
 print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
