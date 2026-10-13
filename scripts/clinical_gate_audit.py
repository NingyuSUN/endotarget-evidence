#!/usr/bin/env python3
"""Source-only clinical candidate audit; never certifies analysis eligibility."""
import collections
import csv
import datetime
import gzip
import hashlib
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def write_tsv(path, rows):
    if not rows:
        raise ValueError("No output rows")
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
def nodes(obj):
    if isinstance(obj, dict):
        yield obj
        for v in obj.values():
            yield from nodes(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from nodes(v)
def main():
    metadata = ROOT / "data/metadata"
    path = metadata / "E-MTAB-14039.sdrf.txt"
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    required = ["Source Name", "Characteristics[individual]", "Characteristics[disease]",
        "Characteristics[age]", "Characteristics[10x technology]", "Characteristics[library_id]",
        "Characteristics[menstrual_cycle_stage]", "Characteristics[menstrual_cycle_stage_fine]"]
    if not rows or not set(required).issubset(rows[0]):
        raise ValueError("SDRF required fields missing")
    if {r["Characteristics[age]"] for r in rows} != {"20 to 45"}:
        raise ValueError("Age source changed: re-audit exact ages before regenerating snapshot")
    if any("hormon" in k.lower() for k in rows[0]):
        raise ValueError("New hormone field requires clinical review")
    nuclei = [(i,r) for i,r in enumerate(rows,2) if r["Characteristics[10x technology]"] == "snRNAseq"]
    donors = collections.defaultdict(list)
    for i,r in nuclei:
        donors[r["Characteristics[individual]"]].append((i,r))
    bio = json.loads((metadata / "E-MTAB-14039.json").read_text())
    bio_records = []
    for obj in nodes(bio):
        if obj.get("type") == "Characteristics Table":
            rec = {a["name"]:a["value"] for a in obj.get("attributes", [])}
            bio_records.append(rec)
    def vals(rr,key):
        return "|".join(sorted({r[key] for _,r in rr}))
    candidates = []
    row_audit = []
    for donor, rr in sorted(donors.items()):
        for line,r in rr:
            same = [b for b in bio_records if b.get("individual")==donor and b.get("library_id")==r["Characteristics[library_id]"] and b.get("10x technology")=="snRNAseq"]
            ok = len(same)==1 and all(same[0].get(k)==r["Characteristics["+k+"]"] for k in ["disease","menstrual_cycle_stage","menstrual_cycle_stage_fine"])
            row_audit.append({
                "sdrf_line":line, "source_name":r["Source Name"], "source_suffix_literal":r["Source Name"].rsplit("_",1)[-1],
                "declared_individual":donor, "library_id":r["Characteristics[library_id]"],
                "disease":r["Characteristics[disease]"], "phase":r["Characteristics[menstrual_cycle_stage]"],
                "phase_fine":r["Characteristics[menstrual_cycle_stage_fine]"],
                "source_suffix_matches_individual":r["Source Name"].rsplit("_",1)[-1]==donor,
                "biostudies_declared_record_match":ok,
                "source":"data/metadata/E-MTAB-14039.sdrf.txt",
            })
        if not any(r["Characteristics[menstrual_cycle_stage]"]=="Proliferative" for _,r in rr):
            continue
        conflicts = [k for k in ["Characteristics[disease]","Characteristics[menstrual_cycle_stage]","Characteristics[age]"] if len({r[k] for _,r in rr})>1]
        candidates.append({
            "cohort":"E-MTAB-14039", "assay":"snRNAseq", "declared_individual":donor,
            "disease_original":vals(rr,"Characteristics[disease]"),
            "phase_original":vals(rr,"Characteristics[menstrual_cycle_stage]"),
            "phase_fine_original":vals(rr,"Characteristics[menstrual_cycle_stage_fine]"),
            "age_original":vals(rr,"Characteristics[age]"), "age_exact_years":"UNKNOWN",
            "age_interpretation":"range_only_not_individual_age",
            "tissue_original":vals(rr,"Characteristics[organism part]"),
            "hormone_status_donor_level":"UNKNOWN",
            "hormone_study_definition":"natural_cycle_requires_no_hormonal_treatment_at_least_3_months; membership_not_mapped",
            "phase_ascertainment_donor_level":"UNKNOWN",
            "phase_study_methods":"histology_two_pathologists_or_transcriptomic_fallback",
            "control_indication_donor_level":"UNKNOWN",
            "recruiting_study_donor_level":"UNKNOWN",
            "source_names":"|".join(r["Source Name"] for _,r in rr),
            "source_suffixes":"|".join(r["Source Name"].rsplit("_",1)[-1] for _,r in rr),
            "library_ids":"|".join(r["Characteristics[library_id]"] for _,r in rr),
            "sdrf_lines":"|".join(str(i) for i,_ in rr),
            "record_count":len(rr),
            "source_suffix_mismatch_records":sum(r["Source Name"].rsplit("_",1)[-1]!=donor for _,r in rr),
            "repeat_relation":"UNKNOWN_technical_or_repeat_biopsy" if len(rr)>1 else "single_observed_record",
            "within_declared_id_conflicts":"|".join(conflicts) or "none_in_observed_fields",
            "barcode_demultiplex_link":"NOT_VERIFIED",
            "candidate_status":"PROLIFERATIVE_METADATA_CANDIDATE_ONLY",
            "final_eligibility":"NOT_VERIFIED",
            "provenance":"data/metadata/E-MTAB-14039.sdrf.txt;data/metadata/E-MTAB-14039.json;data/reference/heca_article.xml#Sec10,Sec11,Sec12,Sec21",
        })
    soft = gzip.open(metadata / "GSE266265_family.soft.gz", "rt").read()
    discovery=[]
    for block in soft.split("^SAMPLE = ")[1:]:
        gsm=block.splitlines()[0].strip()
        attrs=collections.defaultdict(list)
        for line in block.splitlines()[1:]:
            if line.startswith("^"):
                break
            if line.startswith("!") and " = " in line:
                k,v=line.split(" = ",1); attrs[k].append(v)
        characteristics={}
        for value in attrs["!Sample_characteristics_ch1"]:
            k,sep,v=value.partition(":")
            if sep:characteristics[k.strip()]=v.strip()
        if set(characteristics) != {"tissue", "endometriosis status", "endometriosis grade", "Sex"}:
            raise ValueError("GEO characteristics changed: re-audit phase/age availability")
        discovery.append({
            "gsm":gsm, "sample_title":"|".join(attrs["!Sample_title"]),
            "tissue":characteristics.get("tissue","UNKNOWN"),
            "disease_original":characteristics.get("endometriosis status","UNKNOWN"),
            "grade_original":characteristics.get("endometriosis grade","UNKNOWN"),
            "sex":characteristics.get("Sex","UNKNOWN"),
            "phase_donor_level":"UNKNOWN", "age_donor_level":"UNKNOWN",
            "hormone_donor_record":"UNKNOWN",
            "hormone_study_evidence":"paper reports all 60 atlas donors excluded exogenous hormones within prior 3 months",
            "participant_independence":"paper reports 60 women; title-to-clinical-table linkage pending",
            "clinical_table_join":"NOT_VERIFIED",
            "relations":"|".join(attrs["!Sample_relation"]),
            "characteristics_original_json":json.dumps(characteristics,ensure_ascii=False,sort_keys=True),
            "source":"data/metadata/GSE266265_family.soft.gz#"+gsm,
            "final_eligibility":"NOT_VERIFIED",
        })
    if len({r["gsm"] for r in discovery}) != len(discovery):
        raise ValueError("duplicate GSM accession")
    write_tsv(ROOT/"manifests/clinical_candidates.tsv",candidates)
    write_tsv(ROOT/"reports/clinical_gate_sdrf_rows.tsv",sorted(row_audit,key=lambda r:r["sdrf_line"]))
    write_tsv(ROOT/"reports/clinical_gate_discovery_gsm.tsv",discovery)
    counts=dict(collections.Counter(r["disease_original"] for r in candidates))
    summary={
        "status":"M0_OPEN_NOT_PASSED", "scope":"existing_metadata_only_before_clinical_XLSX_join",
        "n_sdrf_records":len(rows), "n_nuclei_records":len(nuclei), "n_declared_nuclei_individuals":len(donors),
        "candidate_declared_ids":len(candidates), "candidate_groups":counts,
        "candidate_records":sum(r["record_count"] for r in candidates),
        "candidate_suffix_mismatch_records":sum(r["source_suffix_mismatch_records"] for r in candidates),
        "candidate_ids_with_suffix_mismatch":sum(r["source_suffix_mismatch_records"]>0 for r in candidates),
        "nuclei_suffix_mismatch_records":sum(not r["source_suffix_matches_individual"] for r in row_audit),
        "biostudies_nuclei_matches":sum(r["biostudies_declared_record_match"] for r in row_audit),
        "candidate_repeated_ids":[r["declared_individual"] for r in candidates if r["record_count"]>1],
        "candidate_exact_age_available":0, "candidate_donor_hormone_available":0,
        "geo_samples":len(discovery), "geo_groups":dict(collections.Counter(r["disease_original"] for r in discovery)),
        "geo_characteristic_keys":sorted(set().union(*(json.loads(r["characteristics_original_json"]) for r in discovery))),
        "geo_phase_available":0, "geo_age_available":0,
        "prohibitions":["no_expression_read","no_disease_effects","no_source_suffix_reassignment","no_age_midpoint_imputation","no_hormone_inference_from_phase"],
    }
    (ROOT/"reports/clinical_gate_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    inputs=[path,metadata/"E-MTAB-14039.json",metadata/"GSE266265_family.soft.gz",ROOT/"data/reference/discovery_article.xml",ROOT/"data/reference/heca_article.xml"]
    outputs=[ROOT/"manifests/clinical_candidates.tsv",ROOT/"reports/clinical_gate_sdrf_rows.tsv",ROOT/"reports/clinical_gate_discovery_gsm.tsv",ROOT/"reports/clinical_gate_summary.json"]
    receipt={
        "generated_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "host_expected":"SRV01", "project_root":str(ROOT), "command":"python3 scripts/clinical_gate_audit.py",
        "inputs":[{"path":str(p.relative_to(ROOT)),"sha256":digest(p),"bytes":p.stat().st_size} for p in inputs],
        "outputs":[{"path":str(p.relative_to(ROOT)),"sha256":digest(p),"bytes":p.stat().st_size} for p in outputs],
        "script_sha256":digest(pathlib.Path(__file__)),
        "clinical_xlsx_join":"not_in_this_source_only_pass",
    }
    (ROOT/"receipts/clinical_gate_audit.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
if __name__ == "__main__":
    main()
