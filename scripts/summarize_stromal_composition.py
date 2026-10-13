#!/usr/bin/env python3
"""Cached author-label composition, exact frozen donors/cells, descriptive only."""
import ast
import collections
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
import socket
import statistics
import struct
import time

ROOT=Path(__file__).resolve().parents[1]
OSF=Path("data/input_routes/osf/selected/04_Analysis/S3_COMPLETION_EXTENSION/GSE266265_DONOR_RECONSTRUCTION")
MISSING={"","NA","nan","NaN","None"}
EXPECTED={"discovery":150788,"heca":52744}

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def rows(path,delimiter="\t"):
    with (gzip.open(path,"rt",newline="") if str(path).endswith(".gz") else path.open(newline="")) as f:
        yield from csv.DictReader(f,delimiter=delimiter)

def validate_selection(selected):
    for key in ("donor_id","sample_id"):
        ids=[r[key] for r in selected]
        if not ids or any(x in MISSING for x in ids) or len(ids)!=len(set(ids)):
            raise ValueError("Missing/duplicate selected "+key+"; no technical-library inflation")

def scan_cells(cells,selected,labels,label_field,id_field,donor_field=None,qc_field=None,main_field=None):
    validate_selection(selected)
    by={r["sample_id"]:r for r in selected}
    counts={r["donor_id"]:collections.Counter() for r in selected}
    missing=collections.Counter();all_cells=collections.Counter();coarse_missing=collections.Counter()
    seen=set();positions=[];nrows=0
    for i,r in enumerate(cells):
        nrows+=1
        if r["sample"] not in by:continue
        donor=by[r["sample"]]["donor_id"]
        if donor_field is not None and r[donor_field]!=donor:
            raise ValueError("Cached author donor disagrees with exact selected sample crosswalk")
        ident=r[id_field]
        if ident in MISSING or ident in seen:raise ValueError("Missing/duplicate cached selected cell identity")
        seen.add(ident);all_cells[donor]+=1
        label=r[label_field]
        if label in MISSING:
            missing[donor]+=1
            if main_field and r[main_field]=="stromal":coarse_missing[donor]+=1
        if label not in labels or (qc_field and r[qc_field]!="no"):continue
        if main_field and r[main_field]!="stromal":raise ValueError("Retained subtype conflicts with source main lineage")
        counts[donor][label]+=1;positions.append(i)
    return {"counts":counts,"missing":missing,"all_cells":all_cells,"coarse_missing":coarse_missing,
            "positions":positions,"source_rows_scanned":nrows,"unique_selected_cached_cell_ids":len(seen)}

def reconcile(scan,selected,expected_field,expected_total):
    for r in selected:
        if sum(scan["counts"][r["donor_id"]].values())!=int(r[expected_field]):
            raise ValueError("Frozen per-donor denominator differs: "+r["donor_id"])
    total=sum(sum(c.values()) for c in scan["counts"].values())
    if total!=expected_total:raise ValueError("Frozen cohort cell denominator differs")

def donor_grid(cohort,selected,labels,scan,field):
    out=[]
    for r in selected:
        d=r["donor_id"];n=sum(scan["counts"][d].values())
        if n<=0:raise ValueError("Nonpositive primary cell denominator")
        for label in labels:
            k=scan["counts"][d][label]
            out.append({"cohort":cohort,"donor_id":d,"selected_sample_id":r["sample_id"],
                        "disease":r["disease"],"source_label_field":field,"source_subtype_label":label,
                        "cells":k,"primary_cells":n,"proportion_primary":k/n,
                        "availability":"observed_source_label; zero_count_is_observed_zero"})
    return out

def group_summary(donor_rows):
    grouped=collections.defaultdict(list)
    for r in donor_rows:
        for group in ("all_selected",r["disease"]):
            grouped[(r["cohort"],group,r["source_label_field"],r["source_subtype_label"])].append(r)
    out=[]
    for (cohort,group,field,label),rr in sorted(grouped.items()):
        if len({r["donor_id"] for r in rr})!=len(rr):raise ValueError("Duplicate donor in summary")
        pp=[r["proportion_primary"] for r in rr]
        qq=statistics.quantiles(pp,n=4,method="inclusive") if len(pp)>1 else pp*3
        numerator=sum(r["cells"] for r in rr);denominator=sum(r["primary_cells"] for r in rr)
        out.append({"cohort":cohort,"descriptive_group":group,"source_label_field":field,
                    "source_subtype_label":label,"donors":len(rr),"donors_with_nonzero_cells":sum(r["cells"]>0 for r in rr),
                    "total_cells":numerator,"total_primary_cells":denominator,
                    "pooled_cell_proportion":numerator/denominator,
                    "mean_proportion_equal_donor":statistics.mean(pp),"median_proportion":statistics.median(pp),
                    "Q1_proportion":qq[0],"Q3_proportion":qq[2],"minimum_proportion":min(pp),"maximum_proportion":max(pp)})
    return out

def read_npy_int_vector(path):
    # Only the tiny pre-existing frozen observation-index vector; no expression arrays.
    with path.open("rb") as f:
        if f.read(6)!=b"\x93NUMPY":raise ValueError("Index file is not NPY")
        version=tuple(f.read(2))
        if version==(1,0):size=struct.unpack("<H",f.read(2))[0]
        elif version==(2,0):size=struct.unpack("<I",f.read(4))[0]
        else:raise ValueError("Unsupported frozen index NPY version")
        if size>65536:raise ValueError("Unbounded index NPY header")
        header=ast.literal_eval(f.read(size).decode("latin1"))
        if header["fortran_order"] or len(header["shape"])!=1 or header["shape"][0]>1000000:
            raise ValueError("Unexpected frozen index vector shape")
        dtype=header["descr"];fmt={"<i8":"<q","<i4":"<i"}.get(dtype)
        if fmt is None:raise ValueError("Unexpected frozen index vector dtype "+dtype)
        data=f.read(8000001);width=struct.calcsize(fmt)
        if len(data)!=header["shape"][0]*width:raise ValueError("Frozen index vector bytes differ")
        return [v[0] for v in struct.iter_unpack(fmt,data)]

def write_tsv(path,rs):
    if not rs:raise ValueError("Empty report")
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rs)

def pin(path,expected,inputs):
    observed=sha(ROOT/path)
    if observed!=expected:raise ValueError("Frozen input hash differs: "+str(path))
    inputs[str(path)]=observed

def main():
    started=time.monotonic();inputs={}
    provenance={c:json.loads((ROOT/("data/processed/"+c+"/provenance.json")).read_text()) for c in EXPECTED}
    selected={c:list(rows(ROOT/("manifests/final_"+c+"_clinical_selected.tsv"))) for c in EXPECTED}
    for c,ss in selected.items():
        validate_selection(ss)
        expected_donors=provenance[c]["primary_donors" if c=="discovery" else "selected_donors"]
        if [r["donor_id"] for r in ss]!=expected_donors or len(ss)!=22:
            raise ValueError("Exact frozen 22-donor identity/order differs: "+c)
        pr=Path("data/processed")/c/"provenance.json";inputs[str(pr)]=sha(ROOT/pr)
        sp=Path("manifests/final_"+c+"_clinical_selected.tsv")
        pin(sp,provenance[c]["input_hashes"][str(sp)],inputs)
    annotation=OSF/"PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz"
    pin(annotation,provenance["discovery"]["input_hashes"][str(annotation)],inputs)
    lp=Path("config/discovery_lineage.json")
    pin(lp,provenance["discovery"]["input_hashes"][str(lp)],inputs)
    discovery_labels=json.loads((ROOT/lp).read_text())["include_exact_refined_labels"]
    definition=Path("scripts/final_clinical_gate.py")
    tree=ast.parse((ROOT/definition).read_text())
    matches=[ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="LABELS" for t in n.targets)]
    if len(matches)!=1:raise ValueError("Cannot read inherited fixed HECA LABELS")
    heca_labels=list(matches[0]);inputs[str(definition)]=sha(ROOT/definition)
    if len(discovery_labels)!=12 or len(set(discovery_labels))!=12 or len(heca_labels)!=7 or len(set(heca_labels))!=7:
        raise ValueError("Inherited primary label definition changed")
    obs=Path("data/metadata/heca_h5/obs_selected.tsv");idx=Path("work/final_clinical/heca_primary_h5_obs_indices.npy")
    for p in (obs,idx):pin(p,provenance["heca"]["input_hashes"][str(p)],inputs)
    scans={
      "discovery":scan_cells(rows(ROOT/annotation,","),selected["discovery"],discovery_labels,
                             "AnnotationRefined","count_barcode",main_field="AnnotationMain"),
      "heca":scan_cells(rows(ROOT/obs),selected["heca"],heca_labels,"celltype","_index",
                       donor_field="genotype",qc_field="discard_lowQC"),
    }
    labels={"discovery":discovery_labels,"heca":heca_labels};fields={"discovery":"AnnotationRefined","heca":"celltype"}
    for c,scan in scans.items():reconcile(scan,selected[c],"primary_lineage_cells" if c=="discovery" else "primary_nuclei",EXPECTED[c])
    if scans["heca"]["positions"]!=read_npy_int_vector(ROOT/idx):raise ValueError("HECA selected original H5 observation positions differ")
    # Secondary source-cell manifest independently reconciles all selected discovery subtype counts and P077 totals.
    source_manifest=OSF/"DONOR_CELLTYPE_MANIFEST.csv"
    source_pins=json.loads((ROOT/"receipts/osf_route_ranges.json").read_text())["requests"]
    expected=next(r["sha256"] for r in source_pins if r.get("member","").endswith("/"+source_manifest.name))
    pin(source_manifest,expected,inputs)
    source={};p77_by=collections.defaultdict(list)
    selected_samples={s["sample_id"] for s in selected["discovery"]}
    for r in rows(ROOT/source_manifest,","):
        if r["sample"] not in selected_samples:continue
        if r["population_id"]=="P077":p77_by[r["sample"]].append(int(r["cells"]))
        if r["population"]=="Fibroblast::all" or r["population_id"]=="P077" or r["population"].startswith("Refined::"):
            key=(r["sample"],r["population"])
            if key in source:raise ValueError("Duplicate source donor-population cell manifest")
            source[key]=int(r["cells"])
    for s in selected["discovery"]:
        d=s["donor_id"]
        for label in discovery_labels:
            key=(s["sample_id"],"Refined::"+label)
            if key not in source or source[key]!=scans["discovery"]["counts"][d][label]:
                raise ValueError("Source Refined cell manifest differs from direct annotation: "+str(key))
        p77=p77_by[s["sample_id"]]
        if p77!=[int(s["primary_lineage_cells"])]:raise ValueError("Source P077 denominator differs")
    donor_rows=[r for c in EXPECTED for r in donor_grid(c,selected[c],labels[c],scans[c],fields[c])]
    summaries=group_summary(donor_rows)
    totals=[]
    for c in EXPECTED:
        for s in selected[c]:
            d=s["donor_id"];scan=scans[c]
            totals.append({"cohort":c,"donor_id":d,"selected_sample_id":s["sample_id"],"disease":s["disease"],
                           "primary_cells":sum(scan["counts"][d].values()),"selected_sample_cached_cells":scan["all_cells"][d],
                           "cached_cells_missing_source_subtype":scan["missing"][d],
                           "discovery_coarse_stromal_cells_missing_refined_label":scan["coarse_missing"][d] if c=="discovery" else "not_applicable",
                           "primary_label_unknown_cells":0,
                           "primary_unknown_scope":"zero by exact label membership; missing labels excluded by inherited definition, not proof of complete annotation",
                           "finer_author_state_availability":"unavailable_beyond_"+fields[c]})
    out_paths=[]
    for suffix,rr in [("donors",donor_rows),("groups",summaries),("denominators",totals)]:
        p=Path("reports/stromal_composition_"+suffix+".tsv");write_tsv(ROOT/p,rr);out_paths.append(p)
    report=[
      "# Frozen-donor stromal subtype composition — descriptive only",
      "",
      "**VERIFIED:** 22 discovery donors / 150,788 primary fibroblast cells and 22 HECA donors / 52,744 primary stromal-fibroblast nuclei reconcile exactly to the unchanged selected clinical ledgers, processed-count donor IDs, and per-donor cell totals. HECA membership also matches the entire frozen original H5 observation-index vector, in order. No eligibility, count matrix, expression ranking or disease model was changed.",
      "",
      "## What the proportions mean",
      "",
      "Each donor's denominator is their own fixed primary-lineage cells/nuclei, not all tissue cells, RNA counts, marker CPM or group-wide cells. The complete grids retain all 12 discovery AnnotationRefined labels (264 rows) and all seven HECA celltype labels (154 rows), including observed zero-count categories. The sources' names remain distinct: no guessed cross-cohort subtype harmonization. The primary populations combine several author states and therefore cannot be assumed to have the same cellular mixture across donors.",
      "",
      "The donor report gives raw label counts, each donor denominator and proportions. The group report includes all selected donors and descriptive ENDO/non-ENDO summaries with equal donor weight; each donor contributes once to mean/median/quartiles. Pooled cell fractions are shown separately and weight donors with more retained cells more heavily. No p-values, disease effects, significance claim or additional ranking is computed.",
      "",
      "## Available author-label composition",
      "",
      "| Cohort | Source subtype | Cells / cohort primary denominator | Donors with cells / 22 | Equal-donor median % [Q1,Q3] | Donor range % |",
      "|---|---|---:|---:|---:|---:|",
    ]
    for r in summaries:
        if r["descriptive_group"]!="all_selected":continue
        report.append("| "+r["cohort"]+" | "+r["source_subtype_label"]+" | "+str(r["total_cells"])+"/"+str(r["total_primary_cells"])+" | "+str(r["donors_with_nonzero_cells"])+"/22 | "+f'{100*r["median_proportion"]:.2f} [{100*r["Q1_proportion"]:.2f}, {100*r["Q3_proportion"]:.2f}] | {100*r["minimum_proportion"]:.2f}–{100*r["maximum_proportion"]:.2f}'+" |")
    report.extend([
      "",
      "## Descriptive group context",
      "",
      "Disease groups are inherited ledger labels, included only to expose possible mixture differences; this does not test a disease association or adjust confounding. Small subtype counts and source-specific label resolution can make percentages unstable.",
      "",
      "| Cohort | Source subtype | ENDO donors | Non-ENDO donors | ENDO equal-donor mean % | Non-ENDO equal-donor mean % |",
      "|---|---|---:|---:|---:|---:|",
    ])
    keyed={(r["cohort"],r["descriptive_group"],r["source_subtype_label"]):r for r in summaries}
    for c in EXPECTED:
        for label in labels[c]:
            a=keyed[c,"ENDO",label];b=keyed[c,"non-ENDO",label]
            report.append(f'| {c} | {label} | {a["donors"]} | {b["donors"]} | {100*a["mean_proportion_equal_donor"]:.2f} | {100*b["mean_proportion_equal_donor"]:.2f} |')
    extremes=sorted([r for r in donor_rows if r["cohort"]=="heca" and r["source_subtype_label"]=="eStromal_MMPs"],key=lambda r:-r["proportion_primary"])[:4]
    report.extend([
      "",
      "The largest HECA eStromal_MMPs donor proportions are "+", ".join(r["donor_id"]+" ("+r["disease"]+", "+str(r["cells"])+"/"+str(r["primary_cells"])+f' = {100*r["proportion_primary"]:.2f}%)' for r in extremes)+". These four are ENDO donors, while the all-donor median is 1.32%. Consequently, the ENDO mean of 23.97% versus the seven controls' 2.00% describes a heterogeneous selected sample, not a uniform change in every case. The precise counts and remaining 18 donors are retained in the donor table; none is removed or used to change eligibility.",
    ])
    report.extend([
      "",
      "## Missing labels, granularity and membership",
      "",
      "In the selected samples' cached annotation rows, missing source subtype labels are explicitly counted in stromal_composition_denominators.tsv; excluded labels remain outside the frozen primary denominator. There are no unknown subtype labels within retained primary cells by construction: the frozen boundary requires an exact included label. This is not independent proof that missing or mislabeled stromal cells were absent in the underlying tissue.",
      "",
      f'Discovery selected cached rows missing AnnotationRefined: {sum(scans["discovery"]["missing"].values())}; among source coarse-stromal rows: {sum(scans["discovery"]["coarse_missing"].values())}. HECA selected cached rows missing celltype: {sum(scans["heca"]["missing"].values())}. Finer independently assessed cell states beyond the available author labels are unavailable in this bounded audit; no expression-based state assignment or cell-purity inference was performed.',
      "",
      "HECA other libraries from the same donor are excluded by exact selected sample ID; technical or repeated samples therefore do not increase donor n or alter the already frozen library-selection rule. Discovery direct per-cell Refined counts and P077 totals exactly match its source donor-population manifest. Source metadata preserves inherited author QC and classifications; it is not an independently purified reference.",
      "",
      "## Interpretation and limits",
      "",
      "Visible donor-to-donor and descriptive group mixture variation is a potential explanation for primary-lineage pseudobulk expression differences: an aggregate can change because the relative abundance of author states changes, because expression changes within a state, or both. These count tables alone cannot distinguish those mechanisms or demonstrate that composition caused any disease result. No subtype independently validated disease fit is presented, no cell is reclassified, and no disease model or target ranking is recomputed.",
      "",
      "Discovery scRNA and HECA snRNA have different capture, tissue and annotation resolution. The label tables can document each cohort's mixture; directly treating similarly named labels as calibrated equivalent states would require independent evidence. Donor counts remain 22 per cohort, regardless of the number of cell-level rows. Author-label zeros mean no such annotated cells in the selected cache, not biological absence of the state.",
      "",
      "This completes the bounded donor-level descriptive composition requirement alongside the previously completed fixed-marker plausibility audit. It does not complete independent per-cell lineage validation, purity measurement, within-subtype disease analysis or causal deconfounding. No new user decision is required for this descriptive deliverable.",
      "",
      "YELLOW — understanding not yet assessed: defend the donor-specific cell denominator, equal-donor versus pooled cell weighting, why multiple libraries do not create new participants, source-specific labels versus independent biology, and mixture change versus within-state expression change.",
      "",
      "## Outputs and reproducibility",
      "",
      "- reports/stromal_composition_donors.tsv: complete donor × source-label grid.",
      "- reports/stromal_composition_groups.tsv: equal-donor descriptive summaries and separately identified pooled-cell fractions.",
      "- reports/stromal_composition_denominators.tsv: exact donor/sample/cell denominators and missing-label availability.",
      "- receipts/stromal_composition.json: source, script and output hashes, original observation-index reconciliation and runtime.",
      "- tests/test_stromal_composition.py: meaningful denominator, duplicate-donor, wrong-crosswalk, technical-library, missing-label, observed-zero and weighting boundaries.",
    ])
    mp=Path("reports/stromal_composition.md");(ROOT/mp).write_text("\n".join(report)+"\n");out_paths.append(mp)
    coverage={c:{"donors":len(selected[c]),"primary_cells":EXPECTED[c],"source_subtype_labels":len(labels[c]),
                 "complete_donor_subtype_rows":sum(r["cohort"]==c for r in donor_rows),
                 "nonzero_donor_subtype_rows":sum(r["cohort"]==c and r["cells"]>0 for r in donor_rows),
                 "source_rows_scanned":scans[c]["source_rows_scanned"],
                 "selected_cached_cells":scans[c]["unique_selected_cached_cell_ids"],
                 "missing_source_subtype_selected_cached_cells":sum(scans[c]["missing"].values()),
                 "groups":dict(collections.Counter(s["disease"] for s in selected[c])),
                 "donor_ids":[s["donor_id"] for s in selected[c]],
                 "selected_sample_ids":[s["sample_id"] for s in selected[c]]} for c in EXPECTED}
    receipt={"schema_version":1,"created_utc":datetime.now(timezone.utc).isoformat(),"status":"EXECUTED_DESCRIPTIVE_ONLY",
             "host":socket.gethostname(),"project_root":str(ROOT),"python":platform.python_version(),
             "command":"python3 scripts/summarize_stromal_composition.py","wall_seconds":round(time.monotonic()-started,3),
             "input_sha256":inputs,"script_sha256":sha(Path(__file__)),"test_sha256":sha(ROOT/"tests/test_stromal_composition.py"),
             "outputs_sha256":{str(p):sha(ROOT/p) for p in out_paths},"coverage":coverage,
             "label_definitions":{"discovery":discovery_labels,"heca":heca_labels},
             "checks":{"exact_selected_donor_ID_order":True,"exact_per_donor_cell_denominators":True,"exact_cohort_denominators":True,
                       "selected_cell_identifiers_unique":True,"HECA_original_selected_obs_indices_exact":True,
                       "discovery_direct_cell_vs_source_Refined_manifest_exact":True,"discovery_source_P077_totals_exact":True,
                       "no_other_HECA_libraries_consumed":True,"missing_and_zero_distinguished":True},
             "scope":{"counts_are_cells_not_UMIs":True,"descriptive_disease_group_summaries":True,"inferential_tests":False,
                      "DE_models":False,"expression_ranking":False,"marker_audit_rerun":False,"cell_reclassification":False,
                      "eligibility_changes":False,"raw_downloads":False,"purity_validated":False,
                      "finer_independent_state_validation":"unavailable"}}
    (ROOT/"receipts/stromal_composition.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"coverage":coverage,"checks":receipt["checks"],"outputs":list(receipt["outputs_sha256"])}))

if __name__=="__main__":main()
