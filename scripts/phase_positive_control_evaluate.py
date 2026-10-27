#!/usr/bin/env python3
"""Phase positive control: evaluate the frozen F1 discovery lists against the HECA SEC_minus_PRO
fit with the v1 external-evaluation definitions, imported unchanged from evaluate_replication.py
(exclusion precedence, agree/evaluable/selected, missing-outcome bounds, coverage rule, G1-B0
contrast, same-direction external q<0.05). Records P1 and P2 against the protocol thresholds
and reads the HECA marker summary for S1. No dossiers and no matched null sets: those are v1
disease-analysis instruments, not positive-control endpoints. Standard library only."""
import argparse,csv,datetime,gzip,hashlib,json,pathlib,platform,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import evaluate_replication as v1

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=HERE.parent);ap.add_argument("--f1-receipt",default="receipts/phase_positive_control_f1_discovery_freeze.json")
    ap.add_argument("--external-dir",required=True);ap.add_argument("--counts",default="data/processed/phase_positive_control/heca/counts_stable.tsv.gz")
    ap.add_argument("--markers-summary",required=True);ap.add_argument("--phase-config",default="config/phase_positive_control.json")
    ap.add_argument("--output-dir",required=True);ap.add_argument("--git-head",required=True);a=ap.parse_args()
    root=a.root.resolve();ext=root/a.external_dir;out=root/a.output_dir;pcfg=json.loads((root/a.phase_config).read_text())
    if out.exists() and any(out.iterdir()):raise ValueError("Output directory must be empty; previous evaluation preserved")
    # ---- F1 freeze: lists, directions, ranks ----
    f1=json.loads((root/a.f1_receipt).read_text())
    if f1.get("phase")!="PHASE_POSITIVE_CONTROL_F1_DISCOVERY_FROZEN_PRE_HECA" or f1.get("heca_phase_effects_computed") is not False:raise ValueError("Invalid F1 freeze")
    f1_check={}
    for p,d in f1["artifacts"].items():
        ok=(root/p).exists() and sha(root/p)==d;f1_check[p]=ok
        if not ok and not p.startswith("docs/"):raise ValueError("Changed F1 artifact: "+p)
    products=f1["products"];lists={};ranks={}
    for m in v1.METHODS:
        full=v1.read_tsv(root/products[m+"_full"],("rank","gene_id","method","discovery_sign","logFC"))
        for rank,row in enumerate(full,1):
            if int(row["rank"])!=rank or row["method"]!=m or int(row["discovery_sign"])!=v1.sign(row["logFC"]):raise ValueError("Invalid frozen ranking or direction: "+m)
        ranks[m]={row["gene_id"]:int(row["rank"]) for row in full}
        for k in (20,50):
            top=v1.read_tsv(root/products[f"{m}_top{k}"])
            if top!=full[:k] or len(top)!=k:raise ValueError("Frozen list is not the exact full-rank prefix")
            lists[(m,k)]=top
    # ---- external (HECA) fit ----
    prov=json.loads((ext/"provenance.json").read_text())
    if prov.get("cohort")!="heca" or prov.get("contrast")!="SEC_minus_PRO" or prov.get("data_kind")!="public_cohort" or prov.get("scientific_effect_run") is not True:raise ValueError("External provenance is not the HECA phase fit")
    if prov.get("f1_freeze_phase")!="PHASE_POSITIVE_CONTROL_F1_DISCOVERY_FROZEN_PRE_HECA":raise ValueError("HECA fit did not verify the F1 freeze")
    outputs={}
    for original,digest in prov["output_sha256"].items():
        name=pathlib.Path(original).name
        if name in outputs or sha(ext/name)!=digest:raise ValueError("External output hash mismatch: "+original)
        outputs[name]=digest
    filter_rows=v1.read_tsv(ext/"expression_filter.tsv",("gene_id","expression_pass"))
    filters={g:v1.boolean(r["expression_pass"]) for g,r in v1.index_unique(filter_rows).items()}
    with gzip.open(root/a.counts,"rt") as f:
        reader=csv.reader(f,delimiter="\t");header=next(reader);counts_ids=[r[0] for r in reader]
    if header[0]!="gene_id" or counts_ids!=[r["gene_id"] for r in filter_rows]:raise ValueError("Expression filter ledger does not cover the exact count IDs/order")
    if len(counts_ids)!=prov["audit"]["input_genes"] or len(header)-1!=prov["audit"]["n_donors"]:raise ValueError("Counts dimensions disagree with external provenance")
    family=v1.read_tsv(ext/"de_all.tsv",v1.DE_COLUMNS);fitted=v1.read_tsv(ext/"de_fitted.tsv",v1.DE_COLUMNS)
    external=v1.validate_family(filters,family,fitted)
    if len(family)!=prov["tested_family_size"] or len(fitted)!=prov["estimable_genes"]:raise ValueError("Family/fit counts disagree with provenance")
    # ---- endpoints ----
    details=[];summaries={}
    for m in v1.METHODS:
        summaries[m]={}
        for k in (20,50):
            records,summary=v1.evaluate_set(lists[(m,k)],filters,external,k)
            details.extend({**r,"method":m,"k":k} for r in records);summaries[m][str(k)]=summary
    primary=v1.primary_contrast(summaries["B0"]["20"],summaries["G1"]["20"])
    ep=pcfg["endpoints"];b0=summaries["B0"]["20"]
    p2={"list":"frozen discovery B0 top-20","same_sign":b0["agree"],"evaluable":b0["evaluable"],"selected":20,"same_sign_and_external_q_lt_0_05":b0["significant_same_direction"],"missing_bounds":b0["missing_bounds"],"coverage_limited":b0["evaluable"]<16,"thresholds":ep["P2"]["pass"],"pass":bool(b0["evaluable"]>=16 and b0["agree"]>=ep["P2"]["pass"]["same_sign_min"] and b0["significant_same_direction"]>=ep["P2"]["pass"]["same_sign_and_external_q_lt_0_05_min"])}
    disc_sig=f1["summary"]["BH_q_lt_0_05"];heca_sig=prov["significant_BH05"]
    p1={"discovery_BH_q_lt_0_05":disc_sig,"heca_BH_q_lt_0_05":heca_sig,"heca_tested_family":prov["tested_family_size"],"thresholds":ep["P1"]["pass"],"pass":bool(disc_sig>=ep["P1"]["pass"]["discovery_min"] and heca_sig>=ep["P1"]["pass"]["heca_min"]),"strong_label_min":ep["P1"]["strong_label_min"],"strong":bool(disc_sig>=ep["P1"]["strong_label_min"] and heca_sig>=ep["P1"]["strong_label_min"])}
    ms=json.loads((root/a.markers_summary).read_text())
    if ms.get("cohort")!="heca":raise ValueError("Marker summary is not the HECA readout")
    s1={"discovery":f1["summary"]["markers_discovery"],"heca":{k:ms[k] for k in ("panel_genes","in_expression_family","expected_sign_agrees","expected_sign_agrees_and_q_lt_0_05","threshold","pass")},"pass":bool(f1["summary"]["markers_discovery"]["pass"] and ms["pass"]),"role":"explanatory, not decisive"}
    if p1["pass"] and p2["pass"]:outcome={"code":"P1_and_P2_pass","statement":"The v1 pipeline detects and transfers a known stromal program at these sample sizes; the v1 disease null is not attributable to pipeline insensitivity. No gene becomes a target because of this result."}
    elif p1["pass"]:outcome={"code":"P1_pass_P2_fail","statement":"Within-cohort detection works in both cohorts; cross-platform transfer of the frozen top-20 is limited, so the v1 transfer endpoint had bounded sensitivity."}
    else:outcome={"code":"P1_fail","statement":"A pipeline, eligibility or count-preparation problem is presumed until explained; no further disease analysis on this pipeline before the cause is documented."}
    out.mkdir(parents=True,exist_ok=True)
    v1.write_tsv(out/"selected_genes.tsv",details,v1.DETAIL_FIELDS);v1.write_tsv(out/"bootstrap_interface.tsv",details,v1.DETAIL_FIELDS)
    summary={"schema_version":1,"cohort":"heca","contrast":"SEC_minus_PRO","protocol":"docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md","definitions":"docs/REPLICATION_EVALUATION.md (v1), applied unchanged","methods":summaries,"primary":primary,"P1":p1,"P2":p2,"S1":s1,"S2":{"status":"pending_paired_bootstrap","affects_pass_fail":False},"outcome":outcome,"top50_role":"secondary","dossiers_and_matched_nulls":"not part of the positive control","scientific_scope":"cross-study transcriptomic reproducibility of a cycle-phase program; not disease evidence, not therapeutic validation"}
    v1.write_json(out/"summary.json",summary)
    provenance={"schema_version":1,"completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"python":platform.python_version(),"git_head":a.git_head,"implementation_sha256":sha(pathlib.Path(__file__)),"v1_definitions_sha256":sha(HERE/"evaluate_replication.py"),"f1_receipt_sha256":sha(root/a.f1_receipt),"f1_artifact_check":f1_check,"external_provenance_sha256":sha(ext/"provenance.json"),"external_outputs":outputs,"markers_summary_sha256":sha(root/a.markers_summary),"output_sha256":{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
    v1.write_json(out/"provenance.json",provenance)
    print(json.dumps({"status":"PHASE_POSITIVE_CONTROL_EVALUATION_COMPLETE","P1":p1,"P2":p2,"S1":s1,"primary":primary,"outcome":outcome,"methods_top20":{m:summaries[m]["20"] for m in v1.METHODS}},indent=2))
if __name__=="__main__":main()
