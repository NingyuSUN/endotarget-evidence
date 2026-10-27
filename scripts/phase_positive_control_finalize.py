#!/usr/bin/env python3
"""Phase positive control: final evaluation receipt. Verifies and hashes the HECA phase fit, the
evaluation, the HECA marker readout and (when present) the S2 bootstrap, applies the protocol's
fixed interpretation rules, and publishes compact copies into the repository tree.
Standard library only; receipt paths are relative to the project root."""
import argparse,datetime,gzip,hashlib,json,pathlib,shutil

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def gz_copy(src,dst):
    with open(src,"rb") as f, open(dst,"wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz:shutil.copyfileobj(f,gz)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--heca-dir",required=True);ap.add_argument("--evaluation-dir",required=True)
    ap.add_argument("--markers",required=True);ap.add_argument("--marker-summary",required=True)
    ap.add_argument("--bootstrap-dir",default=None);ap.add_argument("--git-head",required=True)
    ap.add_argument("--f1-receipt",default="receipts/phase_positive_control_f1_discovery_freeze.json")
    ap.add_argument("--receipt",default="receipts/phase_positive_control_evaluation.json");ap.add_argument("--publish-dir",default="reports/phase_positive_control_f2");a=ap.parse_args()
    root=a.root.resolve();hd=root/a.heca_dir;ev=root/a.evaluation_dir;pub=root/a.publish_dir;rel=lambda p:str(pathlib.Path(p).resolve().relative_to(root))
    hp=read(hd/"provenance.json");ep=read(ev/"provenance.json");summary=read(ev/"summary.json");ms=read(root/a.marker_summary);f1=read(root/a.f1_receipt)
    if hp.get("cohort")!="heca" or hp.get("contrast")!="SEC_minus_PRO" or hp.get("f1_freeze_phase")!="PHASE_POSITIVE_CONTROL_F1_DISCOVERY_FROZEN_PRE_HECA":raise ValueError("HECA provenance invalid")
    for original,digest in hp["output_sha256"].items():
        if sha(hd/pathlib.Path(original).name)!=digest:raise ValueError("Changed HECA output: "+original)
    for name,digest in ep["output_sha256"].items():
        if sha(ev/name)!=digest:raise ValueError("Changed evaluation output: "+name)
    if ep["external_provenance_sha256"]!=sha(hd/"provenance.json"):raise ValueError("Evaluation did not consume this HECA fit")
    if ms.get("cohort")!="heca" or ep["markers_summary_sha256"]!=sha(root/a.marker_summary):raise ValueError("Marker summary mismatch")
    s2={"status":"not_run","affects_pass_fail":False}
    if a.bootstrap_dir:
        bd=root/a.bootstrap_dir;bs=read(bd/"summary.json");proc=read(bd/"process.json")
        if proc.get("status")!="COMPLETED":raise ValueError("Bootstrap did not complete")
        fixed=summary["primary"]["concordance_difference"]
        s2={"status":"completed","draws_required":bs["total_required"],"valid":bs["valid"],"invalid":bs["invalid"],"interval_suppressed":bs["interval_suppressed"],"conditional_95_percent_interval_G1_minus_B0":bs.get("interval"),"fixed_sample_difference":fixed,"failure_reasons":bs.get("failure_reasons"),"coverage_adequate_both":summary["primary"]["coverage_adequate_both"],"affects_pass_fail":False}
        iv=bs.get("interval")
        s2["direction_statement"]=("interval suppressed; no statement" if bs["interval_suppressed"] else "G1 below B0 across the interval" if iv[1]<0 else "G1 above B0 across the interval" if iv[0]>0 else "interval includes zero; no superiority statement")
        if not summary["primary"]["coverage_adequate_both"]:s2["direction_statement"]="coverage-limited; no superiority statement"
    if pub.exists() and any(pub.iterdir()):raise ValueError("Publish directory is not empty: "+str(pub))
    pub.mkdir(parents=True,exist_ok=True);published={}
    gz_copy(hd/"de_all.tsv",pub/"heca_de_all.tsv.gz");published["heca_de_all.tsv.gz"]=rel(pub/"heca_de_all.tsv.gz")
    for name in ["expression_filter.tsv","normalization.tsv","design.tsv","provenance.json"]:shutil.copyfile(hd/name,pub/("heca_de_"+name));published["heca_de_"+name]=rel(pub/("heca_de_"+name))
    for name in ["selected_genes.tsv","summary.json","provenance.json"]:shutil.copyfile(ev/name,pub/("evaluation_"+name));published["evaluation_"+name]=rel(pub/("evaluation_"+name))
    shutil.copyfile(root/a.markers,pub/"markers_heca.tsv");shutil.copyfile(root/a.marker_summary,pub/"markers_heca_summary.json");published["markers_heca.tsv"]=rel(pub/"markers_heca.tsv");published["markers_heca_summary.json"]=rel(pub/"markers_heca_summary.json")
    if a.bootstrap_dir:
        for name in ["paired_differences.tsv","summary.json","draw_freeze.json"]:shutil.copyfile(root/a.bootstrap_dir/name,pub/("bootstrap_"+name));published["bootstrap_"+name]=rel(pub/("bootstrap_"+name))
    products={"heca_de_all":rel(hd/"de_all.tsv"),"heca_de_fitted":rel(hd/"de_fitted.tsv"),"heca_expression_filter":rel(hd/"expression_filter.tsv"),"heca_provenance":rel(hd/"provenance.json"),"evaluation_selected_genes":rel(ev/"selected_genes.tsv"),"evaluation_summary":rel(ev/"summary.json"),"markers_heca":a.markers,"markers_heca_summary":a.marker_summary}
    if a.bootstrap_dir:products.update({"bootstrap_summary":rel(root/a.bootstrap_dir/"summary.json"),"bootstrap_paired_differences":rel(root/a.bootstrap_dir/"paired_differences.tsv"),"bootstrap_draw_ledger":rel(root/a.bootstrap_dir/"draw_ledger.tsv")})
    inputs={p:d for p,d in hp["input_sha256"].items()}
    for p in ["scripts/phase_positive_control_fit_cohort.R","scripts/phase_positive_control_evaluate.py","scripts/evaluate_replication.py","scripts/phase_positive_control_markers.py","scripts/phase_positive_control_bootstrap.R","scripts/resample_donor_models.R","scripts/phase_positive_control_finalize.py","config/phase_positive_control.json",a.f1_receipt]:inputs[p]=sha(root/p)
    artifacts={**inputs,**{p:sha(root/p) for p in products.values()},**{v:sha(root/v) for v in published.values()}}
    v1_before={"data/processed/discovery/counts_stable.tsv.gz":sha(root/"data/processed/discovery/counts_stable.tsv.gz"),"data/processed/heca/counts_stable.tsv.gz":sha(root/"data/processed/heca/counts_stable.tsv.gz")}
    receipt={"schema_version":1,"phase":"PHASE_POSITIVE_CONTROL_EVALUATED","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"git_head":a.git_head,"protocol":"docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md","f0_freeze":"receipts/phase_positive_control_freeze.json","f1_freeze":a.f1_receipt,"artifacts":artifacts,"products":products,"published_copies":published,
      "heca_fit":{"tested_family_size":hp["tested_family_size"],"estimable_genes":hp["estimable_genes"],"BH_q_lt_0_05":hp["significant_BH05"],"BH_q_lt_0_05_up":hp.get("significant_BH05_up"),"BH_q_lt_0_05_down":hp.get("significant_BH05_down"),"design":hp["audit"],"filter_audit":hp["filter_audit"],"captured_warnings":hp.get("captured_warnings"),"f0_artifact_check_on_analysis_host":hp.get("f0_artifact_check"),"f1_artifact_check_on_analysis_host":hp.get("f1_artifact_check")},
      "endpoints":{"P1":summary["P1"],"P2":summary["P2"],"S1":summary["S1"],"S2":s2,"primary_contrast":summary["primary"],"methods":summary["methods"]},"outcome":summary["outcome"],"v1_matrices_sha256_at_evaluation":v1_before,"v1_artifacts_written":False,"heca_phase_effects_computed":True,"software":{"fit":hp.get("software"),"evaluation_python":ep.get("python")},"interpretation":"cross-study transcriptomic reproducibility of a cycle-phase program as a pipeline positive control; not disease evidence, not targets, not therapeutic validation"}
    out=root/a.receipt;out.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ("phase","git_head","heca_fit","endpoints","outcome")},indent=2,default=str))
if __name__=="__main__":main()
