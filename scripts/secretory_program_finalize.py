#!/usr/bin/env python3
"""Finalize the secretory-phase program-level disease contrast: verify the freeze on the published
checkout and every cohort output's provenance, apply the fixed pair-level evidence label, publish
compact copies under reports/secretory_program/ and write the results receipt. Standard library
only; receipt paths are relative to the project root."""
import argparse,datetime,gzip,hashlib,json,pathlib,shutil

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def gz_copy(src,dst):
    with open(src,"rb") as f, open(dst,"wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz:shutil.copyfileobj(f,gz)

COHORTS=("discovery","heca")
FREEZE_PHASE="SECRETORY_PROGRAM_RULES_FROZEN_PRE_ANALYSIS"

def pair_label(summaries,rules):
    labels={c:summaries[c]["label"] for c in COHORTS}
    if any(l in ("NOT_INTERPRETABLE","NOT_ESTIMABLE") for l in labels.values()):return "NOT_EVALUABLE"
    m={c:summaries[c]["M1_primary"] for c in COHORTS}
    a=rules["replication_alpha_each"];b=rules["per_cohort_alpha"]
    if all(m[c]["estimate"]<0 and m[c]["p"]<a for c in COHORTS):return "replicated_attenuation"
    if all(m[c]["estimate"]>0 and m[c]["p"]<a for c in COHORTS):return "replicated_enhancement"
    n=sum(1 for c in COHORTS if m[c]["p"]<b)
    return "single_cohort_signal" if n==1 else ("no_detectable_difference" if n==0 else "single_cohort_signal_in_each_direction")

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--work",default="work/secretory_program");ap.add_argument("--git-head",required=True)
    ap.add_argument("--host-cleanup",default=None,help="optional JSON written by the host clean-up step")
    ap.add_argument("--freeze",default="receipts/secretory_program_freeze.json");a=ap.parse_args()
    root=a.root.resolve();cfg=read(root/"config/secretory_program.json");work=root/a.work
    pub=root/cfg["outputs"]["publish_dir"];rel=lambda p:str(pathlib.Path(p).resolve().relative_to(root))
    freeze=read(root/a.freeze)
    if freeze["phase"]!=FREEZE_PHASE:raise ValueError("wrong freeze")
    freeze_sha=sha(root/a.freeze)
    for p,digest in freeze["artifacts"].items():
        if sha(root/p)!=digest:raise ValueError("frozen artifact changed on the published checkout: "+p)
    if pub.exists() and any(pub.iterdir()):raise ValueError("publish directory is not empty: "+str(pub))
    pub.mkdir(parents=True,exist_ok=True)
    exempt=set(freeze.get("host_verification_exempt",[]))
    artifacts=dict(freeze["artifacts"]);published={};cohorts={};summaries={}
    for c in COHORTS:
        base=work/c;prov=read(base/"provenance.json")
        if prov["synthetic"] or prov["cohort"]!=c or prov["analysis"]!="SECRETORY_PROGRAM_DISEASE_CONTRAST":raise ValueError("provenance invalid: "+c)
        if prov["freeze"]["sha256"]!=freeze_sha:raise ValueError("output was not produced under this freeze: "+c)
        bad=[p for p,ok in prov["freeze"]["artifact_check"].items() if not ok and not p.startswith("docs/") and p not in exempt]
        if bad:raise ValueError("host freeze check failed for "+", ".join(bad))
        if not prov["identity_check"]["passed"]:raise ValueError("identity checks did not pass: "+c)
        for k,digest in prov["input_sha256"].items():
            key=[kk for kk,v in cfg["inputs"][c].items() if v==k]
            other=cfg["signature"]["source"][c]
            pinned=cfg["inputs"][c].get(key[0]+"_sha256") if key else (cfg["inputs"][other]["phase_fit_table_sha256"] if k==cfg["inputs"][other]["phase_fit_table"] else None)
            if pinned!=digest:raise ValueError("input hash on the host differs from the pin: "+k)
        for p,digest in prov["output_sha256"].items():
            q=base/p
            if not q.exists() or sha(q)!=digest:raise ValueError("output changed or missing: "+str(q))
        s=read(base/"summary.json");summaries[c]=s
        if s["cohort"]!=c or s["synthetic"]:raise ValueError("summary invalid: "+c)
        for fn in ("interaction_de_all.tsv",):
            gz_copy(base/fn,pub/f"{c}_{fn}.gz");published[f"{c}_{fn}.gz"]=rel(pub/f"{c}_{fn}.gz")
        for fn in ("donor_scores.tsv","signature_genes.tsv","leave_one_donor_out.tsv","interaction_top20_by_abs_t.tsv","identity_check.json","summary.json","provenance.json"):
            shutil.copyfile(base/fn,pub/f"{c}_{fn}");published[f"{c}_{fn}"]=rel(pub/f"{c}_{fn}")
        if (base/"scores.png").exists():shutil.copyfile(base/"scores.png",pub/f"{c}_scores.png");published[f"{c}_scores.png"]=rel(pub/f"{c}_scores.png")
        for fn in sorted(base.iterdir()):
            if fn.is_file():artifacts[rel(fn)]=sha(fn)
        m=s["M1_primary"]
        cohorts[c]={"label":s["label"],"gate":{"pass":s["validity_gate"]["pass"],"auc":s["validity_gate"]["auc"],"phase_p":s["validity_gate"]["model"].get("p")},
                    "theta":{k:m.get(k) for k in ("status","estimate","se","t","p","ci_low","ci_high","residual_df","mde_log","estimate_in_secretory_sd_units")},
                    "delta_M2":{k:s["M2_secretory"].get(k) for k in ("status","estimate","p","ci_low","ci_high","residual_df")},
                    "delta_M2_subphase":{k:s["M2_subphase"].get(k) for k in ("status","estimate","p","ci_low","ci_high","residual_df")},
                    "signature":s["signature"],"donors":s["donors"]["by_phase_disease"],"gene_level":{k:s["gene_level"].get(k) for k in ("genes_BH_q_lt_0_05","fry","interaction_projection_on_signature","mean_signed_interaction_log2fc")},
                    "leave_one_donor_out":s["leave_one_donor_out"],"provenance_sha256":sha(base/"provenance.json")}
    for v in published.values():artifacts[v]=sha(root/v)
    pair=pair_label(summaries,cfg["evidence_labels"])
    host=read(root/a.host_cleanup) if a.host_cleanup else None
    receipt={"schema_version":1,"phase":"SECRETORY_PROGRAM_ANALYSIS_COMPLETED","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
             "git_head":a.git_head,"protocol":cfg["protocol"],"freeze":{"path":a.freeze,"sha256":freeze_sha,"git_head":freeze["git_head"]},
             "artifacts":artifacts,"published_copies":published,"cohorts":cohorts,"pair_label":pair,"label_rules":cfg["evidence_labels"],"host_cleanup":host,
             "attestations":{"computed_once":True,"labels_applied_as_frozen":True,"new_disease_fits":True,"v1_artifacts_modified":False,"positive_control_artifacts_modified":False},
             "interpretation":"a program-level association in eutopic endometrium between endometriosis and the secretory stromal program; labels follow the frozen rules; no gene-level claim and no target nomination"}
    out=root/cfg["outputs"]["results_receipt"];out.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"pair_label":pair,"cohorts":{c:{"label":cohorts[c]["label"],"theta":cohorts[c]["theta"],"gate":cohorts[c]["gate"]} for c in COHORTS}},indent=2))
if __name__=="__main__":main()
