#!/usr/bin/env python3
"""Finalize the minimum detectable effect analysis: verify the freeze and every output's
provenance, compute the prespecified bridge statistics between the phase program and the
disease MDE, publish compact copies under reports/mde_analysis/ and write the results
receipt. Standard library only; receipt paths are relative to the project root."""
import argparse,csv,datetime,gzip,hashlib,json,pathlib,shutil,statistics

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def rows(p):
    with open(p,newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def num(v):
    try:return float(v)
    except (TypeError,ValueError):return None
def gz_copy(src,dst):
    with open(src,"rb") as f, open(dst,"wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz:shutil.copyfileobj(f,gz)

FITS=("v1_discovery","v1_heca","phase_discovery","phase_heca");COHORTS=("discovery","heca")

def verify_outputs(prov,base):
    for p,digest in prov["output_sha256"].items():
        q=base/pathlib.Path(p).name
        if not q.exists() or sha(q)!=digest:raise ValueError("output changed or missing: "+p)

def bridge(disease_rows,phase_rows,scenarios):
    est=lambda r:r["fit_estimable"]=="TRUE"
    d={r["gene_id"]:r for r in disease_rows if est(r)}
    program=[r for r in phase_rows if est(r) and num(r["adj.P.Val"]) is not None and num(r["adj.P.Val"])<0.05]
    out={"phase_program_genes":len(program),"disease_family_estimable":len(d),
         "median_abs_log2fc_phase_program":statistics.median(abs(num(r["logFC"])) for r in program) if program else None,"scenarios":{}}
    matched=[r for r in program if r["gene_id"] in d]
    out["phase_program_genes_in_disease_family"]=len(matched)
    for s in scenarios:
        col="mde_"+s
        vals=[num(r[col]) for r in d.values() if num(r[col]) is not None]
        med=statistics.median(vals) if vals else None
        above_med=sum(1 for r in program if med is not None and abs(num(r["logFC"]))>med)
        gm=[(abs(num(r["logFC"])),num(d[r["gene_id"]][col])) for r in matched]
        gm=[(a,b) for a,b in gm if b is not None]
        out["scenarios"][s]={"disease_median_mde_log2":med,"program_genes_above_disease_median_mde":above_med,
                             "share_above_disease_median_mde":above_med/len(program) if program else None,
                             "gene_matched_n":len(gm),"gene_matched_above_own_disease_mde":sum(1 for a,b in gm if a>b),
                             "gene_matched_share":(sum(1 for a,b in gm if a>b)/len(gm)) if gm else None}
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--work",default="work/mde_analysis");ap.add_argument("--git-head",required=True)
    ap.add_argument("--host-cleanup",default=None,help="optional JSON written by the host clean-up step")
    ap.add_argument("--freeze",default="receipts/mde_analysis_freeze.json");a=ap.parse_args()
    root=a.root.resolve();cfg=read(root/"config/mde_analysis.json");work=root/a.work
    pub=root/cfg["outputs"]["publish_dir"];rel=lambda p:str(pathlib.Path(p).resolve().relative_to(root))
    freeze=read(root/a.freeze)
    if freeze["phase"]!="MDE_ANALYSIS_RULES_FROZEN_PRE_COMPUTATION":raise ValueError("wrong freeze")
    freeze_sha=sha(root/a.freeze)
    for p,digest in freeze["artifacts"].items():
        if sha(root/p)!=digest:raise ValueError("frozen artifact changed on the published checkout: "+p)
    if pub.exists() and any(pub.iterdir()):raise ValueError("publish directory is not empty: "+str(pub))
    pub.mkdir(parents=True,exist_ok=True)
    artifacts={p:d for p,d in freeze["artifacts"].items()};published={};fits={};curves={}
    for name in FITS:
        base=work/"fits"/name;prov=read(base/"provenance.json")
        if prov["synthetic"] or prov["fit_name"]!=name or prov["freeze"]["sha256"]!=freeze_sha or prov["analysis"]!="MDE_PART_A_FROM_SAVED_FIT":raise ValueError("Part A provenance invalid: "+name)
        if prov["model"]["sha256"]!=(cfg["fits"][name].get("model_sha256") or prov["model"]["chain"]["sha256"]):raise ValueError("model hash not the pinned one: "+name)
        exempt=set(freeze.get("host_verification_exempt",[]))
        bad=[p for p,ok in prov["freeze"]["artifact_check"].items() if not ok and not p.startswith("docs/") and p not in exempt]
        if bad:raise ValueError("host freeze check failed for "+", ".join(bad))
        verify_outputs(prov,base)
        s=read(base/"summary.json");fits[name]={"summary":s,"provenance_sha256":sha(base/"provenance.json"),
            "host_doc_artifact_mismatches":[p for p,ok in prov["freeze"]["artifact_check"].items() if not ok]}
        gz_copy(base/"mde_per_gene.tsv",pub/f"{name}_mde_per_gene.tsv.gz");published[f"{name}_mde_per_gene.tsv.gz"]=rel(pub/f"{name}_mde_per_gene.tsv.gz")
        for fn in ("mde_by_expression_bin.tsv","summary.json","provenance.json"):
            shutil.copyfile(base/fn,pub/f"{name}_{fn}");published[f"{name}_{fn}"]=rel(pub/f"{name}_{fn}")
        if (base/"mde_ecdf.png").exists():shutil.copyfile(base/"mde_ecdf.png",pub/f"{name}_mde_ecdf.png");published[f"{name}_mde_ecdf.png"]=rel(pub/f"{name}_mde_ecdf.png")
        for fn in sorted(base.iterdir()):artifacts[rel(fn)]=sha(fn)
    for cohort in COHORTS:
        base=work/"curve"/cohort;prov=read(base/"provenance.json")
        if prov["synthetic"] or prov["cohort"]!=cohort or prov["freeze"]["sha256"]!=freeze_sha or prov["analysis"]!="MDE_PART_B_SUBSAMPLING_CURVE":raise ValueError("Part B provenance invalid: "+cohort)
        if not prov["identity_check"]["passed"]:raise ValueError("identity check did not pass: "+cohort)
        if read(base/"process.json")["status"]!="COMPLETED":raise ValueError("curve run did not complete: "+cohort)
        verify_outputs(prov,base)
        s=read(base/"summary.json");summ=rows(base/"curve_summary.tsv")
        curves[cohort]={"summary":s,"levels":[{k:(num(r[k]) if k not in ("design","invalid_reasons","cohort") else r[k]) for k in r} for r in summ],"provenance_sha256":sha(base/"provenance.json")}
        gz_copy(base/"curve_draws.tsv",pub/f"curve_{cohort}_draws.tsv.gz");published[f"curve_{cohort}_draws.tsv.gz"]=rel(pub/f"curve_{cohort}_draws.tsv.gz")
        gz_copy(base/"draw_ledger.tsv",pub/f"curve_{cohort}_draw_ledger.tsv.gz");published[f"curve_{cohort}_draw_ledger.tsv.gz"]=rel(pub/f"curve_{cohort}_draw_ledger.tsv.gz")
        for fn in ("curve_summary.tsv","summary.json","identity_check.json","draw_freeze.json","provenance.json"):
            shutil.copyfile(base/fn,pub/f"curve_{cohort}_{fn}");published[f"curve_{cohort}_{fn}"]=rel(pub/f"curve_{cohort}_{fn}")
        if (base/"curve.png").exists():shutil.copyfile(base/"curve.png",pub/f"curve_{cohort}.png");published[f"curve_{cohort}.png"]=rel(pub/f"curve_{cohort}.png")
        for fn in sorted(base.iterdir()):
            if fn.is_file():artifacts[rel(fn)]=sha(fn)
    bridges={c:bridge(rows(work/"fits"/f"v1_{c}"/"mde_per_gene.tsv"),rows(work/"fits"/f"phase_{c}"/"mde_per_gene.tsv"),cfg["bridge"]["disease_mde_reference_scenarios"]+["bh_k1","bh_k10","bh_k1000"]) for c in COHORTS}
    for v in published.values():artifacts[v]=sha(root/v)
    host=read(root/a.host_cleanup) if a.host_cleanup else None
    receipt={"schema_version":1,"phase":"MDE_ANALYSIS_COMPLETED","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
             "git_head":a.git_head,"protocol":cfg["protocol"],"freeze":{"path":a.freeze,"sha256":freeze_sha,"git_head":freeze["git_head"]},
             "artifacts":artifacts,"published_copies":published,"part_a":fits,"part_b":curves,"bridge":bridges,"host_cleanup":host,
             "attestations":{"computed_once":True,"new_disease_fits":False,"v1_artifacts_modified":False,"pass_fail_endpoints":False},
             "interpretation":"design sensitivity of the v1 analysis under realized variance estimates and a descriptive degradation curve of a known program; not evidence about any gene, not a sample-size decision"}
    out=root/cfg["outputs"]["results_receipt"];out.write_text(json.dumps(receipt,indent=2)+"\n")
    brief={n:{"family":f["summary"]["family_size"],"median_mde":{s:f["summary"]["scenarios"][s]["quantiles_log2"]["q0.5"] for s in f["summary"]["scenarios"]},
              "ci_within_1_per_gene":f["summary"]["ci_containment"]["per_gene"]["within_1"],"ci_within_1_bonferroni":f["summary"]["ci_containment"]["family_bonferroni"]["within_1"]} for n,f in fits.items()}
    print(json.dumps({"part_a":brief,"part_b":{c:curves[c]["summary"]["valid_draws"] for c in COHORTS},"bridge":{c:{s:bridges[c]["scenarios"][s]["share_above_disease_median_mde"] for s in ("alpha_0.05","bh_k100")} for c in COHORTS}},indent=2))
if __name__=="__main__":main()
