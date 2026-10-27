#!/usr/bin/env python3
"""Freeze the rules of the minimum detectable effect analysis before any MDE value is
computed on real fits: hash the protocol, the configuration, the code, the committed
receipts that pin the input fits, and the positive-control inputs; cross-check every pinned
hash in config/mde_analysis.json against the committed receipts; write the freeze receipt.
Standard library only. Paths in the receipt are relative to the project root."""
import argparse,datetime,hashlib,json,pathlib

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())

CODE=["scripts/mde_lib.R","scripts/mde_from_fits.R","scripts/mde_subsample_curve.R","scripts/mde_synthetic_fixture.R",
      "scripts/mde_freeze.py","scripts/mde_finalize.py","scripts/donor_de_lib.R",
      "requirements/statistical-runtime.explicit.txt"]
# Published receipts whose bytes differ from the analysis host's originals only by the
# documented placeholder substitution; the host run uses them through the hex values
# cross-checked below, so their hashes are recorded but verified on the published checkout only.
HOST_EXEMPT=["receipts/discovery_primary_execution.json","receipts/discovery_rank_freeze.json","receipts/heca_primary_execution.json",
             "receipts/phase_positive_control_f1_discovery_freeze.json","receipts/phase_positive_control_evaluation.json"]
RECEIPTS=["receipts/discovery_primary_execution.json","receipts/discovery_rank_freeze.json","receipts/heca_primary_execution.json",
          "receipts/phase_positive_control_f1_discovery_freeze.json","receipts/phase_positive_control_evaluation.json",
          "reports/phase_positive_control_f1/de_provenance.json","reports/phase_positive_control_f2/heca_de_provenance.json"]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--git-head",required=True,help="commit of the published checkout holding this protocol and code")
    ap.add_argument("--supersedes",default=None,help="commit of a previous freeze receipt that this one replaces, with --supersedes-reason")
    ap.add_argument("--supersedes-reason",default=None)
    ap.add_argument("--receipt",default="receipts/mde_analysis_freeze.json");a=ap.parse_args()
    root=a.root.resolve();cfg=read(root/"config/mde_analysis.json")
    if cfg["protocol"]!="docs/MDE_ANALYSIS_PROTOCOL.md":raise ValueError("protocol path differs")
    artifacts={}
    for p in [cfg["protocol"],"config/mde_analysis.json",cfg["estimator_config"],cfg["phase_config"]]+CODE+RECEIPTS:
        if not (root/p).exists():raise FileNotFoundError(p)
        artifacts[p]=sha(root/p)
    # pinned fit chains versus the committed receipts
    checks={}
    d=read(root/"receipts/discovery_primary_execution.json")
    checks["v1_discovery_provenance_sha256"]=cfg["fits"]["v1_discovery"]["provenance_sha256"]==d["DE_provenance_sha256"]
    r=read(root/"receipts/discovery_rank_freeze.json")
    checks["v1_discovery_model_sha256"]=cfg["fits"]["v1_discovery"]["model_sha256"]==r["artifacts"].get(cfg["fits"]["v1_discovery"]["model"])
    h=read(root/"receipts/heca_primary_execution.json")
    checks["v1_heca_model_sha256"]=cfg["fits"]["v1_heca"]["model_sha256"]==h["output_sha256"].get(cfg["fits"]["v1_heca"]["model"])
    f1=read(root/"reports/phase_positive_control_f1/de_provenance.json")
    checks["phase_discovery_model_sha256"]=cfg["fits"]["phase_discovery"]["model_sha256"]==f1["output_sha256"].get(cfg["fits"]["phase_discovery"]["model"])
    f2=read(root/"reports/phase_positive_control_f2/heca_de_provenance.json")
    checks["phase_heca_model_sha256"]=cfg["fits"]["phase_heca"]["model_sha256"]==f2["output_sha256"].get(cfg["fits"]["phase_heca"]["model"])
    # positive-control inputs: local bytes must equal the hashes the earlier freezes recorded
    sources={"discovery":read(root/"receipts/phase_positive_control_f1_discovery_freeze.json")["artifacts"],
             "heca":read(root/"receipts/phase_positive_control_evaluation.json")["artifacts"]}
    for cohort in ("discovery","heca"):
        for key in ("counts","clinical"):
            p=cfg["subsampling"][key][cohort]
            if p not in sources[cohort]:raise ValueError("earlier receipt does not pin "+p)
            digest=sha(root/p);artifacts[p]=digest
            checks[f"{cohort}_{key}_matches_earlier_receipt"]=digest==sources[cohort][p]
    if not all(checks.values()):raise ValueError("pinned hash cross-check failed: "+json.dumps(checks))
    for p in (cfg["outputs"]["results_receipt"],cfg["outputs"]["publish_dir"]):
        if (root/p).exists():raise ValueError("results already exist; the freeze must precede computation: "+p)
    if (a.supersedes is None)!=(a.supersedes_reason is None):raise ValueError("--supersedes and --supersedes-reason go together")
    receipt={"schema_version":1,"phase":"MDE_ANALYSIS_RULES_FROZEN_PRE_COMPUTATION",
             "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"git_head":a.git_head,
             "protocol":cfg["protocol"],"config":"config/mde_analysis.json","artifacts":artifacts,
             "host_verification_exempt":HOST_EXEMPT,
             "host_verification_note":"documentation and the listed published receipts are hashed here and verified on the published checkout; the analysis host keeps their pre-substitution bytes and verifies every other artifact",
             "pinned_fit_chains":{k:{kk:vv for kk,vv in v.items() if kk!="cohort"} for k,v in cfg["fits"].items()},
             "pinned_hash_cross_checks":checks,
             "scenarios":{"power":cfg["power"],"nominal_alpha":cfg["nominal_alpha"],"bh_q":cfg["bh_q"],"bh_scenarios_k":cfg["bh_scenarios_k"],
                          "ci_containment_log2_thresholds":cfg["ci_containment_log2_thresholds"]},
             "subsampling":{k:cfg["subsampling"][k] for k in ("draws_per_level","seed","rng","replace","equal_arms","levels_per_arm","designs","invalid_draw_policy","refilter","renormalize","identity_check")},
             "supersedes":None if a.supersedes is None else {"freeze_commit":a.supersedes,"reason":a.supersedes_reason},
             "attestations":{"mde_values_computed_on_real_fits_before_freeze":False,"code_tested_on_synthetic_fixture_only":True,
                             "new_disease_fits":False,"v1_artifacts_modified":False,"pass_fail_endpoints":False,"interpretation_rules_fixed":True},
             "interpretation":"bounds on what the v1 design could detect and a descriptive degradation curve of a known program; not evidence about any gene"}
    out=root/a.receipt;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"phase":receipt["phase"],"git_head":a.git_head,"artifacts":len(artifacts),"checks":checks},indent=2))
if __name__=="__main__":main()
