#!/usr/bin/env python3
"""Freeze the rules of the secretory-phase program-level disease contrast before any donor is
scored or any disease-within-secretory quantity is computed: hash the protocol, the configuration,
the code, the committed receipts that pin the inputs, the positive-control count matrices and
ledgers and the published phase-fit copies; cross-check every pinned input hash in
config/secretory_program.json against the committed receipts; write the freeze receipt.
Standard library only. Paths in the receipt are relative to the project root."""
import argparse,datetime,hashlib,json,pathlib

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())

CODE=["scripts/secretory_program_lib.R","scripts/secretory_program_analysis.R","scripts/secretory_program_synthetic_fixture.R",
      "scripts/secretory_program_freeze.py","scripts/secretory_program_finalize.py","scripts/mde_lib.R","scripts/donor_de_lib.R",
      "requirements/statistical-runtime.explicit.txt"]
RECEIPTS=["receipts/phase_positive_control_f1_discovery_freeze.json","receipts/phase_positive_control_evaluation.json"]
# Verified on the published checkout only: the receipts keep pre-substitution bytes on the analysis
# host and the published phase-fit copies exist only in the published checkout (the host verifies
# the work-tree originals against the pinned hashes).
HOST_EXEMPT=RECEIPTS+["reports/phase_positive_control_f1/de_all.tsv.gz","reports/phase_positive_control_f2/heca_de_all.tsv.gz"]
# receipt key that pins each input, per cohort (HECA normalization and design through published byte copies)
PINS={"discovery":{k:(RECEIPTS[0],None) for k in ("counts","clinical","phase_fit_table","expression_filter","normalization","design","published_phase_fit_table")},
      "heca":{"counts":(RECEIPTS[1],None),"clinical":(RECEIPTS[1],None),"phase_fit_table":(RECEIPTS[1],None),"expression_filter":(RECEIPTS[1],None),
              "normalization":(RECEIPTS[1],"reports/phase_positive_control_f2/heca_de_normalization.tsv"),
              "design":(RECEIPTS[1],"reports/phase_positive_control_f2/heca_de_design.tsv"),"published_phase_fit_table":(RECEIPTS[1],None)}}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--git-head",required=True,help="commit of the published checkout holding this protocol and code")
    ap.add_argument("--supersedes",default=None);ap.add_argument("--supersedes-reason",default=None)
    ap.add_argument("--receipt",default="receipts/secretory_program_freeze.json");a=ap.parse_args()
    root=a.root.resolve();cfg=read(root/"config/secretory_program.json")
    if cfg["protocol"]!="docs/SECRETORY_PROGRAM_PROTOCOL.md":raise ValueError("protocol path differs")
    artifacts={}
    for p in [cfg["protocol"],"config/secretory_program.json",cfg["estimator_config"],cfg["phase_config"]]+CODE+RECEIPTS:
        if not (root/p).exists():raise FileNotFoundError(p)
        artifacts[p]=sha(root/p)
    checks={};pinned={}
    for cohort,pins in PINS.items():
        inp=cfg["inputs"][cohort];pinned[cohort]={}
        for key,(receipt,via) in pins.items():
            path=inp[key];digest=inp[key+"_sha256"];rec=read(root/receipt)["artifacts"]
            lookup=via or path
            if lookup not in rec:raise ValueError(f"{receipt} does not pin {lookup}")
            checks[f"{cohort}_{key}_pinned_by_receipt"]=rec[lookup]==digest
            pinned[cohort][key]={"path":path,"sha256":digest,"receipt":receipt,"receipt_key":lookup}
            if path.startswith(("data/","reports/")):   # present in the published checkout: bytes must match the pin
                if not (root/path).exists():raise FileNotFoundError(path)
                artifacts[path]=sha(root/path);checks[f"{cohort}_{key}_local_bytes_match"]=artifacts[path]==digest
    if not all(checks.values()):raise ValueError("pinned hash cross-check failed: "+json.dumps(checks,indent=1))
    for p in (cfg["outputs"]["results_receipt"],cfg["outputs"]["publish_dir"]):
        if (root/p).exists():raise ValueError("results already exist; the freeze must precede the analysis: "+p)
    if (a.supersedes is None)!=(a.supersedes_reason is None):raise ValueError("--supersedes and --supersedes-reason go together")
    receipt={"schema_version":1,"phase":"SECRETORY_PROGRAM_RULES_FROZEN_PRE_ANALYSIS",
             "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"git_head":a.git_head,
             "protocol":cfg["protocol"],"config":"config/secretory_program.json","artifacts":artifacts,
             "host_verification_exempt":HOST_EXEMPT,
             "host_verification_note":"documentation, the listed receipts and the published phase-fit copies are hashed here and verified on the published checkout; the analysis host verifies every other artifact and the pinned work-tree inputs below",
             "pinned_inputs":pinned,"pinned_hash_cross_checks":checks,
             "signature":cfg["signature"],"score":cfg["score"],"models":cfg["models"],"validity_gate":cfg["validity_gate"],
             "gene_level":cfg["gene_level"],"evidence_labels":cfg["evidence_labels"],
             "supersedes":None if a.supersedes is None else {"freeze_commit":a.supersedes,"reason":a.supersedes_reason},
             "attestations":{"donor_scores_computed_before_freeze":False,"disease_within_secretory_quantity_computed_before_freeze":False,
                             "phase_fit_disease_coefficient_examined":False,"code_tested_on_synthetic_fixture_only":True,
                             "new_disease_fits":True,"v1_artifacts_modified":False,"positive_control_artifacts_modified":False,
                             "evidence_labels_fixed":True,"interpretation_rules_fixed":True},
             "interpretation":"a program-level association test in eutopic endometrium with two small secretory arms; labels are fixed in advance; no gene-level claim and no target nomination follows from it"}
    out=root/a.receipt;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"phase":receipt["phase"],"git_head":a.git_head,"artifacts":len(artifacts),"checks_passed":sum(checks.values()),"checks":len(checks)},indent=2))
if __name__=="__main__":main()
