#!/usr/bin/env python3
"""Freeze tested replication implementation before first external disease effect."""
import argparse,datetime,hashlib,json,pathlib,subprocess
from freeze_analysis import verify_manifest,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
NEW=["replication.nf","replication.config","scripts/evaluate_replication.py",
     "scripts/resample_donor_models.R","scripts/run_discovery_lodo.py","scripts/freeze_replication.py",
     "tests/test_replication_evaluation.py","tests/test_resampling.R","tests/test_lodo.py",
     "tests/make_replication_synthetic.py","tests/check_replication_nextflow.py","docs/REPLICATION_EVALUATION.md","docs/RESAMPLING.md",
     "docs/REPLICATION_EXECUTION.md","receipts/replication_evaluation_tests.json",
     "receipts/resampling_tests.json","receipts/replication_nextflow_tests.json",
     "receipts/replication_independent_review.json","receipts/replication_fix_adjudication.json",
     "config/replication_run.json"]
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output",required=True,type=pathlib.Path)
    a=ap.parse_args();out=a.output.resolve()
    if not out.is_relative_to(ROOT/"receipts"):raise ValueError("freeze must be in receipts")
    if out.exists():raise ValueError("preserve existing freeze")
    if subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip():
        raise ValueError("commit intentional files before freezing")
    if (ROOT/"work/statistical/heca_primary").exists():raise ValueError("external run already exists")
    old=verify_manifest(ROOT/"receipts/discovery_rank_freeze.json","DISCOVERY_FROZEN_PRE_EXTERNAL")
    if old.get("external_disease_effects_computed") is not False:raise ValueError("discovery freeze invalid")
    verify_manifest(ROOT/"receipts/statistical_pre_effect_freeze.json","STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT")
    for receipt in ("replication_evaluation_tests.json","resampling_tests.json","replication_nextflow_tests.json"):
        if not (ROOT/"receipts"/receipt).is_file():raise ValueError("missing test evidence")
    names=sorted(set(old["artifacts"])|set(NEW)|{"receipts/discovery_rank_freeze.json","receipts/statistical_pre_effect_freeze.json"})
    obj={"schema_version":1,"phase":"REPLICATION_IMPLEMENTATION_FROZEN_PRE_EXTERNAL",
         "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
         "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
         "artifacts":{n:sha(ROOT/n) for n in names},
         "external_disease_effects_computed":False,
         "planned_execution":"Nextflow native runtime on SRV01; Docker execution not verified",
         "claim":"auditable pre-external checkpoint, not prospective preregistration"}
    with out.open("x") as f:json.dump(obj,f,indent=2);f.write("\n")
    print(json.dumps({"output":str(out.relative_to(ROOT)),"artifacts":len(names),"phase":obj["phase"]}))
if __name__=="__main__":main()
