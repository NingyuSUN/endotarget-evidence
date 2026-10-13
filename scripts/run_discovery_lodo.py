#!/usr/bin/env python3
"""Leave each discovery donor out; refit, refilter and rebuild frozen-network ranks."""
import argparse
import csv
import datetime
import hashlib
import json
import os
import pathlib
import subprocess
import sys
import time

from rank_targets import digest, execute, table_rows

ROOT=pathlib.Path(__file__).resolve().parents[1]
METHODS=("B0","B1","G1")

def top20_metrics(reference,fold):
    if len(reference)!=20 or len(set(reference))!=20 or len(fold)!=20 or len(set(fold))!=20:
        raise ValueError("Both top20 lists must have 20 unique genes")
    overlap=len(set(reference)&set(fold))
    return {"intersection":overlap,"reference_n":20,"fold_n":20,"retention":overlap/20,"jaccard":overlap/len(set(reference)|set(fold))}

def atomic_json(path,value):
    temporary=path.with_name(path.name+f".partial.{os.getpid()}")
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)

def verified_freeze(path,phase):
    freeze=json.loads(path.read_text())
    if freeze.get("phase")!=phase or not freeze.get("artifacts"):
        raise ValueError("Wrong freeze phase/schema")
    for item,sha in freeze["artifacts"].items():
        if digest(item)!=sha: raise ValueError(f"Frozen artifact changed: {item}")
    return freeze

def checkpoint_binding(paths,synthetic):
    return {"synthetic":synthetic,"sha256":{str(p):digest(p) for p in paths},
            "python":sys.version,"executable":sys.executable,
            "policy":"22_folds_no_replacement_refilter_frozen_network"}

def run(args):
    os.chdir(ROOT)
    synthetic=args.synthetic_test
    if synthetic:
        for p in (args.counts,args.clinical,args.baseline_dir,args.edges,args.output_dir):
            if not p.resolve().is_relative_to((ROOT/"work/synthetic").resolve()):
                raise ValueError("Synthetic inputs/outputs must remain in work/synthetic")
    else:
        dfreeze=verified_freeze(args.discovery_freeze,"DISCOVERY_FROZEN_PRE_EXTERNAL")
        pre=verified_freeze(args.pre_effect_freeze,"STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT")
        expected={"counts":dfreeze["products"]["discovery_counts"],"clinical":dfreeze["products"]["discovery_clinical"],"edges":dfreeze["products"]["network"]}
        for key,p in expected.items():
            if getattr(args,key).resolve()!=(ROOT/p).resolve():raise ValueError("Production requires frozen discovery input paths")
        for method in METHODS:
            if (args.baseline_dir/f"{method}_top20.tsv").resolve()!=(ROOT/dfreeze["products"][f"{method}_top20"]).resolve():raise ValueError("Production baseline must be full-discovery frozen list")
        if args.config.resolve()!=(ROOT/"config/statistical_analysis.json").resolve() or args.ranking_config.resolve()!=(ROOT/"config/ranking_algorithm.json").resolve():raise ValueError("Production configs must be frozen")
    clinical=list(table_rows(args.clinical,["donor_id","disease","primary_lineage_cells"]))
    donors=[r["donor_id"] for r in clinical]
    if len(donors)!=len(set(donors)) or (not synthetic and len(donors)!=22):raise ValueError("Unique original discovery people required")
    baseline={m:[r["gene_id"] for r in table_rows(args.baseline_dir/f"{m}_top20.tsv",["gene_id"])] for m in METHODS}
    for m in METHODS:top20_metrics(baseline[m],baseline[m])
    paths=[args.counts,args.clinical,args.edges,args.config,args.ranking_config,pathlib.Path(__file__),ROOT/"scripts/resample_donor_models.R",ROOT/"scripts/donor_de_lib.R",ROOT/"scripts/rank_targets.py",ROOT/"scripts/statistical_runtime_env.sh",ROOT/"requirements/statistical-runtime.explicit.txt"]
    paths.extend(args.baseline_dir/f"{m}_top20.tsv" for m in METHODS)
    if not synthetic: paths.extend((args.discovery_freeze,args.pre_effect_freeze))
    binding=checkpoint_binding(paths,synthetic)
    out=args.output_dir
    out.mkdir(parents=True,exist_ok=True)
    bound=out/"binding.json"
    if bound.exists():
        if json.loads(bound.read_text())!=binding:raise ValueError("Unsafe LODO resume: code/input/runtime changed")
    else:
        if any(out.iterdir()):raise ValueError("Nonempty unbound LODO outputs preserved")
        atomic_json(bound,binding)
    fingerprint=digest(bound)
    process={"pid":os.getpid(),"host":os.uname().nodename,"started_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"status":"RUNNING","binding_sha256":fingerprint,"exit_status":None,"workers":1,"BLAS_threads":1}
    atomic_json(out/"process.json",process)
    rows=[]
    try:
        for index,donor in enumerate(donors,1):
            fold=out/f"fold_{index:02d}"
            receipt=fold/"checkpoint.json"
            if receipt.exists():
                result=json.loads(receipt.read_text())
                if result.get("binding_sha256")!=fingerprint or result.get("dropped_donor")!=donor:raise ValueError("Unsafe fold checkpoint")
                for p,sha in result.get("output_sha256",{}).items():
                    if digest(fold/p)!=sha:raise ValueError(f"Fold checkpoint output changed: {p}")
            else:
                if fold.exists():raise ValueError(f"Incomplete fold preserved; inspect and move partial {fold} before retry")
                fold.mkdir()
                remaining=[r for r in clinical if r["donor_id"]!=donor]
                ncase=sum(r["disease"]=="ENDO" for r in remaining);ncontrol=sum(r["disease"]=="non-ENDO" for r in remaining)
                result={"dropped_donor":donor,"binding_sha256":fingerprint,"ENDO":ncase,"non_ENDO":ncontrol,"valid":False,"failed_fold_replaced":False}
                started=time.monotonic()
                try:
                    if min(ncase,ncontrol)<8 or any(float(r["primary_lineage_cells"])<50 for r in remaining):raise ValueError("Original group/cell eligibility minimum failed")
                    command=["Rscript","scripts/resample_donor_models.R","--mode","lodo-fold","--counts",str(args.counts),"--clinical",str(args.clinical),"--config",str(args.config),"--drop-donor",donor,"--output-dir",str(fold/"de"),"--pre-effect-freeze",str(args.pre_effect_freeze),"--discovery-freeze",str(args.discovery_freeze)]
                    if synthetic:command.extend(("--synthetic-test","true"))
                    env=os.environ.copy();env["OPENBLAS_NUM_THREADS"]="1";env["OMP_NUM_THREADS"]="1"
                    with (fold/"fit.log").open("w") as log:
                        child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env)
                        result["fit_pid"]=child.pid
                        process["active_fold"]=index;process["child_pid"]=child.pid
                        atomic_json(out/"process.json",process)
                        code=child.wait();result["fit_exit_status"]=code
                    if code:raise ValueError(f"Frozen donor fit failed with exit {code}; see fit.log")
                    rank_receipt=execute(fold/"de/de_fitted.tsv",args.edges,fold/"ranking",args.ranking_config,de_all_path=fold/"de/de_all.tsv",synthetic_test=synthetic)
                    result["metrics"]={m:top20_metrics(baseline[m],[r["gene_id"] for r in table_rows(fold/"ranking"/f"{m}_top20.tsv",["gene_id"])]) for m in METHODS}
                    result["universe_genes"]=rank_receipt["counts"]["final_universe_genes"]
                    result["valid"]=True;result["reason"]=""
                except (ValueError,OSError) as exc: result["reason"]=str(exc)
                result["elapsed_seconds"]=time.monotonic()-started
                result["output_sha256"]={str(p.relative_to(fold)):digest(p) for p in fold.rglob("*") if p.is_file()}
                atomic_json(receipt,result)
            for m in METHODS:
                metric=result.get("metrics",{}).get(m,{})
                rows.append({"fold":index,"dropped_donor":donor,"method":m,"valid":result["valid"],"reason":result.get("reason",""),"ENDO":result["ENDO"],"non_ENDO":result["non_ENDO"],"intersection":metric.get("intersection","NA"),"retention":metric.get("retention","NA"),"jaccard":metric.get("jaccard","NA")})
            print("LODO_CHECKPOINT",f"{index}/{len(donors)}",donor,"valid=",result["valid"],flush=True)
        with (out/"fold_metrics.tsv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");writer.writeheader();writer.writerows(rows)
        if checkpoint_binding(paths,synthetic)!=binding:raise ValueError("Input/code/runtime changed during LODO")
        valid=sum(r["valid"] for r in rows[::3])
        summary={"total_folds":len(donors),"valid_folds":valid,"failed_folds":len(donors)-valid,"failed_folds_replaced":False,"binding_sha256":fingerprint,"synthetic":synthetic,"scientific_validation":False,"metrics_meaning":"fixed20_retention_and_intersection_over_union_discovery_stability"}
        atomic_json(out/"summary.json",summary)
        process["status"]="COMPLETED";process["exit_status"]=0
        return summary
    except BaseException:
        process["status"]="INTERRUPTED_OR_FAILED";process["exit_status"]=1
        raise
    finally:
        process["finished_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
        atomic_json(out/"process.json",process)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts",type=pathlib.Path,default=pathlib.Path("data/processed/discovery/counts_stable.tsv.gz"))
    p.add_argument("--clinical",type=pathlib.Path,default=pathlib.Path("manifests/final_discovery_clinical_selected.tsv"))
    p.add_argument("--edges",type=pathlib.Path,default=pathlib.Path("data/reference/string/direct_evidence_gene_edges.tsv.gz"))
    p.add_argument("--baseline-dir",type=pathlib.Path,default=pathlib.Path("work/statistical/discovery_primary/ranking"))
    p.add_argument("--config",type=pathlib.Path,default=pathlib.Path("config/statistical_analysis.json"))
    p.add_argument("--ranking-config",type=pathlib.Path,default=pathlib.Path("config/ranking_algorithm.json"))
    p.add_argument("--pre-effect-freeze",type=pathlib.Path,default=pathlib.Path("receipts/statistical_pre_effect_freeze.json"))
    p.add_argument("--discovery-freeze",type=pathlib.Path,default=pathlib.Path("receipts/discovery_rank_freeze.json"))
    p.add_argument("--output-dir",type=pathlib.Path,required=True)
    p.add_argument("--synthetic-test",action="store_true")
    args=p.parse_args()
    try:print(json.dumps(run(args),indent=2))
    except (ValueError,OSError) as exc:p.exit(1,f"LODO failed: {exc}\n")
if __name__=="__main__":main()
