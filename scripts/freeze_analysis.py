#!/usr/bin/env python3
"""Create or verify immutable pre-effect/discovery checkpoints; all paths stay in project."""
import argparse, csv, datetime, hashlib, json, pathlib, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):h.update(block)
    return h.hexdigest()
def rel(path):
    return str(pathlib.Path(path).resolve().relative_to(ROOT))
def read(path): return json.loads(pathlib.Path(path).read_text())
def rows(path):
    with pathlib.Path(path).open() as f:return list(csv.DictReader(f,delimiter="\t"))
def verify_manifest(path, phase=None):
    obj=read(path)
    if phase is not None and obj.get("phase")!=phase:raise ValueError("Wrong freeze phase")
    if not obj.get("artifacts"):raise ValueError("Empty artifact map")
    for p,digest in obj["artifacts"].items():
        if sha(ROOT/p)!=digest:raise ValueError("Changed frozen artifact: "+p)
    return obj
def check_outputs(directory, receipt):
    for name, item in receipt["outputs"].items():
        if sha(pathlib.Path(directory)/name)!=item["sha256"]:
            raise ValueError("Changed output: "+name)
def validate_discovery(de_dir, rank_dir, null_dir):
    de_dir,rank_dir,null_dir=map(pathlib.Path,(de_dir,rank_dir,null_dir))
    dp,rp,np=read(de_dir/"provenance.json"),read(rank_dir/"provenance.json"),read(null_dir/"status.json")
    if dp.get("data_kind")!="public_cohort" or dp.get("cohort")!="discovery":
        raise ValueError("Only real discovery cohort effects can be frozen")
    if rp.get("status")!="DISCOVERY_RANKING_ARTIFACT_NOT_COMPLETE_EXTERNAL_FREEZE" or rp.get("external_data_used") is not False:
        raise ValueError("Invalid discovery ranking provenance")
    if np.get("status")!="DISCOVERY_NULL_MEMBERSHIPS_FROZEN_NOT_COMPLETE_EXTERNAL_FREEZE" or np.get("external_data_used") is not False:
        raise ValueError("Invalid null provenance")
    check_outputs(rank_dir,rp);check_outputs(null_dir,np)
    for original,digest in dp["output_sha256"].items():
        if sha(de_dir/pathlib.Path(original).name)!=digest:raise ValueError("Changed DE output")
    if rp["inputs"]["discovery_de"]["sha256"]!=sha(de_dir/"de_fitted.tsv"):
        raise ValueError("Ranking used different fitted DE")
    if rp["inputs"]["discovery_de_all"]["sha256"]!=sha(de_dir/"de_all.tsv"):
        raise ValueError("Ranking used different testing family")
    for key in ["universe.tsv","induced_edges.tsv","provenance.json",
                "B0.tsv","B0_top20.tsv","B1.tsv","B1_top20.tsv","G1.tsv","G1_top20.tsv"]:
        if np["inputs"][key]["sha256"]!=sha(rank_dir/key):raise ValueError("Nulls used different ranking: "+key)
    if np["inputs"]["discovery_de"]["sha256"]!=sha(de_dir/"de_fitted.tsv"):
        raise ValueError("Nulls used different discovery DE")
    universe=[r["gene_id"] for r in rows(rank_dir/"universe.tsv")]
    if len(universe)<50 or len(set(universe))!=len(universe):raise ValueError("Invalid comparative universe")
    for method in ["B0","B1","G1"]:
        full=rows(rank_dir/f"{method}.tsv")
        if len(full)!=len(universe) or {r["gene_id"] for r in full}!=set(universe):
            raise ValueError("Unequal method universe")
        for k in [20,50]:
            if rows(rank_dir/f"{method}_top{k}.tsv")!=full[:k]:
                raise ValueError("Top list is not the frozen full-rank prefix")
        state=np["methods"][method]
        members=rows(null_dir/f"{method}_null_membership.tsv")
        if state["status"]=="ESTIMABLE":
            if state["sampled_sets"]!=200 or len(members)!=4000:raise ValueError("Incomplete matched sets")
        elif state["status"]=="NOT_ESTIMABLE":
            if state["sampled_sets"]!=0 or members or not state["reason"]:raise ValueError("Missing explicit null failure")
        else:raise ValueError("Unknown matched set status")
    if {r["gene_id"] for r in rows(null_dir/"all_U_directions.tsv")}!=set(universe):
        raise ValueError("Incomplete all-U prediction directions")
    return {"universe_size":len(universe),
            "null_status":{m:np["methods"][m]["status"] for m in ["B0","B1","G1"]}}
def base_paths():
    fixed=["docs/RESEARCH_PROTOCOL.md","docs/PRE_EFFECT_METHOD_DECISIONS.md",
           "docs/STATISTICAL_IMPLEMENTATION_DECISIONS.md","config/statistical_analysis.json",
           "config/ranking_algorithm.json","requirements/statistical-runtime.explicit.txt",
           "requirements/statistical-runtime.yml",
           "receipts/runtime_statistics_install.json","receipts/runtime_statistics_R_smoke.json",
           "receipts/m0_input_gate_decision.json","main.nf","nextflow.config",
           "data/reference/string/direct_evidence_gene_edges.tsv.gz"]
    fixed += [str(p.relative_to(ROOT)) for p in sorted((ROOT/"scripts").glob("*")) if p.is_file() and p.suffix in (".R",".py",".sh")]
    for cohort in ["discovery","heca"]:
        fixed += [f"data/processed/{cohort}/counts_stable.tsv.gz",
                  f"manifests/final_{cohort}_clinical_selected.tsv",
                  f"reports/final_clinical_design_{cohort}_primary.tsv"]
    return sorted(set(fixed))
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase",choices=["implementation","discovery"])
    ap.add_argument("--output",type=pathlib.Path)
    ap.add_argument("--verify",type=pathlib.Path)
    ap.add_argument("--pre-effect-freeze",type=pathlib.Path)
    ap.add_argument("--de-dir",type=pathlib.Path)
    ap.add_argument("--rank-dir",type=pathlib.Path)
    ap.add_argument("--null-dir",type=pathlib.Path)
    a=ap.parse_args()
    if a.verify:
        obj=verify_manifest(a.verify);print(json.dumps({"verified":True,"phase":obj["phase"],"artifacts":len(obj["artifacts"])}));return
    if not a.output or not a.phase:ap.error("--phase and --output required")
    if a.output.exists():raise ValueError("Existing freeze preserved; verify it or use a new dated path")
    if subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip():
        raise ValueError("Commit reviewed code/receipts before freezing")
    m0=read(ROOT/"receipts/m0_input_gate_decision.json")
    if m0.get("input_feasibility_gate_passed") is not True:raise ValueError("M0 input gate not passed")
    for p,h in m0["evidence_sha256"].items():
        if sha(ROOT/p)!=h:raise ValueError("M0 evidence changed: "+p)
    paths=base_paths();products={};validation={}
    if a.phase=="discovery":
        if not all([a.pre_effect_freeze,a.de_dir,a.rank_dir,a.null_dir]):ap.error("Discovery output directories and pre-effect freeze required")
        verify_manifest(a.pre_effect_freeze,"STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT")
        validation=validate_discovery(a.de_dir,a.rank_dir,a.null_dir)
        products={"discovery_de_all":rel(a.de_dir/"de_all.tsv"),
                  "discovery_de_fitted":rel(a.de_dir/"de_fitted.tsv"),
                  "universe":rel(a.rank_dir/"universe.tsv"),"induced_edges":rel(a.rank_dir/"induced_edges.tsv"),
                  "all_U_directions":rel(a.null_dir/"all_U_directions.tsv"),
                  "null_bins":rel(a.null_dir/"universe_bins.tsv"),
                  "null_status":rel(a.null_dir/"status.json"),
                  "null_memberships":rel(a.null_dir/"null_memberships.tsv"),
                  "statistical_config":"config/statistical_analysis.json",
                  "runtime_lock":"requirements/statistical-runtime.explicit.txt",
                  "network":"data/reference/string/direct_evidence_gene_edges.tsv.gz",
                  "discovery_counts":"data/processed/discovery/counts_stable.tsv.gz",
                  "discovery_design":"reports/final_clinical_design_discovery_primary.tsv",
                  "discovery_clinical":"manifests/final_discovery_clinical_selected.tsv"}
        for method in ["B0","B1","G1"]:
            for suffix,file_suffix in [("full",""),("top20","_top20"),("top50","_top50")]:
                products[f"{method}_{suffix}"]=rel(a.rank_dir/f"{method}{file_suffix}.tsv")
        for directory in [a.de_dir,a.rank_dir,a.null_dir]:
            paths += [rel(p) for p in sorted(directory.iterdir()) if p.is_file()]
        paths += [rel(a.pre_effect_freeze)]
    obj={"schema_version":1,"phase":"STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT" if a.phase=="implementation" else "DISCOVERY_FROZEN_PRE_EXTERNAL",
         "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
         "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
         "artifacts":{p:sha(ROOT/p) for p in sorted(set(paths))},
         "products":products,"validation":validation,
         "external_disease_effects_computed":False,
         "claim":"auditable retrospective checkpoint, not preregistration or external scientific validation"}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("x") as f:json.dump(obj,f,indent=2);f.write("\n")
    print(json.dumps({"phase":obj["phase"],"artifacts":len(obj["artifacts"]),"output":str(a.output)}))
if __name__=="__main__":main()
