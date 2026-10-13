#!/usr/bin/env python3
"""Verify completed synthetic replication, cached resume and pre-execution rejection gates."""
import argparse,csv,datetime,hashlib,json,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def readtrace(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fixture",type=pathlib.Path,required=True)
    ap.add_argument("--first-trace",type=pathlib.Path,required=True)
    a=ap.parse_args();folder=a.fixture.resolve()
    if not folder.is_relative_to(ROOT/"work/synthetic"):raise ValueError("synthetic fixture required")
    initial=readtrace(a.first_trace)
    if len(initial)!=2 or any(x["status"]!="COMPLETED" or x["exit"]!="0" for x in initial):
        raise ValueError("first run did not complete both stages")
    first=json.loads((folder/"results_first/evaluation/summary.json").read_text())
    params=json.loads((folder/"params.json").read_text())
    if params["synthetic_test"] is not True:raise ValueError("synthetic mode required")
    checks=[];now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    def execute(label,modifications,outdir,expected=None,resume=False):
        selected={**params,**modifications};param=folder/(label+"_params.json")
        param.write_text(json.dumps(selected,indent=2)+"\n")
        trace=folder/(label+"_trace.tsv")
        cmd=["nextflow","-C","replication.config","-log",str(folder/(label+"_nextflow.log")),
             "run","replication.nf","-profile","srv01","-work-dir",str(folder/"nfwork"),
             "-params-file",str(param),"--outdir",str(outdir),"-with-trace",str(trace)]
        if resume:cmd+=["-resume"]
        p=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
        (folder/(label+"_console.log")).write_text(p.stdout)
        if expected is None:
            if p.returncode:raise ValueError(label+" failed: "+p.stdout[-2000:])
            rows=readtrace(trace)
            if len(rows)!=2 or any(r["status"]!="CACHED" for r in rows):raise ValueError("resume was not fully cached")
        elif p.returncode==0 or expected not in p.stdout:
            raise ValueError("negative gate did not fail as expected: "+label+" "+p.stdout[-1000:])
        checks.append({"name":label,"passed":True,"exit_code":p.returncode,"command":cmd,
                       "expected_error":expected,"console":str((folder/(label+"_console.log")).relative_to(ROOT))})
        print(label+" PASS",flush=True)
    execute("resume",{},folder/"results_resume",resume=True)
    for name in ("selected_genes.tsv","bootstrap_interface.tsv","null_sets.tsv","dossier_pool.tsv","dossier_selection.tsv"):
        if sha(folder/"results_first/evaluation"/name)!=sha(folder/"results_resume/evaluation"/name):
            raise ValueError("cached output changed: "+name)
    execute("missing_freeze",{"evaluation_freeze":None},folder/"rejected_missing",
            "--evaluation_freeze and --outdir required")
    altered=json.loads((ROOT/params["evaluation_freeze"]).read_text())
    altered["artifacts"]["scripts/evaluate_replication.py"]="0"*64
    bad=folder/"corrupted_evaluation_freeze.json";bad.write_text(json.dumps(altered)+"\n")
    execute("corrupt_hash",{"evaluation_freeze":str(bad.relative_to(ROOT))},folder/"rejected_corrupt","Frozen artifact changed:")
    execute("wrong_phase",{"synthetic_test":False},folder/"rejected_phase","Replication freeze phase/schema mismatch")
    execute("existing_output",{},folder/"results_first","Output must be new")
    receipt={"status":"PASS","data_kind":"synthetic_test_only","scientific_validation":False,
             "started_utc":now,"completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
             "first_trace":str(a.first_trace),"first_stages":initial,"checks":checks,
             "fixture":str(folder.relative_to(ROOT)),
             "code_sha256":{n:sha(ROOT/n) for n in ["replication.nf","replication.config",
              "scripts/evaluate_replication.py","tests/make_replication_synthetic.py",
              "tests/check_replication_nextflow.py"]},
             "scope":"actual native Nextflow execution; no real external effects or container execution"}
    p=ROOT/"receipts/replication_nextflow_tests.json"
    with p.open("x") as f:json.dump(receipt,f,indent=2);f.write("\n")
    print(str(p.relative_to(ROOT)))
if __name__=="__main__":main()
