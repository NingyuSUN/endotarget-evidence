#!/usr/bin/env bash
# Full synthetic model -> ranking -> matched-null chain plus strict rejection/cache checks.
set -euo pipefail
cd /home/analyst/te
source scripts/statistical_runtime_env.sh
verify_only=false
if [[ "${1:-}" == "--verify-only" ]]; then
  verify_only=true
  fixture="$2"
else
  fixture="$(python tests/make_nextflow_synthetic.py)"
fi
export ENDOTARGET_SYNTHETIC_FIXTURE="$fixture"
export ENDOTARGET_NEXTFLOW_VERIFY_ONLY="$verify_only"
run_flow() {
  local output_name="$1" trace_name="$2"
  shift 2
  nextflow -log "$fixture/${output_name}.log" run main.nf -profile srv01 \
    -params-file "$fixture/params.json" --outdir "$fixture/$output_name" \
    -work-dir "$fixture/nfwork" -with-trace "$fixture/$trace_name.tsv" "$@" \
    > "$fixture/${output_name}.stdout" 2>&1
}
if [[ "$verify_only" == false ]]; then
run_flow results_first trace_first
run_flow results_resume trace_resume -resume
# Alter one synthetic file without updating freeze. No process may be submitted.
printf '\n' >> "$fixture/clinical.tsv"
if run_flow rejected_mutation trace_rejected -resume; then
  echo "FAIL: changed frozen input accepted" >&2
  exit 1
fi
# Restore exact original bytes; no production file is changed.
python - <<'PY'
import os,pathlib
p=pathlib.Path(os.environ["ENDOTARGET_SYNTHETIC_FIXTURE"])/"clinical.tsv"
content=p.read_bytes()
assert content.endswith(b"\n\n")
p.write_bytes(content[:-1])
PY
if nextflow -log "$fixture/missing_freeze.log" run main.nf -profile srv01 \
  --outdir "$fixture/rejected_missing_freeze" > "$fixture/missing_freeze.stdout" 2>&1; then
  echo "FAIL: missing freeze accepted" >&2
  exit 1
fi
# Changing the synthetic-test switch must reject the synthetic freeze before any task.
if nextflow -log "$fixture/wrong_phase.log" run main.nf -profile srv01 \
  -params-file "$fixture/params.json" --synthetic_test false --outdir "$fixture/rejected_wrong_phase" \
  > "$fixture/wrong_phase.stdout" 2>&1; then
  echo "FAIL: synthetic manifest authorized production" >&2
  exit 1
fi
fi
python - <<'PY'
import json,os,pathlib
folder=pathlib.Path(os.environ["ENDOTARGET_SYNTHETIC_FIXTURE"])
frozen=json.loads((folder/"freeze.json").read_text())
del frozen["artifacts"]["scripts/rank_targets.py"]
(folder/"omitted_code_freeze.json").write_text(json.dumps(frozen,indent=2)+"\n")
PY
if nextflow -log "$fixture/omitted_code.log" run main.nf -profile srv01 \
  -params-file "$fixture/params.json" --pre_effect_freeze "${fixture#/home/analyst/te/}/omitted_code_freeze.json" \
  --outdir "$fixture/rejected_omitted_code" > "$fixture/omitted_code.stdout" 2>&1; then
  echo "FAIL: incomplete code freeze accepted" >&2
  exit 1
fi
python - <<'PY'
import csv,datetime,hashlib,json,os,pathlib,platform,subprocess
root=pathlib.Path.cwd()
folder=pathlib.Path(os.environ["ENDOTARGET_SYNTHETIC_FIXTURE"])
read=lambda p:json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def trace(name):
    with (folder/(name+".tsv")).open() as f:return list(csv.DictReader(f,delimiter="\t"))
first,resume=trace("trace_first"),trace("trace_resume")
assert len(first)==len(resume)==3,(first,resume)
assert all(row["status"]=="COMPLETED" and row["exit"]=="0" for row in first),first
assert all(row["status"]=="CACHED" for row in resume),resume
assert "Frozen artifact changed:" in (folder/"rejected_mutation.stdout").read_text()
assert "--pre_effect_freeze is required" in (folder/"missing_freeze.stdout").read_text()
assert "Freeze phase/schema does not authorize" in (folder/"wrong_phase.stdout").read_text()
assert "Required input/code/runtime omitted from freeze: scripts/rank_targets.py" in (folder/"omitted_code.stdout").read_text()
de=read(folder/"results_first/de/provenance.json")
rank=read(folder/"results_first/ranking/provenance.json")
nulls=read(folder/"results_first/nulls/status.json")
assert de["data_kind"]=="synthetic_test_only",de
assert rank["status"]=="SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE",rank
assert de["tested_family_size"]==1500 and rank["counts"]["final_universe_genes"]==1500
for method in ("B0","B1","G1"):
    status=nulls["methods"][method]
    assert status["status"]=="ESTIMABLE" and status["sampled_sets"]==200 and status["genes_per_set"]==20,(method,status)
    with (folder/"results_first/nulls"/(method+"_null_membership.tsv")).open() as handle:
        members=list(csv.DictReader(handle,delimiter="\t"))
    assert len(members)==4000,(method,len(members))
for stage in ("de","ranking","nulls"):
    for path in (folder/"results_first"/stage).iterdir():
        assert path.is_file()
        assert sha(path)==sha(folder/"results_resume"/stage/path.name),(stage,path)
frozen=read(folder/"freeze.json")
stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
receipt={"status":"PASS","data_kind":"synthetic","utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
 "hostname":subprocess.check_output(["hostname"],text=True).strip(),
 "project":str(root),"fixture":str(folder),"fixture_sha256":sha(folder/"fixture.json"),
 "freeze_sha256":sha(folder/"freeze.json"),"freeze_phase":frozen["phase"],"frozen_artifacts":frozen["artifacts"],
 "fixture_summary":read(folder/"fixture.json"),"first_trace":first,"resume_trace":resume,
 "all_three_processes_completed":True,"all_three_processes_cached_on_resume":True,
 "changed_frozen_input_rejected_before_tasks":True,"missing_freeze_rejected":True,"synthetic_freeze_cannot_authorize_production":True,
 "omitted_required_code_rejected_before_tasks":True,"resume_outputs_byte_identical":True,"donor_model_audit":de["audit"],"ranking_counts":rank["counts"],
 "null_status":nulls,"scientific_validation":False,"real_cohort_effects_computed":False,
 "docker_execution_verified":False,"AWS_execution_verified":False,
 "command":"bash tests/smoke_nextflow.sh",
 "verification_replay_command":("bash tests/smoke_nextflow.sh --verify-only "+str(folder)) if os.environ["ENDOTARGET_NEXTFLOW_VERIFY_ONLY"]=="true" else None,
 "historical_first_smoke_verification_failure":{"fixture":"/home/analyst/te/work/synthetic/20261013T095810Z-0dd4f9ec","cause":"harness expected synthetic but R explicitly emits synthetic_test_only","resolution":"correct verifier label and reverify preserved completed/cached outputs without recomputation","processes_already_successful":True},
 "python":platform.python_version(),
 "nextflow_version":subprocess.check_output(["nextflow","-version"],text=True).strip()}
out=root/"receipts"/("nextflow_synthetic_"+stamp+".json")
out.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
print("PASS full synthetic Nextflow chain, frozen-input rejection, 3/3 cached resume")
print(out)
PY
