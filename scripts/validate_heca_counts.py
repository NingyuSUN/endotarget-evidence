#!/usr/bin/env python3
"""Independent emitted-table/source-ledger consistency check; no model effects."""
import collections,csv,gzip,hashlib,json,pathlib
R=pathlib.Path("/home/analyst/te");O=R/"data/processed/heca"
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
 return h.hexdigest()
def table(p):
 with (gzip.open(p,"rt",newline="") if p.suffix==".gz" else p.open(newline="")) as f:return list(csv.DictReader(f,delimiter="\t"))
def counts(path,donors):
 rows=table(path);assert list(rows[0])==["gene_id"]+donors
 ids=[x["gene_id"] for x in rows];assert len(set(ids))==len(ids)
 by={}
 for row in rows:
  values=[int(row[d]) for d in donors];assert min(values)>=0 and max(values)<2**53
  by[row["gene_id"]]=values
 return by
def main():
 rec=json.loads((R/"receipts/heca_counts_preparation.json").read_text())
 assert rec["source_integer_counts"] and not rec["normalized"] and not rec["rounded"] and not rec["inverse_reconstruction_used"] and not rec["disease_effects_computed"]
 assert rec["execution_script_sha256_at_start"]==sha(R/"scripts/prepare_heca_counts.py")
 for p,h in rec["input_hashes"].items():assert sha(R/p)==h,p
 for x in rec["outputs"]:assert (R/x["path"]).stat().st_size==x["bytes"] and sha(R/x["path"])==x["sha256"],x["path"]
 donors=rec["selected_donors"];source=counts(O/"source_counts.tsv.gz",donors);stable=counts(O/"counts_stable.tsv.gz",donors)
 mapping=table(R/"manifests/heca_features_mapping.tsv");assert {r["source_feature_id"] for r in mapping}==set(source)
 expected={r["stable_gene_id"]:source[r["source_feature_id"]] for r in mapping if r["decision"]=="include"}
 assert stable==expected and len(stable)==sum(r["decision"]=="include" for r in mapping)
 ledger=table(O/"selected_donor_count_ledger.tsv");assert [r["donor_id"] for r in ledger]==donors
 cell=table(O/"raw_vs_author_H5_cell_totals.tsv.gz");assert len(cell)==52744 and len({x["cell_id"] for x in cell})==52744
 ct=collections.Counter();ctraw=collections.Counter();ctorig=collections.Counter()
 for x in cell:
  ct[x["donor_id"]]+=1;ctraw[x["donor_id"]]+=int(x["raw_36601_feature_UMIs"]);ctorig[x["donor_id"]]+=int(x["author_H5_original_n_counts"])
 totals=[]
 for j,row in enumerate(ledger):
  d=row["donor_id"];raw=sum(v[j] for v in source.values());kept=sum(v[j] for v in stable.values());excluded=int(row["mapping_excluded_UMIs"])
  assert raw==int(row["selected_raw_UMIs"])==ctraw[d] and raw==kept+excluded
  assert ctorig[d]==int(row["author_H5_original_UMIs"]) and ct[d]==int(row["matched_primary_nuclei"])==int(row["primary_nuclei"]) and ct[d]>=50
  totals.append({"donor_id":d,"source_UMIs":raw,"stable_retained_UMIs":kept,"mapping_excluded_UMIs":excluded,"selected_nuclei":ct[d],"author_H5_original_UMIs":ctorig[d]})
 downloads=json.loads((R/"receipts/heca_counts_raw_downloads.json").read_text())
 assert len(downloads["archives"])==22 and downloads["exact_archive_bytes"]==1582044098
 probe=json.loads((R/"receipts/heca_counts_raw_route_probe.json").read_text())
 probed=next(x for x in downloads["archives"] if x["sample_id"]==probe["probe"]["sample_id"]);assert probed["etag"]==probe["etag"]
 rawpath=R/probed["path"]
 with rawpath.open("rb") as f:assert hashlib.sha256(f.read(probe["downloaded_bytes"])).hexdigest()==probe["sha256"]
 result={"status":"PASS_TABLE_HASH_INTEGER_SOURCE_MAPPING_CELL_AND_COUNT_MASS_CONSISTENCY","source_genes":len(source),"stable_genes":len(stable),"donors":len(donors),"primary_nuclei":len(cell),"source_UMIs_total":sum(t["source_UMIs"] for t in totals),"stable_UMIs_total":sum(t["stable_retained_UMIs"] for t in totals),"donor_totals":totals,"preparation_receipt_sha256":sha(R/"receipts/heca_counts_preparation.json"),"disease_effects_computed":False,"validation_is_scientific_or_external_confirmation":False}
 (R/"receipts/heca_counts_validation.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps({k:v for k,v in result.items() if k!="donor_totals"},indent=2))
if __name__=="__main__":main()

