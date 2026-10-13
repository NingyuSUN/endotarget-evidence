#!/usr/bin/env python3
"""Pre-effect counts-layer audits. Default is offline; .X is known normalized-log data.
Never use inverse feasibility as approved counts. Additional numeric probe is bounded.
All commands must execute from /home/analyst/te on SRV01.
"""
import argparse
def cached_features():
    import sys,pathlib,json,csv,hashlib,collections,h5py
    sys.path.insert(0,str(pathlib.Path(__file__).parent));import heca_h5_metadata as m
    m.OUT=pathlib.Path("/home/analyst/te/data/metadata/heca_h5/heca_count_probe");m.BLOCK=262144
    r=m.BoundedRange();m.LIMIT=r.cached
    out=m.OUT/"features_cached.tsv"
    with h5py.File(r,"r") as f:
     v=f["var"];names=m.clean(v["_index"][:]);sets=[set() for x in names];used=[];failed=[]
     for k in sorted([k for k in v if k.startswith("gene_ids-")],key=lambda x:int(x.split("-")[-1])):
      try:a=m.read_obs(v[k],v["__categories"][k])
      except RuntimeError:failed.append(k);continue
      used.append(k)
      for i,x in enumerate(a):
       if x and str(x)!="nan":sets[i].add(str(x))
     with out.open("w") as h:
      w=csv.writer(h,delimiter="\t");w.writerow(["feature_index","source_symbol","source_ensembl_ids","source_ensembl_id_count","library_gene_id_columns_read","library_gene_id_column_count"])
      w.writerows((i,n,";".join(sorted(s)),len(s),";".join(used),51) for i,(n,s) in enumerate(zip(names,sets)))
     result={"n_features":len(names),"unique_symbols":len(set(names)),"gene_id_columns_read":used,"gene_id_columns_unread":failed,"source_id_set_sizes":dict(collections.Counter(len(s) for s in sets)),"path":str(out),"sha256":hashlib.sha256(out.read_bytes()).hexdigest(),"new_network_bytes":r.transferred}
     pathlib.Path("/home/analyst/te/receipts/heca_counts_cached_features.json").write_text(json.dumps(result,indent=2)+"\n")
     print(json.dumps(result,indent=2))
    
def numeric_probe():
    import sys,pathlib,json,hashlib,time,csv,h5py,numpy as np,resource
    sys.path.insert(0,"scripts");import heca_h5_metadata as m
    m.OUT=pathlib.Path("data/metadata/heca_h5/heca_count_probe");m.BLOCK=262144
    r=m.BoundedRange();initial=r.cached;m.LIMIT=initial+10*1024*1024;r.global_start=time.time()
    sel=np.load("work/final_clinical/heca_primary_h5_obs_indices.npy"); chosen=sel[:5].tolist()
    result={"source_url":m.URL,"source_size":m.SIZE,"initial_cache_bytes":initial,"additional_probe_budget_bytes":10*1024*1024,"status":"NOT_COMPLETED"}
    try:
     with h5py.File(r,"r") as f:
      x=f["X"];ptr=x["indptr"][:]
      assert ptr[0]==0 and np.all(ptr[1:]>=ptr[:-1]) and ptr[-1]==x["data"].shape[0]
      np.save(m.OUT/"X_indptr.npy",ptr)
      sample=[]
      for start in [0,len(x["data"])//2,len(x["data"])-1024]:
       a=x["data"][start:start+1024]
       sample.append({"start":int(start),"n":len(a),"finite":bool(np.isfinite(a).all()),"nonnegative":bool((a>=0).all()),"integer":bool((a==np.floor(a)).all()),"min":float(a.min()),"max":float(a.max()),"head":a[:10].tolist(),"sha256":hashlib.sha256(a.tobytes()).hexdigest()})
      totals=[]
      for row in chosen:
       s,e=int(ptr[row]),int(ptr[row+1]);a=x["data"][s:e]; ix=x["indices"][s:e];observed=float(f["obs/n_counts"][row])
       totals.append({"obs_row":row,"nnz":e-s,"sum":float(a.astype(np.float64).sum()),"author_obs_n_counts":observed,"matches":float(a.astype(np.float64).sum())==observed,"finite":bool(np.isfinite(a).all()),"nonnegative":bool((a>=0).all()),"integer":bool((a==np.floor(a)).all()),"indices_in_bounds":bool((ix>=0).all() and (ix<33578).all()),"unique_gene_indices":len(set(ix.tolist()))==len(ix)})
      chunk=x["data"].chunks[0]; ranges=sorted(set((int(ptr[i])//chunk,int(max(ptr[i],ptr[i+1]-1))//chunk) for i in sel))
      touched=set()
      for s,e in ranges:touched.update(range(s,e+1))
      allchunks=(len(x["data"])+chunk-1)//chunk
      sizes={}
      for name in ["data","indices"]:
       ds=x[name];bs=[ds.id.get_chunk_info_by_coord((int(i)*chunk,)).size for i in touched]
       sizes[name]={"selected_chunks":len(touched),"all_chunks":allchunks,"selected_chunk_storage_bytes":sum(bs),"full_storage_bytes":ds.id.get_storage_size()}
      result.update({"status":"bounded_numeric_and_geometry_verified","selected_obs_rows":len(sel),"nnz_selected":int(np.sum(ptr[sel+1].astype(np.int64)-ptr[sel].astype(np.int64))),"probes":sample,"row_totals":totals,"indptr_path":str(m.OUT/"X_indptr.npy"),"indptr_sha256":hashlib.sha256(ptr.tobytes()).hexdigest(),"selected_chunk_geometry":sizes,"etag":r.etag})
    finally:
     result["network"]={"new_bytes":r.transferred,"wall_seconds":time.time()-r.global_start,"source_cache_bytes":r.cached,"peak_RSS_KiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"range_responses":r.log}
     pathlib.Path("receipts/heca_counts_numeric_probe.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="network"},indent=2))
    
def inverse_feasibility():
    import sys,pathlib,time,json,h5py,numpy as np,hashlib
    sys.path.insert(0,"scripts");import heca_h5_metadata as m
    m.OUT=pathlib.Path("data/metadata/heca_h5/heca_count_probe");m.BLOCK=262144;r=m.BoundedRange();m.LIMIT=r.cached
    rows=[x["obs_row"] for x in json.load(open("receipts/heca_counts_numeric_probe.json"))["row_totals"]]
    out=[]
    with h5py.File(r,"r") as f:
     ptr=np.load(m.OUT/"X_indptr.npy")
     for row in rows:
      a=f["X/data"][int(ptr[row]):int(ptr[row+1])].astype(np.float64); n=float(f["obs/n_counts"][row])
      estimate=np.expm1(a)*n/10000.0;nearest=np.rint(estimate);delta=np.abs(estimate-nearest)
      reconstructed=np.log1p(nearest*10000/n).astype(np.float32)
      out.append({"obs_row":row,"nonzero_values":len(a),"normalization_target_assumption":10000,"nearest_integer_distance_max":float(delta.max()),"nearest_integer_distance_median":float(np.median(delta)),"raw_estimate_sum":float(estimate.sum()),"nearest_sum":float(nearest.sum()),"author_original_n_counts":n,"nearest_sum_equals_original_total":bool(nearest.sum()==n),"reconstructed_log_float32_matches_fraction":float(np.mean(reconstructed==a.astype(np.float32))),"zero_or_negative_reconstructed":int(np.sum(nearest<=0)),"status":"FIVE_ROW_FEASIBILITY_ONLY_NOT_ACCEPTED_COUNTS"})
    rec={"method_proposed":"expm1(log_normalized_X)*obs_original_n_counts/10000;nearest_integer_proposal_requires_explicit_method_decision","source_notebook":"main_analysis/rna3_clustering_nuclei.ipynb code cell20 at pinned244c282c20c3efd53707e3eb9f58bbd1cbd56ad7","new_network_bytes":r.transferred,"rows":out,"full_dataset_reconstruction_verified":False,"counts_accepted":False}
    pathlib.Path("receipts/heca_counts_inverse_feasibility.json").write_text(json.dumps(rec,indent=2)+"\n");print(json.dumps(rec,indent=2))
    

def self_test():
    import numpy as np
    raw=np.array([1,2,7,19],dtype=np.float32)
    logged=np.log1p(raw/np.float32(raw.sum()/10000))
    assert not np.all(logged==np.floor(logged))
    assert np.any(raw!=logged)
    recovered=np.expm1(logged.astype(np.float64))*float(raw.sum())/10000
    assert np.max(np.abs(recovered-np.rint(recovered)))<1e-3
    assert np.rint(recovered).sum()==raw.sum()
    print("PASS: normalized/log values fail raw-integer gate; inverse fixture is feasibility only")
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["cached-features","numeric-probe","inverse-feasibility","self-test"],default="self-test")
    a=ap.parse_args()
    {"cached-features":cached_features,"numeric-probe":numeric_probe,"inverse-feasibility":inverse_feasibility,"self-test":self_test}[a.mode]()
if __name__=="__main__": main()

