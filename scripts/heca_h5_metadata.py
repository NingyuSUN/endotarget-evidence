#!/usr/bin/env python3
"""Bounded obs-only anonymous HDF5 probe; expression datasets never opened."""
import argparse, collections, csv, hashlib, io, json, pathlib, platform, time
import h5py, numpy as np, requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=pathlib.Path("/home/analyst/te"); OUT=ROOT/"data/metadata/heca_h5"
URL="https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_nuclei.h5ad"
SIZE=8739186814; LIMIT=150*1024*1024; BLOCK=65536
def clean(x):
    if isinstance(x,bytes): return x.decode()
    if isinstance(x,np.generic): return clean(x.item())
    if isinstance(x,np.ndarray): return [clean(y) for y in x]
    if isinstance(x,(tuple,list)): return [clean(y) for y in x]
    if isinstance(x,dict): return {str(k):clean(v) for k,v in x.items()}
    return x
class BoundedRange(io.RawIOBase):
    def __init__(self):
        self.pos=0; self.started=time.monotonic(); self.session=requests.Session()
        self.session.mount("https://",HTTPAdapter(max_retries=Retry(total=2,backoff_factor=0.25,status_forcelist=[429,500,502,503,504])))
        self.log=[]; self.transferred=0; self.etag=None
        OUT.mkdir(parents=True,exist_ok=True); self.cache=OUT/"ranges"; self.cache.mkdir(exist_ok=True)
        self.cached=sum(x.stat().st_size for x in self.cache.glob("*.bin"))
        self.startfile=OUT/"probe_start.txt"
        if not self.startfile.exists(): self.startfile.write_text(str(min([p.stat().st_mtime for p in self.cache.glob("*.bin")] or [time.time()])))
        self.global_start=float(self.startfile.read_text())
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self,off,whence=0):
        p=off if whence==0 else self.pos+off if whence==1 else SIZE+off
        if not 0<=p<=SIZE: raise ValueError("seek outside source")
        self.pos=p; return p
    def chunk(self,start):
        end=min(start+BLOCK,SIZE)-1; p=self.cache/f"{start}-{end}.bin"
        if p.exists(): return p.read_bytes()
        if time.time()-self.global_start>=600: raise RuntimeError("600 second wall cap")
        if self.cached+end-start+1>LIMIT: raise RuntimeError("150 MiB aggregate source cap")
        t=time.monotonic()
        r=self.session.get(URL,headers={"Range":f"bytes={start}-{end}","Accept-Encoding":"identity"},timeout=(15,30),stream=True)
        if r.status_code!=206 or r.headers.get("Content-Range")!=f"bytes {start}-{end}/{SIZE}":
            r.close(); raise RuntimeError(f"invalid partial response: {r.status_code}, {r.headers.get('Content-Range')}")
        tag=r.headers.get("ETag")
        if self.etag is not None and tag!=self.etag: raise RuntimeError("ETag changed")
        self.etag=tag; data=r.raw.read(end-start+2); r.close()
        if len(data)!=end-start+1: raise RuntimeError("range length mismatch")
        p.write_bytes(data); self.transferred+=len(data); self.cached+=len(data)
        self.log.append({"start":start,"end":end,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),"etag":tag,"last_modified":r.headers.get("Last-Modified"),"seconds":round(time.monotonic()-t,3)})
        return data
    def read(self,n=-1):
        if n<0: raise RuntimeError("unbounded read forbidden")
        n=min(n,SIZE-self.pos)
        if n>8*1024*1024: raise RuntimeError("file read >8 MiB forbidden")
        data=bytearray()
        while n:
            start=(self.pos//BLOCK)*BLOCK; c=self.chunk(start)
            off=self.pos-start; take=min(n,len(c)-off); data.extend(c[off:off+take]); self.pos+=take; n-=take
        return bytes(data)
    def readinto(self,b):
        d=self.read(len(b)); b[:len(d)]=d; return len(d)
    def receipt(self,status):
        manifest=[{"file":p.name,"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(self.cache.glob("*.bin"))]
        old=ROOT/"receipts/heca_h5_ranges.json"; prev=json.loads(old.read_text()) if old.exists() else {}
        result={"source_url":URL,"expected_size":SIZE,"status":status,"network_budget_bytes":LIMIT,"wall_budget_seconds":600,"aggregate_exploration_elapsed_seconds":round(time.time()-self.global_start,3),"network_elapsed_seconds":round(time.monotonic()-self.started,3),"new_transferred_bytes":self.transferred,"total_cached_source_bytes":self.cached,"new_range_responses":self.log,"total_response_count":prev.get("total_response_count",0)+len(self.log),"cache_manifest_sha256":hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest(),"versions":{"python":platform.python_version(),"numpy":np.__version__,"h5py":h5py.__version__,"requests":requests.__version__},"expression_datasets_opened":[]}
        (OUT/"range_manifest.json").write_text(json.dumps(manifest,indent=2))
        old.write_text(json.dumps(result,indent=2)+"\n")
        (OUT/f"responses_{int(time.time())}.json").write_text(json.dumps(self.log,indent=2))
def read_obs(o,categories=None):
    if categories is not None:
        cats=clean(categories[:]); codes=o[:]
        if np.any(codes < -1) or np.any(codes >= len(cats)): raise ValueError("category code out of range")
        return np.array([cats[c] if c>=0 else None for c in codes],dtype=object)
    if isinstance(o,h5py.Group):
        cats=clean(o["categories"][:]); codes=o["codes"][:]
        return np.array([cats[c] if c>=0 else None for c in codes],dtype=object)
    return np.asarray(clean(o[:]),dtype=object)

PRIMARY={"eStromal","eStromal_MMPs","eStromal_cycling","dStromal_early","dStromal_mid","dStromal_late","Fibroblast_basalis"}
def audit_export():
    def rows(p): return list(csv.DictReader(p.open(newline=""),delimiter="\t" if p.suffix==".tsv" else ","))
    author=rows(ROOT/"data/metadata/heca_author/metadata-libraries_deconvoluted.csv")
    donors=rows(ROOT/"data/metadata/heca_author/metadata-donors.csv")
    bykey=collections.defaultdict(list)
    for a in author: bykey[a["Library_genotype"]].append(a)
    bydonor={d["Donor"]:d for d in donors}
    sdrf={s["source_name"]:s for s in rows(ROOT/"reports/clinical_gate_sdrf_rows.tsv")}
    fields=["library","genotype","Library_genotype","Group","Age","Hormonal treatment","Binary Stage","Stage","Biopsy_type","Tissue_sampled","discard_lowQC"]
    samples={}; labels=collections.Counter(); samplelabels=collections.Counter(); keys=set(); invalid=collections.Counter(); n=0
    with (OUT/"obs_selected.tsv").open(newline="") as fh:
        for x in csv.DictReader(fh,delimiter="\t"):
            n+=1; sample=x["sample"]; labels[x["celltype"]]+=1; samplelabels[(sample,x["celltype"])]+=1
            s=samples.setdefault(sample,{"nuclei":0,"primary_stroma":0,"HOXA13":0,"mesenchymal":0,"fields":{k:set() for k in fields}})
            s["nuclei"]+=1; s["primary_stroma"]+=x["celltype"] in PRIMARY
            s["HOXA13"]+=x["celltype"]=="HOXA13"; s["mesenchymal"]+=x["lineage"]=="Mesenchymal"
            for k in fields: s["fields"][k].add(x[k])
            idx=x.get("_index","")
            if idx:
                if idx in keys: invalid["duplicate_library_barcode"]+=1
                keys.add(idx)
                if not idx.startswith(x["library"]+"_"): invalid["index_library_prefix"]+=1
                barcode=idx[len(x["library"])+1:]
                if len(barcode)!=16 or set(barcode)-set("ACGT"): invalid["barcode_format"]+=1
            else: invalid["index_not_exported"]+=1
            if sample!=x["Library_genotype"]: invalid["sample_vs_library_genotype"]+=1
    out=[]; counters=collections.Counter(); dc=collections.defaultdict(collections.Counter)
    for sample,s in sorted(samples.items()):
        vals={k:";".join(sorted(v)) for k,v in s["fields"].items()}
        arows=bykey[sample]; a=arows[0] if len(arows)==1 else {}
        d=bydonor.get(vals["genotype"],{}); sr=sdrf.get(sample,{})
        aa=a.get("Donor","UNKNOWN"); decl=sr.get("declared_individual","UNKNOWN")
        h5d="endometriosis" if vals["Group"]=="Endo_Superficial" else "normal" if vals["Group"]=="Ctrl_Superficial" else "UNKNOWN"
        counters["sample_author_unique"]+=len(arows)==1
        counters["h5_genotype_author_donor_match"]+=vals["genotype"]==aa
        counters["h5_genotype_sdrf_declared_match"]+=vals["genotype"]==decl
        counters["h5_disease_sdrf_disease_match"]+=h5d==sr.get("disease")
        counters["h5_phase_sdrf_phase_match"]+=vals["Binary Stage"]==sr.get("phase")
        counters["h5_group_author_group_match"]+=vals["Group"]==a.get("Group")
        if any(len(v)!=1 for v in s["fields"].values()): counters["samples_variable_clinical_identity_fields"]+=1
        rawage=d.get("Age","")
        try: ageeq=float(vals["Age"])==float(rawage)
        except ValueError: ageeq=False
        counters["h5_age_author_age_match"]+=ageeq
        if vals["Hormonal treatment"]!="nan" and vals["Binary Stage"]=="Proliferative": counters["proliferative_explicit_treatment_conflict"]+=1
        for key in ["nuclei","primary_stroma","HOXA13","mesenchymal"]: dc[vals["genotype"]][key]+=s[key]
        row={"h5_sample":sample,**vals,**{k:s[k] for k in ["nuclei","primary_stroma","HOXA13","mesenchymal"]},
             "author_deconv_Donor":aa,"author_deconv_genotyped":a.get("genotyped","UNKNOWN"),"author_deconv_StudyName":a.get("StudyName","UNKNOWN"),"sdrf_declared_individual":decl,"sdrf_disease":sr.get("disease","UNKNOWN"),"sdrf_phase":sr.get("phase","UNKNOWN"),"sdrf_identity_matches":vals["genotype"]==decl,
             "author_Age":rawage,"author_Hormonal_treatment":d.get("Hormonal treatment","UNKNOWN"),"author_Binary_Stage":d.get("Binary Stage","UNKNOWN"),"author_Binary_Stage_freezed":d.get("Binary Stage_freezed","UNKNOWN"),"author_histology_phase":d.get("MensPhase3_pathologist","UNKNOWN"),"author_expression_phase":d.get("Binary Stage from scRNAseq data","UNKNOWN"),"author_Endometrium_histology":d.get("Endometrium histology","UNKNOWN"),
             "clinical_eligibility":"NOT_ADJUDICATED","primary_mapping_status":"PROPOSED_AUTHOR_LABEL_ONLY"}
        out.append(row)
    def write_tsv(path,data):
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t"); w.writeheader(); w.writerows(data)
    write_tsv(ROOT/"reports/heca_h5_sample_crosswalk.tsv",out)
    write_tsv(ROOT/"reports/heca_h5_donor_counts.tsv",[{"h5_genotype":k,**v,"sample_count":sum(r["genotype"]==k for r in out),"clinical_eligibility":"NOT_ADJUDICATED"} for k,v in sorted(dc.items())])
    write_tsv(ROOT/"reports/heca_h5_author_label_counts.tsv",[{"h5_sample":s,"author_celltype":l,"nuclei":v,"primary_mapping_proposed":l in PRIMARY} for (s,l),v in sorted(samplelabels.items())])
    candidates=rows(ROOT/"manifests/clinical_candidates.tsv")
    candidateids={x["declared_individual"] for x in candidates}
    mismatched=[r for r in out if not r["sdrf_identity_matches"]]
    prolif=[r for r in out if r["Binary Stage"]=="Proliferative"]
    groups=collections.defaultdict(set)
    for r in prolif: groups[r["Group"]].add(r["genotype"])
    files=[ROOT/"scripts/heca_h5_metadata.py",OUT/"obs_selected.tsv",ROOT/"reports/heca_h5_schema.json",ROOT/"reports/heca_h5_sample_crosswalk.tsv",ROOT/"reports/heca_h5_donor_counts.tsv",ROOT/"reports/heca_h5_author_label_counts.tsv",ROOT/"data/metadata/heca_author/metadata-libraries_deconvoluted.csv",ROOT/"data/metadata/heca_author/metadata-donors.csv",ROOT/"reports/clinical_gate_sdrf_rows.tsv"]
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    result={"status":"author_annotation_crosswalk_verified_clinical_eligibility_not_adjudicated","n_obs":n,"samples":len(samples),"genotype_ids":len(dc),"libraries":len({r["library"] for r in out}),"barcode_unique":len(keys),"validation_errors":dict(invalid),"comparisons":dict(counters),"sdrf_identity_mismatch_samples":len(mismatched),"sdrf_candidate_ids":len(candidateids),"candidate_declared_ids_affected":sorted(candidateids & {r["sdrf_declared_individual"] for r in mismatched}),"h5_proliferative_samples":len(prolif),"h5_proliferative_unique_genotypes_by_group":{k:sorted(v) for k,v in groups.items()},"author_celltype_counts":dict(labels),"primary_author_labels_proposed":sorted(PRIMARY),"primary_stroma_total":sum(v["primary_stroma"] for v in dc.values()),"HOXA13_separate_total":labels["HOXA13"],"mesenchymal_total":sum(v["mesenchymal"] for v in dc.values()),"expression_datasets_opened":[],"clinical_eligibility":"NOT_ADJUDICATED","hashes":hashes}
    (ROOT/"receipts/heca_h5_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    return result
def self_test():
    assert "uSMCs" not in PRIMARY and all(not x.startswith("ePV") for x in PRIMARY)
    r=BoundedRange()
    for f in [lambda:r.read(-1),lambda:r.seek(SIZE+1),lambda:r.read(8*1024*1024+1)]:
        try: f()
        except (ValueError,RuntimeError): pass
        else: raise AssertionError("negative boundary accepted")
    fresh=next(s for s in range(0,SIZE,BLOCK) if not (r.cache/f"{s}-{min(s+BLOCK,SIZE)-1}.bin").exists())
    r.global_start=time.time(); r.cached=LIMIT
    try: r.chunk(fresh)
    except RuntimeError as e: assert "150 MiB" in str(e)
    else: raise AssertionError("network byte budget accepted")
    r.cached=0; r.global_start=time.time()-601
    try: r.chunk(fresh)
    except RuntimeError as e: assert "600 second" in str(e)
    else: raise AssertionError("wall budget accepted")
    assert r.transferred==0
    print("PASS: byte budget, wall budget, unbounded read, out-of-source seek, >8MiB read rejected without network; mural/perivascular labels excluded")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--export",nargs="*"); ap.add_argument("--audit",action="store_true"); ap.add_argument("--self-test",action="store_true"); a=ap.parse_args()
    if a.self_test: self_test(); return
    if a.audit: audit_export(); return
    r=BoundedRange(); status="failed"
    try:
        with h5py.File(r,"r") as f:
            obs=f["obs"]; info={"source_url":URL,"source_size":SIZE,"obs_attributes":clean(dict(obs.attrs)),"columns":{}}
            for k in obs.keys():
                o=obs[k]
                if isinstance(o,h5py.Group):
                    if "categories" in o and "codes" in o:
                        info["columns"][k]={"encoding":clean(dict(o.attrs)),"categories":clean(o["categories"][:]),"codes_shape":list(o["codes"].shape),"codes_dtype":str(o["codes"].dtype)}
                    else:
                        info["columns"][k]={"encoding":clean(dict(o.attrs)),"children":{kk:{"shape":list(o[kk].shape),"dtype":str(o[kk].dtype),"head":clean(o[kk][:5])} for kk in o}}
                else:
                    info["columns"][k]={"shape":list(o.shape),"dtype":str(o.dtype),"head":clean(o[:5])}
                    if "__categories" in obs and k in obs["__categories"]:
                        info["columns"][k]["categories"]=clean(obs["__categories"][k][:])
            (OUT/"obs_schema.json").write_text(json.dumps(info,indent=2)+"\n")
            (ROOT/"reports/heca_h5_schema.json").write_text(json.dumps(info,indent=2)+"\n")
            print(json.dumps(info,indent=2),flush=True)
            if a.export is not None:
                chosen=a.export
                if not chosen or any(k not in obs for k in chosen): raise ValueError("explicit existing obs fields required")
                arrays={k:read_obs(obs[k],obs["__categories"][k] if "__categories" in obs and k in obs["__categories"] else None) for k in chosen}; n=len(next(iter(arrays.values())))
                if any(len(v)!=n for v in arrays.values()): raise ValueError("obs lengths differ")
                with (OUT/"obs_selected.tsv").open("w") as fh:
                    w=csv.writer(fh,delimiter="\t"); w.writerow(chosen); w.writerows(zip(*(arrays[k] for k in chosen)))
                summary={"n_obs":n,"columns":chosen,"marginals":{k:dict(collections.Counter(str(x) for x in v)) for k,v in arrays.items() if len(set(str(x) for x in v))<=1000}}
                (OUT/"obs_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
                print(json.dumps({"export_rows":n,"columns":chosen}))
            status="ok"
    finally: r.receipt(status)
if __name__=="__main__": main()

