#!/usr/bin/env python3
"""Public raw-UMI archive intake and fixed-author-cell donor aggregation; no DE."""
import argparse,collections,csv,gzip,hashlib,inspect,io,json,os,pathlib,platform,resource,shutil,tarfile,time
EXECUTION_SCRIPT_SHA=None
from concurrent.futures import ThreadPoolExecutor
import numpy as np,requests
ROOT=pathlib.Path("/home/analyst/te"); OUT=ROOT/"data/processed/heca"
RAW=ROOT/"data/metadata/heca_h5/heca_count_raw_archives"
DISK_RESERVE=20*1024**3; RSS_LIMIT_KIB=12*1024**2
def sha(path):
    h=hashlib.sha256()
    with pathlib.Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def tsv(path):return list(csv.DictReader(pathlib.Path(path).open(newline=""),delimiter="\t"))
def write_tsv(path,rows,fields):
    path=pathlib.Path(path);path.parent.mkdir(exist_ok=True,parents=True)
    if str(path).endswith(".gz"):
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz:
                with io.TextIOWrapper(gz,encoding="utf-8",newline="") as f:
                    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
    else:
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def guard(extra=0):
    free=shutil.disk_usage(ROOT).free
    if free-extra<DISK_RESERVE:raise RuntimeError("20 GiB disk reserve would be violated")
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>RSS_LIMIT_KIB:raise RuntimeError("12 GiB RSS guard exceeded")
    return free
def inventory():return json.loads((ROOT/"reports/heca_counts_raw_inventory.json").read_text())
def fetch(item):
    guard(item["expected_bytes"]);RAW.mkdir(exist_ok=True,parents=True)
    p=RAW/item["path"];receipt=RAW/(item["path"]+".receipt.json")
    if p.exists():
        if not receipt.exists():raise ValueError("Existing archive lacks receipt")
        old=json.loads(receipt.read_text())
        if p.stat().st_size!=item["expected_bytes"] or sha(p)!=old["sha256"]:raise ValueError("Cached archive differs from pinned receipt")
        return {**old,"reused":True}
    tmp=p.with_suffix(p.suffix+".partial")
    if tmp.exists():raise ValueError("Partial archive needs explicit resume audit")
    started=time.monotonic();r=requests.get(item["url"],stream=True,timeout=(15,90))
    r.raise_for_status();content=r.headers.get("Content-Length")
    if content is not None and int(content)!=item["expected_bytes"]:r.close();raise ValueError("Archive Content-Length changed")
    received=0;hasher=hashlib.sha256()
    try:
        with tmp.open("wb") as f:
            for b in r.iter_content(1024*1024):
                received+=len(b)
                if received>item["expected_bytes"]:raise ValueError("Archive exceeds official size")
                hasher.update(b);f.write(b);guard()
        if received!=item["expected_bytes"]:raise ValueError("Archive shorter than official size")
        tmp.replace(p)
    finally:r.close()
    rec={**item,"path":str(p.relative_to(ROOT)),"bytes":received,"sha256":hasher.hexdigest(),"etag":r.headers.get("ETag"),"last_modified":r.headers.get("Last-Modified"),"final_url":r.url,"http_status":r.status_code,"wall_seconds":round(time.monotonic()-started,3),"reused":False}
    receipt.write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps({"downloaded":item["sample_id"],"bytes":received,"seconds":rec["wall_seconds"]}),flush=True)
    return rec
def download():
    items=inventory();guard(sum(x["expected_bytes"] for x in items))
    started=time.monotonic()
    with ThreadPoolExecutor(max_workers=2) as ex:records=list(ex.map(fetch,items))
    rec={"status":"PUBLIC_RAW_ARCHIVES_DOWNLOADED_NOT_YET_COUNT_VALIDATED","source_inventory_sha256":sha(ROOT/"reports/heca_counts_raw_inventory.json"),"archives":records,"exact_archive_bytes":sum(x["expected_bytes"] for x in items),"actual_new_bytes":sum(x["bytes"] for x in records if not x["reused"]),"wall_seconds":round(time.monotonic()-started,3),"disk_free_bytes_after":guard(),"versions":{"python":platform.python_version(),"numpy":np.__version__,"requests":requests.__version__},"source_scope":"22 exact author-selected public processed raw-count archives; no controlled reads"}
    (ROOT/"receipts/heca_counts_raw_downloads.json").write_text(json.dumps(rec,indent=2)+"\n")
    return records

def join_barcodes(raw,selected,library):
    if len(selected)!=len(set(selected)):raise ValueError("Duplicate selected library+barcode")
    normalized=[]
    for barcode in raw:
        if barcode.startswith(library+"_"):
            token=barcode[len(library)+1:]
        else:token=barcode.split("-")[0]
        if len(token)!=16 or set(token)-set("ACGT"):raise ValueError("Unrecognized raw barcode structure")
        normalized.append(library+"_"+token)
    if len(normalized)!=len(set(normalized)):raise ValueError("Raw barcode normalization creates identity ambiguity")
    missing=set(selected)-set(normalized)
    if missing:raise ValueError(f"Selected author barcodes absent from raw matrix: {len(missing)}")
    use=set(selected)
    return np.array([x in use for x in normalized],dtype=bool),normalized
def check_features(observed,reference):
    if observed!=reference:raise ValueError("Archive source feature identities/order differ; no guessed realignment")
    if len({x[0] for x in observed})!=len(observed):raise ValueError("Duplicate original source gene IDs")
def aggregate_mtx(stream,n_features,selected_mask):
    header=stream.readline().strip()
    if header!=b"%%MatrixMarket matrix coordinate integer general":raise ValueError("Expected source integer MatrixMarket counts")
    dims=stream.readline()
    while dims.startswith(b"%"):dims=stream.readline()
    ng,nc,declared=[int(x) for x in dims.split()]
    if ng!=n_features or nc!=len(selected_mask):raise ValueError("Matrix dimensions differ from source features/barcodes")
    if nc>100000 or ng>100000:raise ValueError("Unexpected source matrix size")
    mask=np.array(selected_mask,dtype=bool);gene=np.zeros(ng,dtype=np.float64);cell=np.zeros(nc,dtype=np.float64)
    seen=0;total=0;selected_total=0;last=-1;buf=b"";iterations=0
    while True:
        block=stream.read(2*1024*1024)
        if not block:
            part=buf;buf=b""
            if not part:break
        else:
            buf+=block;cut=buf.rfind(b"\n")
            if cut<0:
                if len(buf)>8*1024*1024:raise ValueError("Unbounded MatrixMarket line")
                continue
            part=buf[:cut+1];buf=buf[cut+1:]
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            a=np.fromstring(part.decode("ascii"),dtype=np.int64,sep=" ")
        if a.size%3:raise ValueError("Malformed matrix triplets")
        a=a.reshape(-1,3);g=a[:,0]-1;c=a[:,1]-1;v=a[:,2]
        if np.any(g<0) or np.any(g>=ng) or np.any(c<0) or np.any(c>=nc) or np.any(v<0):raise ValueError("Invalid sparse coordinate or negative count")
        keys=c*ng+g
        if len(keys) and (keys[0]<=last or np.any(keys[1:]<=keys[:-1])):raise ValueError("Duplicate/unordered source coordinates require source audit")
        if len(keys):last=int(keys[-1])
        use=mask[c];vv=v[use];gg=g[use]
        gene+=np.bincount(gg,weights=vv,minlength=ng)
        cell+=np.bincount(c,weights=v,minlength=nc)
        seen+=len(a);total+=int(v.sum(dtype=np.int64));selected_total+=int(vv.sum(dtype=np.int64))
        if np.any(gene>=2**53) or np.any(cell>=2**53):raise ValueError("Exact float64 integer range exceeded")
        iterations+=1
        if iterations%100==0:guard()
    if seen!=declared:raise ValueError("Matrix nnz differs from declared header")
    if int(gene.sum())!=selected_total or int(cell[mask].sum())!=selected_total or int(cell.sum())!=total:raise ValueError("Count mass conservation failed")
    if np.any(gene!=np.floor(gene)) or np.any(cell!=np.floor(cell)):raise ValueError("Noninteger aggregation")
    return {"gene_counts":gene,"cell_counts":cell,"declared_nnz":declared,"observed_nnz":seen,"all_total":total,"selected_total":selected_total,"source_shape":[ng,nc]}
def archive_metadata(item):
    p=RAW/item["path"]
    if not p.exists() or p.stat().st_size!=item["expected_bytes"]:raise ValueError("Exact raw archive unavailable")
    rec=json.loads((RAW/(item["path"]+".receipt.json")).read_text())
    if sha(p)!=rec["sha256"]:raise ValueError("Raw archive hash changed")
    # Reading all outer gzip bytes validates its CRC without extracting members.
    with gzip.open(p,"rb") as f:
        for _ in iter(lambda:f.read(1024*1024),b""):pass
    with tarfile.open(p,"r:gz") as tar:
        members=tar.getmembers()
        for member in members:
            path=pathlib.PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts or member.issym() or member.islnk():raise ValueError("Unsafe or ambiguous tar member")
        lookup={}
        for tail in ["matrix.mtx.gz","barcodes.tsv.gz","features.tsv.gz"]:
            hits=[x for x in members if x.isfile() and pathlib.PurePosixPath(x.name).name==tail]
            if len(hits)!=1:raise ValueError("Ambiguous source archive member "+tail)
            if pathlib.PurePosixPath(hits[0].name).parent.name!=item["sample_id"]:raise ValueError("Archive internal sample identity differs")
            lookup[tail]=hits[0].name
        with gzip.open(tar.extractfile(lookup["features.tsv.gz"]),"rt") as f:features=list(csv.reader(f,delimiter="\t"))
        with gzip.open(tar.extractfile(lookup["barcodes.tsv.gz"]),"rt") as f:barcodes=[x.strip() for x in f]
    if any(len(x)!=3 or x[2]!="Gene Expression" for x in features):raise ValueError("Unexpected source feature modality")
    return features,barcodes,lookup,rec
def author_selected():
    path=ROOT/"manifests/final_heca_clinical_selected.tsv"
    if sha(path)!="1a0ba5813046f9cb6e18862e17ca7fe64c31a39605d048309f1e3522447046cc":raise ValueError("Final clinical selection changed")
    clinical=tsv(path);items=inventory()
    if [x["sample_id"] for x in clinical]!=[x["sample_id"] for x in items] or len({x["donor_id"] for x in clinical})!=len(clinical):raise ValueError("Clinical ledger and raw inventory differ")
    primary={"eStromal_MMPs","eStromal","eStromal_cycling","dStromal_early","dStromal_mid","dStromal_late","Fibroblast_basalis"}
    by={r["sample_id"]:[] for r in clinical};pos=[];byrow={}
    with (ROOT/"data/metadata/heca_h5/obs_selected.tsv").open() as f:
        for i,row in enumerate(csv.DictReader(f,delimiter="\t")):
            if row["sample"] in by and row["celltype"] in primary and row["discard_lowQC"]=="no":
                by[row["sample"]].append(row["_index"]);pos.append(i);byrow[row["_index"]]=i
    expected=ROOT/"work/final_clinical/heca_primary_h5_obs_indices.npy"
    if sha(expected)!="d02ed9b8cee2c789c4f7aec7a3792a8f9444307163531933b632845cecf864ef" or not np.array_equal(np.asarray(pos),np.load(expected)):raise ValueError("Own fixed lineage differs from selected original H5 indices")
    for row in clinical:
        if len(by[row["sample_id"]])!=int(row["primary_nuclei"]) or len(by[row["sample_id"]])<50:raise ValueError("Selected primary nuclei minimum mismatch")
    return clinical,by,byrow

def checkpoint_pins(item,features,selected,download):
    return {"archive_sha256":download["sha256"],"selection_sha256":sha(ROOT/"manifests/final_heca_clinical_selected.tsv"),"selected_barcodes_sha256":hashlib.sha256(json.dumps(sorted(selected),separators=(",",":")).encode()).hexdigest(),"source_features_sha256":hashlib.sha256(json.dumps(features,separators=(",",":")).encode()).hexdigest(),"parser_sha256":hashlib.sha256(inspect.getsource(aggregate_mtx).encode()).hexdigest(),"checkpoint_schema_version":2,"sample_id":item["sample_id"],"donor_id":item["donor_id"]}
def validate_checkpoint(meta,saved,pins,npz,ng,mask):
    if any(meta.get(k)!=v for k,v in pins.items()):raise ValueError("Checkpoint source/barcode/feature/parser pin mismatch")
    if meta.get("npz_sha256")!=sha(npz):raise ValueError("Checkpoint NPZ hash changed")
    gene=saved["gene_counts"];cell=saved["cell_counts"]
    if gene.shape!=(ng,) or cell.shape!=(len(mask),):raise ValueError("Checkpoint array dimensions differ")
    for a in [gene,cell]:
        if not np.isfinite(a).all() or np.any(a<0) or np.any(a!=np.floor(a)) or np.any(a>=2**53):raise ValueError("Checkpoint numeric integrity failed")
    if meta["observed_nnz"]!=meta["declared_nnz"] or meta["source_shape"]!=[ng,len(mask)]:raise ValueError("Checkpoint source geometry mismatch")
    if int(gene.sum())!=meta["selected_total"] or int(cell[mask].sum())!=meta["selected_total"] or int(cell.sum())!=meta["all_total"]:raise ValueError("Checkpoint count mass mismatch")
    return True

def aggregate(wait_for_downloads=False):
    started=time.monotonic();OUT.mkdir(exist_ok=True,parents=True);check=OUT/"checkpoints";check.mkdir(exist_ok=True)
    clinical,by,byrow=author_selected();donors=[x["donor_id"] for x in clinical]
    reference=None;columns=[];ledgers=[];recorded=[];cell_comparison=[]
    original_h5_counts=np.load(ROOT/"data/metadata/heca_h5/heca_count_probe/original_obs_n_counts.npy")
    for item,row in zip(inventory(),clinical):
        tick=time.monotonic()
        while wait_for_downloads and (not (RAW/item["path"]).exists() or not (RAW/(item["path"]+".receipt.json")).exists()):
            if time.monotonic()-tick>1800:raise RuntimeError("Archive download wait exceeded 30 minutes")
            time.sleep(5)
        features,barcodes,members,download=archive_metadata(item)
        if reference is None:reference=features
        check_features(features,reference)
        mask,normalized=join_barcodes(barcodes,by[item["sample_id"]],item["library_id"])
        if int(mask.sum())!=int(row["primary_nuclei"]):raise ValueError("Selected source barcode count differs from clinical minimum ledger")
        npz=check/(item["sample_id"]+".npz");jsonpath=check/(item["sample_id"]+".json")
        pins=checkpoint_pins(item,features,by[item["sample_id"]],download)
        existing_meta=json.loads(jsonpath.read_text()) if jsonpath.exists() else {}
        reuse=npz.exists() and existing_meta.get("checkpoint_schema_version")==2
        if reuse:
            with np.load(npz) as saved:
                validate_checkpoint(existing_meta,saved,pins,npz,len(features),mask)
                result={**existing_meta,"gene_counts":saved["gene_counts"].copy(),"cell_counts":saved["cell_counts"].copy()}
        else:
            with tarfile.open(RAW/item["path"],"r:gz") as tar:
                with gzip.open(tar.extractfile(members["matrix.mtx.gz"]),"rb") as matrix:
                    result=aggregate_mtx(matrix,len(features),mask)
            meta={k:v for k,v in result.items() if k not in ["gene_counts","cell_counts"]}
            meta.update(pins)
            np.savez_compressed(npz,gene_counts=result["gene_counts"],cell_counts=result["cell_counts"])
            meta["npz_sha256"]=sha(npz);meta["execution_script_sha256"]=EXECUTION_SCRIPT_SHA
            jsonpath.write_text(json.dumps(meta,indent=2)+"\n")
            with np.load(npz) as saved:validate_checkpoint(meta,saved,pins,npz,len(features),mask)
        columns.append(result["gene_counts"])
        for index,key in enumerate(normalized):
            if mask[index]:
                original=int(original_h5_counts[byrow[key]]);raw_count=int(result["cell_counts"][index])
                cell_comparison.append({"donor_id":item["donor_id"],"sample_id":item["sample_id"],"cell_id":key,"original_h5_row":byrow[key],"raw_36601_feature_UMIs":raw_count,"author_H5_original_n_counts":original,"raw_minus_author_H5":raw_count-original})
        ledger={**row,"raw_source_cells":len(barcodes),"matched_primary_nuclei":int(mask.sum()),"original_source_gene_rows":len(features),"selected_raw_UMIs":result["selected_total"],"raw_total_all_source_cells":result["all_total"],"matrix_nnz":result["observed_nnz"],"source_archive_sha256":download["sha256"],"source_integer_counts":True,"normalized":False,"inverted_or_rounded":False,"barcode_complete_match":True}
        ledgers.append(ledger)
        recorded.append({"sample_id":item["sample_id"],"donor_id":item["donor_id"],"archive_sha256":download["sha256"],"members":members,"source_feature_sha256":hashlib.sha256(json.dumps(features).encode()).hexdigest(),"wall_seconds":round(time.monotonic()-tick,3),"source_shape":result["source_shape"]})
        print(json.dumps({"aggregated":item["sample_id"],"selected_nuclei":int(mask.sum()),"UMIs":result["selected_total"],"seconds":recorded[-1]["wall_seconds"]}),flush=True);guard()
    matrix=np.asarray(columns,dtype=np.float64).T
    source_ids=[x[0] for x in reference]
    from prepare_discovery_counts import map_features,table
    annotation=json.loads((ROOT/"config/gene_mapping_reference.json").read_text());annpath=ROOT/annotation["path"]
    if sha(annpath)!=annotation["sha256"]:raise ValueError("Pinned HGNC annotation differs")
    mappings=map_features(source_ids,table(annpath))
    for row,feature in zip(mappings,reference):row["source_symbol"]=feature[1];row["source_symbol_agrees_HGNC"]=feature[1]==row["approved_symbol"]
    include=[i for i,x in enumerate(mappings) if x["decision"]=="include"];exclude=[i for i,x in enumerate(mappings) if x["decision"]!="include"]
    retained=matrix[include,:];lost=matrix[exclude,:]
    if not np.array_equal(matrix.sum(axis=0),retained.sum(axis=0)+lost.sum(axis=0)):raise ValueError("Mapping count mass failed")
    source_rows=[{"gene_id":g,**{d:str(int(matrix[i,j])) for j,d in enumerate(donors)}} for i,g in enumerate(source_ids)]
    stable_rows=[{"gene_id":mappings[i]["stable_gene_id"],**{d:str(int(matrix[i,j])) for j,d in enumerate(donors)}} for i in include]
    stable_rows.sort(key=lambda x:x["gene_id"])
    if len(stable_rows)!=len({r["gene_id"] for r in stable_rows}):raise ValueError("Stable gene mapping is nonunique")
    write_tsv(OUT/"source_counts.tsv.gz",source_rows,["gene_id"]+donors)
    write_tsv(OUT/"counts_stable.tsv.gz",stable_rows,["gene_id"]+donors)
    write_tsv(ROOT/"manifests/heca_features_mapping.tsv",mappings,list(mappings[0]))
    write_tsv(OUT/"source_features.tsv",[{"gene_id":x[0],"source_symbol":x[1],"feature_type":x[2]} for x in reference],["gene_id","source_symbol","feature_type"])
    for j,r in enumerate(ledgers):
        r["mapping_retained_UMIs"]=int(retained[:,j].sum());r["mapping_excluded_UMIs"]=int(lost[:,j].sum());r["mapping_excluded_fraction"]=float(lost[:,j].sum()/matrix[:,j].sum())
        mine=[x for x in cell_comparison if x["donor_id"]==r["donor_id"]]
        r["author_H5_original_UMIs"]=sum(x["author_H5_original_n_counts"] for x in mine)
        r["raw_minus_author_H5_UMIs"]=sum(x["raw_minus_author_H5"] for x in mine)
    write_tsv(OUT/"selected_donor_count_ledger.tsv",ledgers,list(ledgers[0]))
    write_tsv(OUT/"raw_vs_author_H5_cell_totals.tsv.gz",cell_comparison,list(cell_comparison[0]))
    differences=np.asarray([x["raw_minus_author_H5"] for x in cell_comparison])
    difference_summary={"cells":len(cell_comparison),"exact_total_matches":int(np.sum(differences==0)),"raw_greater":int(np.sum(differences>0)),"raw_less":int(np.sum(differences<0)),"min":int(differences.min()),"median":float(np.median(differences)),"max":int(differences.max()),"total_difference":int(differences.sum()),"author_H5_original_UMIs":sum(x["author_H5_original_n_counts"] for x in cell_comparison),"interpretation":"Different source feature filtering; no forced adjustment or inversion"}
    discovery={r["stable_gene_id"] for r in tsv(ROOT/"manifests/discovery_features_mapping.tsv") if r["decision"]=="include"}
    common=sorted(discovery & {r["stable_gene_id"] for r in mappings if r["decision"]=="include"})
    write_tsv(OUT/"shared_stable_features.tsv",[{"gene_id":g} for g in common],["gene_id"])
    inputs=[ROOT/"manifests/final_heca_clinical_selected.tsv",ROOT/"work/final_clinical/heca_primary_h5_obs_indices.npy",ROOT/"data/metadata/heca_h5/obs_selected.tsv",ROOT/"reports/heca_counts_raw_inventory.json",ROOT/"config/gene_mapping_reference.json",ROOT/"config/gene_mapping_policy.json",ROOT/"scripts/prepare_heca_counts.py",ROOT/"scripts/prepare_discovery_counts.py",ROOT/"data/metadata/heca_h5/heca_count_probe/original_obs_n_counts.npy"]
    outputs=[OUT/"source_counts.tsv.gz",OUT/"counts_stable.tsv.gz",OUT/"selected_donor_count_ledger.tsv",OUT/"source_features.tsv",ROOT/"manifests/heca_features_mapping.tsv",OUT/"shared_stable_features.tsv",OUT/"raw_vs_author_H5_cell_totals.tsv.gz"]
    rec={"status":"RAW_COUNTS_PREPARED_NO_DE_M0_STILL_OPEN","source":"22 exact public E-MTAB-14039 integer CellRanger processed matrix archives","selected_donors":donors,"donor_groups":dict(collections.Counter(r["disease"] for r in clinical)),"selected_primary_nuclei":sum(len(v) for v in by.values()),"original_source_gene_rows":len(source_ids),"mapped_stable_rows":len(include),"mapping_decisions":dict(collections.Counter(x["decision"] for x in mappings)),"common_with_discovery_stable_features":len(common),"common_features_are_not_frozen_expression_or_graph_universe":True,"source_integer_counts":True,"finite_nonnegative_counts":True,"source_count_order_and_coordinate_uniqueness":True,"exact_integer_mass_conservation":True,"normalized":False,"inverse_reconstruction_used":False,"rounded":False,"all_author_selected_barcodes_matched":True,"original_H5_log_object_used_for_counts":False,"archive_validation":"source byte size/SHA256, outer gzip CRC, inner gzip CRC and tar/member/sample identity verified","source_features_all_archives_identical":True,"hgnc_annotation_sha256":annotation["sha256"],"hgnc_scope_limitation":"approved named-gene universe only; unmapped/ambiguous/duplicate-target rows excluded","raw_vs_author_H5_cell_totals":difference_summary,"execution_script_sha256_at_start":EXECUTION_SCRIPT_SHA,"source_filter_difference":"raw archives retain 36601 CellRanger features; H5var33578 after author gene filtering; source totals not forced equal to H5obs/n_counts","disease_effects_computed":False,"runtime":{"host":"SRV01","project":str(ROOT),"python":platform.python_version(),"numpy":np.__version__,"requests":requests.__version__,"peak_RSS_KiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"wall_seconds":round(time.monotonic()-started,3),"disk_free_bytes_after":guard(),"execution_model_request":"gpt-6.1-sol high","model_evidence":"parent model_dispatches.json delegation configuration, not backend attestation"},"input_hashes":{str(p.relative_to(ROOT)):sha(p) for p in inputs},"archive_audits":recorded,"outputs":[{"path":str(p.relative_to(ROOT)),"bytes":p.stat().st_size,"sha256":sha(p)} for p in outputs]}
    if sha(pathlib.Path(__file__))!=EXECUTION_SCRIPT_SHA:raise ValueError("On-disk script changed during execution; final receipt withheld")
    (OUT/"provenance.json").write_text(json.dumps(rec,indent=2)+"\n");(ROOT/"receipts/heca_counts_preparation.json").write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps({k:v for k,v in rec.items() if k not in ["archive_audits","input_hashes","outputs"]},indent=2),flush=True)

def main():
    global EXECUTION_SCRIPT_SHA
    EXECUTION_SCRIPT_SHA=sha(pathlib.Path(__file__))
    ap=argparse.ArgumentParser();ap.add_argument("--download-only",action="store_true");ap.add_argument("--wait-for-downloads",action="store_true");a=ap.parse_args()
    if a.download_only:download()
    else:aggregate(a.wait_for_downloads)
if __name__=="__main__":main()

