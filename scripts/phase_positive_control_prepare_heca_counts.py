#!/usr/bin/env python3
"""Phase positive control: public raw-UMI archive intake for the secretory-arm HECA donors
and fixed donor aggregation for all 42 HECA donors (22 v1 proliferative + 20 secretory). No DE.

Reuses the v1 MatrixMarket parser, barcode join, feature check and checkpoint validation
from prepare_heca_counts.py unchanged. Writes only to a new output tree; the frozen v1
archives, checkpoints and outputs are read, never written. The 22 v1 donor columns are
re-derived from the frozen v1 archives and must equal the frozen v1 counts exactly.
Run from the project root, e.g. `python3 -I scripts/phase_positive_control_prepare_heca_counts.py --download-only`
then without the flag."""
import argparse,collections,csv,gzip,hashlib,inspect,json,pathlib,platform,resource,shutil,sys,tarfile,time
from concurrent.futures import ThreadPoolExecutor
import numpy as np,requests
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import prepare_heca_counts as v1  # pure helpers only; its module-level paths are not used
from prepare_discovery_counts import map_features,table

PRIMARY={"eStromal","eStromal_MMPs","eStromal_cycling","dStromal_early","dStromal_mid","dStromal_late","Fibroblast_basalis"}
V1_SELECTION_SHA="1a0ba5813046f9cb6e18862e17ca7fe64c31a39605d048309f1e3522447046cc"
V1_OBS_INDEX_SHA="d02ed9b8cee2c789c4f7aec7a3792a8f9444307163531933b632845cecf864ef"
V1_COUNTS_SHA="77121542de07e346ea2112f80bfee6248993c42f90c393ec3154b0a9de529af0"
DISK_RESERVE=20*1024**3;RSS_LIMIT_KIB=12*1024**2
EXECUTION_SCRIPT_SHA=None
sha=v1.sha;tsv=v1.tsv;write_tsv=v1.write_tsv

def guard(root,extra=0):
    free=shutil.disk_usage(root).free
    if free-extra<DISK_RESERVE:raise RuntimeError("20 GiB disk reserve would be violated")
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>RSS_LIMIT_KIB:raise RuntimeError("12 GiB RSS guard exceeded")
    return free
def selected_rows(a):
    rows=[r for r in tsv(a.root/a.ledger) if r["selected_primary"]=="True"]
    if len(rows)!=len({r["donor_id"] for r in rows}) or len(rows)!=len({r["sample_id"] for r in rows}):raise ValueError("Ledger selection is not one sample per donor")
    return rows
def archive_paths(a,row):
    if row["counts_available"]=="yes":p=a.root/a.v1_raw_dir/(row["sample_id"]+".tar.gz")
    else:p=a.root/a.raw_dir/row["sample_id"]/(row["sample_id"]+".tar.gz")
    return p,p.with_name(p.name+".receipt.json")

def fetch(a,row,item):
    p,receipt=archive_paths(a,row);p.parent.mkdir(parents=True,exist_ok=True);expected=int(item["bytes"])
    if p.exists():
        if not receipt.exists():raise ValueError("Existing archive lacks receipt: "+p.name)
        old=json.loads(receipt.read_text())
        if p.stat().st_size!=expected or sha(p)!=old["sha256"]:raise ValueError("Cached archive differs from its receipt: "+p.name)
        return {**old,"reused":True}
    tmp=p.with_suffix(p.suffix+".partial")
    if tmp.exists():raise ValueError("Partial archive needs explicit resume audit: "+str(tmp))
    guard(a.root,expected);started=time.monotonic()
    r=requests.get(item["url"],stream=True,timeout=(15,90));r.raise_for_status()
    content=r.headers.get("Content-Length")
    if content is not None and int(content)!=expected:r.close();raise ValueError("Content-Length differs from the download manifest: "+p.name)
    received=0;hasher=hashlib.sha256()
    try:
        with tmp.open("wb") as f:
            for b in r.iter_content(1024*1024):
                received+=len(b)
                if received>expected:raise ValueError("Archive exceeds manifest size: "+p.name)
                hasher.update(b);f.write(b)
        if received!=expected:raise ValueError("Archive shorter than manifest size: "+p.name)
        tmp.replace(p)
    finally:r.close()
    rec={"donor_id":row["donor_id"],"sample_id":row["sample_id"],"library_id":row["library_id"],"url":item["url"],"expected_bytes":expected,"bytes":received,"sha256":hasher.hexdigest(),"etag":r.headers.get("ETag"),"last_modified":r.headers.get("Last-Modified"),"final_url":r.url,"http_status":r.status_code,"wall_seconds":round(time.monotonic()-started,3),"reused":False}
    receipt.write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps({"downloaded":row["sample_id"],"bytes":received,"seconds":rec["wall_seconds"]}),flush=True)
    return rec
def download(a):
    rows=[r for r in selected_rows(a) if r["counts_available"]=="no"]
    manifest={r["sample_id"]:r for r in tsv(a.root/a.download_manifest)}
    items=[]
    for r in rows:
        m=manifest.get(r["sample_id"])
        if m is None or m["http_status"]!="200" or not m["bytes"].isdigit():raise ValueError("Download manifest lacks a sized entry for "+r["sample_id"])
        if m["donor_id"]!=r["donor_id"]:raise ValueError("Download manifest donor differs for "+r["sample_id"])
        items.append(m)
    guard(a.root,sum(int(x["bytes"]) for x in items));started=time.monotonic()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:records=list(ex.map(lambda rm:fetch(a,rm[0],rm[1]),zip(rows,items)))
    rec={"status":"PHASE_POSITIVE_CONTROL_HECA_ARCHIVES_DOWNLOADED_NOT_YET_COUNT_VALIDATED","ledger":a.ledger,"ledger_sha256":sha(a.root/a.ledger),"download_manifest_sha256":sha(a.root/a.download_manifest),"archives":records,"manifest_bytes":sum(int(x["bytes"]) for x in items),"actual_new_bytes":sum(x["bytes"] for x in records if not x["reused"]),"wall_seconds":round(time.monotonic()-started,3),"disk_free_bytes_after":guard(a.root),"versions":{"python":platform.python_version(),"numpy":np.__version__,"requests":requests.__version__},"source_scope":"public E-MTAB-14039 processed raw-count archives of the frozen secretory-arm samples; no controlled reads","expression_values_read":False}
    out=a.root/a.download_receipt;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps({k:v for k,v in rec.items() if k!="archives"},indent=2),flush=True)

def archive_metadata(p,receipt_path,sample_id):
    if not receipt_path.exists():raise ValueError("Archive receipt missing: "+p.name)
    rec=json.loads(receipt_path.read_text())
    if not p.exists() or p.stat().st_size!=int(rec["bytes"]):raise ValueError("Exact raw archive unavailable: "+p.name)
    if sha(p)!=rec["sha256"]:raise ValueError("Raw archive hash changed: "+p.name)
    with gzip.open(p,"rb") as f:  # reading all outer gzip bytes validates its CRC without extracting members
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
            if pathlib.PurePosixPath(hits[0].name).parent.name!=sample_id:raise ValueError("Archive internal sample identity differs: "+sample_id)
            lookup[tail]=hits[0].name
        with gzip.open(tar.extractfile(lookup["features.tsv.gz"]),"rt") as f:features=list(csv.reader(f,delimiter="\t"))
        with gzip.open(tar.extractfile(lookup["barcodes.tsv.gz"]),"rt") as f:barcodes=[x.strip() for x in f]
    if any(len(x)!=3 or x[2]!="Gene Expression" for x in features):raise ValueError("Unexpected source feature modality")
    return features,barcodes,lookup,rec
def selection(a,rows):
    by={r["sample_id"]:[] for r in rows};byrow={};pos_v1=[]
    v1_samples={r["sample_id"] for r in rows if r["counts_available"]=="yes"};n_obs=0
    with (a.root/a.obs).open() as f:
        for i,row in enumerate(csv.DictReader(f,delimiter="\t")):
            n_obs=i+1
            if row["sample"] in by and row["celltype"] in PRIMARY and row["discard_lowQC"]=="no":
                by[row["sample"]].append(row["_index"]);byrow[row["_index"]]=i
                if row["sample"] in v1_samples:pos_v1.append(i)
    for r in rows:
        if len(by[r["sample_id"]])!=int(r["primary_nuclei"]) or len(by[r["sample_id"]])<50:raise ValueError("Selected primary nuclei mismatch: "+r["sample_id"])
    v1_index=a.root/"work/final_clinical/heca_primary_h5_obs_indices.npy";v1_check="not_available"
    if v1_index.exists():
        if sha(v1_index)!=V1_OBS_INDEX_SHA:raise ValueError("Frozen v1 H5 index file differs from its pin")
        if not np.array_equal(np.asarray(pos_v1),np.load(v1_index)):raise ValueError("Re-derived v1 lineage positions differ from the frozen v1 H5 indices")
        v1_check="identical"
    return by,byrow,n_obs,v1_check

def aggregate(a):
    started=time.monotonic();out=a.root/a.out;out.mkdir(parents=True,exist_ok=True);check=out/"checkpoints";check.mkdir(exist_ok=True)
    v1_selection=a.root/"manifests/final_heca_clinical_selected.tsv"
    if sha(v1_selection)!=V1_SELECTION_SHA:raise ValueError("Frozen v1 clinical selection differs from its pin")
    v1_counts=a.root/"data/processed/heca/counts_stable.tsv.gz";v1_counts_sha_before=sha(v1_counts)
    if v1_counts_sha_before!=V1_COUNTS_SHA:raise ValueError("Frozen v1 count matrix differs from its pin")
    rows=selected_rows(a);rows=[r for r in rows if r["counts_available"]=="yes"]+[r for r in rows if r["counts_available"]=="no"]
    v1_order=[r["sample_id"] for r in tsv(v1_selection)]
    if [r["sample_id"] for r in rows if r["counts_available"]=="yes"]!=v1_order:raise ValueError("v1 donors must appear in the frozen v1 order")
    donors=[r["donor_id"] for r in rows];by,byrow,n_obs,v1_index_check=selection(a,rows)
    original_h5_counts=np.load(a.root/a.h5_counts)
    if original_h5_counts.shape!=(n_obs,):raise ValueError("Original H5 per-nucleus count vector length differs from the observation table")
    reference=None;columns=[];ledgers=[];recorded=[];cell_comparison=[];parser_sha=hashlib.sha256(inspect.getsource(v1.aggregate_mtx).encode()).hexdigest()
    for row in rows:
        tick=time.monotonic();s=row["sample_id"];p,rp=archive_paths(a,row)
        features,barcodes,members,download=archive_metadata(p,rp,s)
        if reference is None:reference=features
        v1.check_features(features,reference)
        mask,normalized=v1.join_barcodes(barcodes,by[s],row["library_id"])
        if int(mask.sum())!=int(row["primary_nuclei"]):raise ValueError("Selected source barcode count differs from the ledger: "+s)
        pins={"archive_sha256":download["sha256"],"selected_barcodes_sha256":hashlib.sha256(json.dumps(sorted(by[s]),separators=(",",":")).encode()).hexdigest(),"source_features_sha256":hashlib.sha256(json.dumps(features,separators=(",",":")).encode()).hexdigest(),"parser_sha256":parser_sha,"checkpoint_schema_version":2,"sample_id":s,"donor_id":row["donor_id"]}
        npz=check/(s+".npz");jsonpath=check/(s+".json")
        existing=json.loads(jsonpath.read_text()) if jsonpath.exists() else {}
        if npz.exists() and existing.get("checkpoint_schema_version")==2:
            with np.load(npz) as saved:
                v1.validate_checkpoint(existing,saved,pins,npz,len(features),mask)
                result={**existing,"gene_counts":saved["gene_counts"].copy(),"cell_counts":saved["cell_counts"].copy()}
        else:
            with tarfile.open(p,"r:gz") as tar:
                with gzip.open(tar.extractfile(members["matrix.mtx.gz"]),"rb") as matrix:result=v1.aggregate_mtx(matrix,len(features),mask)
            meta={k:v for k,v in result.items() if k not in ["gene_counts","cell_counts"]};meta.update(pins)
            np.savez_compressed(npz,gene_counts=result["gene_counts"],cell_counts=result["cell_counts"])
            meta["npz_sha256"]=sha(npz);meta["execution_script_sha256"]=EXECUTION_SCRIPT_SHA
            jsonpath.write_text(json.dumps(meta,indent=2)+"\n")
            with np.load(npz) as saved:v1.validate_checkpoint(meta,saved,pins,npz,len(features),mask)
        columns.append(result["gene_counts"])
        for index,key in enumerate(normalized):
            if mask[index]:
                original=int(original_h5_counts[byrow[key]]);raw_count=int(result["cell_counts"][index])
                cell_comparison.append({"donor_id":row["donor_id"],"sample_id":s,"cell_id":key,"original_h5_row":byrow[key],"raw_36601_feature_UMIs":raw_count,"author_H5_original_n_counts":original,"raw_minus_author_H5":raw_count-original})
        ledgers.append({**row,"raw_source_cells":len(barcodes),"matched_primary_nuclei":int(mask.sum()),"original_source_gene_rows":len(features),"selected_raw_UMIs":result["selected_total"],"raw_total_all_source_cells":result["all_total"],"matrix_nnz":result["observed_nnz"],"source_archive_sha256":download["sha256"],"source_integer_counts":True,"normalized":False,"inverted_or_rounded":False,"barcode_complete_match":True})
        recorded.append({"sample_id":s,"donor_id":row["donor_id"],"phase":row["phase"],"archive_path":str(p.relative_to(a.root)),"archive_sha256":download["sha256"],"members":members,"wall_seconds":round(time.monotonic()-tick,3),"source_shape":result["source_shape"]})
        print(json.dumps({"aggregated":s,"phase":row["phase"],"selected_nuclei":int(mask.sum()),"UMIs":result["selected_total"],"seconds":recorded[-1]["wall_seconds"]}),flush=True);guard(a.root)
    matrix=np.asarray(columns,dtype=np.float64).T;source_ids=[x[0] for x in reference]
    annotation=json.loads((a.root/"config/gene_mapping_reference.json").read_text());annpath=a.root/annotation["path"]
    if sha(annpath)!=annotation["sha256"]:raise ValueError("Pinned HGNC annotation differs")
    mappings=map_features(source_ids,table(annpath))
    for m,feature in zip(mappings,reference):m["source_symbol"]=feature[1];m["source_symbol_agrees_HGNC"]=feature[1]==m["approved_symbol"]
    include=[i for i,x in enumerate(mappings) if x["decision"]=="include"];exclude=[i for i,x in enumerate(mappings) if x["decision"]!="include"]
    retained=matrix[include,:];lost=matrix[exclude,:]
    if not np.array_equal(matrix.sum(axis=0),retained.sum(axis=0)+lost.sum(axis=0)):raise ValueError("Mapping count mass failed")
    source_rows=[{"gene_id":g,**{d:str(int(matrix[i,j])) for j,d in enumerate(donors)}} for i,g in enumerate(source_ids)]
    stable_rows=[{"gene_id":mappings[i]["stable_gene_id"],**{d:str(int(matrix[i,j])) for j,d in enumerate(donors)}} for i in include]
    stable_rows.sort(key=lambda x:x["gene_id"])
    if len(stable_rows)!=len({r["gene_id"] for r in stable_rows}):raise ValueError("Stable gene mapping is nonunique")
    # Identity check against the frozen v1 matrix: same gene rows, same 22 donor columns, same values.
    with gzip.open(v1_counts,"rt") as f:
        reader=csv.DictReader(f,delimiter="\t");v1_rows=list(reader);v1_cols=[c for c in reader.fieldnames if c!="gene_id"]
    v1_donor_ids=[r["donor_id"] for r in rows if r["counts_available"]=="yes"]
    same_genes=[r["gene_id"] for r in v1_rows]==[r["gene_id"] for r in stable_rows]
    same_values=same_genes and v1_cols==v1_donor_ids and all(o[d]==n[d] for o,n in zip(v1_rows,stable_rows) for d in v1_cols)
    v1_mapping=a.root/"manifests/heca_features_mapping.tsv";same_mapping=tsv(v1_mapping)==[{k:str(v) for k,v in m.items()} for m in mappings] if v1_mapping.exists() else "not_available"
    write_tsv(out/"source_counts.tsv.gz",source_rows,["gene_id"]+donors);write_tsv(out/"counts_stable.tsv.gz",stable_rows,["gene_id"]+donors)
    write_tsv(out/"heca_features_mapping.tsv",mappings,list(mappings[0]));write_tsv(out/"source_features.tsv",[{"gene_id":x[0],"source_symbol":x[1],"feature_type":x[2]} for x in reference],["gene_id","source_symbol","feature_type"])
    for j,r in enumerate(ledgers):
        r["mapping_retained_UMIs"]=int(retained[:,j].sum());r["mapping_excluded_UMIs"]=int(lost[:,j].sum());r["mapping_excluded_fraction"]=float(lost[:,j].sum()/matrix[:,j].sum())
        mine=[x for x in cell_comparison if x["donor_id"]==r["donor_id"]]
        r["author_H5_original_UMIs"]=sum(x["author_H5_original_n_counts"] for x in mine);r["raw_minus_author_H5_UMIs"]=sum(x["raw_minus_author_H5"] for x in mine)
    write_tsv(out/"selected_donor_count_ledger.tsv",ledgers,list(ledgers[0]));write_tsv(out/"raw_vs_author_H5_cell_totals.tsv.gz",cell_comparison,list(cell_comparison[0]))
    differences=np.asarray([x["raw_minus_author_H5"] for x in cell_comparison])
    difference_summary={"cells":len(cell_comparison),"exact_total_matches":int(np.sum(differences==0)),"raw_greater":int(np.sum(differences>0)),"raw_less":int(np.sum(differences<0)),"min":int(differences.min()),"median":float(np.median(differences)),"max":int(differences.max()),"total_difference":int(differences.sum()),"author_H5_original_UMIs":sum(x["author_H5_original_n_counts"] for x in cell_comparison),"interpretation":"Different source feature filtering; no forced adjustment or inversion"}
    discovery={r["stable_gene_id"] for r in tsv(a.root/"manifests/discovery_features_mapping.tsv") if r["decision"]=="include"}
    common=sorted(discovery & {r["stable_gene_id"] for r in mappings if r["decision"]=="include"});write_tsv(out/"shared_stable_features.tsv",[{"gene_id":g} for g in common],["gene_id"])
    inputs=[a.root/a.ledger,a.root/a.obs,a.root/a.h5_counts,v1_selection,v1_counts,a.root/"config/gene_mapping_reference.json",a.root/"config/gene_mapping_policy.json",a.root/"config/phase_positive_control.json",HERE/"prepare_heca_counts.py",HERE/"prepare_discovery_counts.py",pathlib.Path(__file__)]
    outputs=[out/"source_counts.tsv.gz",out/"counts_stable.tsv.gz",out/"selected_donor_count_ledger.tsv",out/"source_features.tsv",out/"heca_features_mapping.tsv",out/"shared_stable_features.tsv",out/"raw_vs_author_H5_cell_totals.tsv.gz"]
    phases=collections.Counter(r["phase"] for r in rows);groups=collections.Counter((r["phase"],r["disease"]) for r in rows)
    rec={"status":"PHASE_POSITIVE_CONTROL_HECA_COUNTS_PREPARED_NO_DE" if same_values else "FAILED_V1_COLUMN_IDENTITY","protocol":"docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md","source":"public E-MTAB-14039 integer CellRanger processed matrix archives; 22 frozen v1 archives read in place plus the secretory-arm archives","selected_donors":donors,"phase_counts":dict(phases),"phase_disease_counts":{f"{p}|{d}":n for (p,d),n in groups.items()},"selected_primary_nuclei":sum(len(v) for v in by.values()),"original_source_gene_rows":len(source_ids),"mapped_stable_rows":len(include),"mapping_decisions":dict(collections.Counter(x["decision"] for x in mappings)),"common_with_discovery_stable_features":len(common),"v1_identity":{"v1_lineage_h5_indices":v1_index_check,"v1_gene_rows_identical":same_genes,"v1_donor_columns_identical":same_values,"v1_feature_mapping_identical":same_mapping,"v1_counts_sha256_before":v1_counts_sha_before,"v1_counts_sha256_after":sha(v1_counts)},"source_integer_counts":True,"finite_nonnegative_counts":True,"exact_integer_mass_conservation":True,"normalized":False,"rounded":False,"all_author_selected_barcodes_matched":True,"archive_validation":"receipt byte size/SHA256, outer gzip CRC, inner gzip CRC and tar/member/sample identity verified","source_features_all_archives_identical":True,"hgnc_annotation_sha256":annotation["sha256"],"raw_vs_author_H5_cell_totals":difference_summary,"disease_or_phase_effects_computed":False,"v1_artifacts_written":False,"runtime":{"host":"SRV01","python":platform.python_version(),"numpy":np.__version__,"peak_RSS_KiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"wall_seconds":round(time.monotonic()-started,3),"disk_free_bytes_after":guard(a.root)},"execution_script_sha256_at_start":EXECUTION_SCRIPT_SHA,"input_hashes":{str(p.relative_to(a.root)):sha(p) for p in inputs},"archive_audits":recorded,"outputs":[{"path":str(p.relative_to(a.root)),"bytes":p.stat().st_size,"sha256":sha(p)} for p in outputs]}
    if sha(pathlib.Path(__file__))!=EXECUTION_SCRIPT_SHA:raise ValueError("On-disk script changed during execution; final receipt withheld")
    (out/"provenance.json").write_text(json.dumps(rec,indent=2)+"\n");r=a.root/a.receipt;r.parent.mkdir(parents=True,exist_ok=True);r.write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps({k:v for k,v in rec.items() if k not in ["archive_audits","input_hashes","outputs"]},indent=2),flush=True)
    if not same_values:raise ValueError("Re-derived v1 donor columns differ from the frozen v1 matrix")

def main():
    global EXECUTION_SCRIPT_SHA;EXECUTION_SCRIPT_SHA=sha(pathlib.Path(__file__))
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root",type=pathlib.Path,default=HERE.parent);ap.add_argument("--ledger",default="manifests/phase_positive_control_heca_candidates.tsv")
    ap.add_argument("--download-manifest",default="manifests/phase_positive_control_heca_download_manifest.tsv")
    ap.add_argument("--raw-dir",default="data/raw/phase_positive_control/heca_archives");ap.add_argument("--v1-raw-dir",default="data/metadata/heca_h5/heca_count_raw_archives")
    ap.add_argument("--obs",default="data/metadata/heca_h5/obs_selected.tsv");ap.add_argument("--h5-counts",default="data/metadata/heca_h5/heca_count_probe/original_obs_n_counts.npy")
    ap.add_argument("--out",default="data/processed/phase_positive_control/heca");ap.add_argument("--receipt",default="receipts/phase_positive_control_heca_counts.json")
    ap.add_argument("--download-receipt",default="receipts/phase_positive_control_heca_downloads.json");ap.add_argument("--workers",type=int,default=3)
    ap.add_argument("--download-only",action="store_true");a=ap.parse_args();a.root=a.root.resolve()
    if a.download_only:download(a)
    else:aggregate(a)
if __name__=="__main__":main()
