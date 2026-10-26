#!/usr/bin/env python3
"""Phase positive control: donor aggregation of the discovery atlas from the streamed GEO object.

Consumes the output of phase_positive_control_rds_stream.py (structure.pickle plus the
retained counts/meta.data payloads), validates the decoded structure against the author
route's recorded dimensions, aggregates the frozen stromal-fibroblast lineage per donor for
the 35 protocol donors in one sequential pass over the CSC payload, applies the v1 feature
rules (irreversible multi-source groups excluded, pinned HGNC mapping), and checks that the
22 v1 donor columns reproduce the frozen v1 discovery matrix. No DE. Writes only to a new
output tree. NumPy and the standard library only."""
import argparse,collections,csv,gzip,hashlib,io,json,pathlib,pickle,platform,resource,sys,time
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import phase_positive_control_rds_stream as rs
from prepare_discovery_counts import map_features,table

AUTHOR_RECORD={"features":34727,"cells":466371,"samples":60,"nonzero":1188298946,"fractional":104078188}
V1_COUNTS_SHA="c286af7cbe4cf4358afdd985fbeb42cbb4583754cc3650f12ddb92c987d8e067"
V1_LEDGER_SHA="24ba9af01c457ff6840796097a94f7f75b22b379d65b37dd0f1c8af15277c08d"
RSS_LIMIT_KIB=12*1024**2;CHUNK_COLUMNS=4000
EXECUTION_SCRIPT_SHA=None

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def tsv(path,delimiter="\t"):
    with open(path,newline="") as f:return list(csv.DictReader(f,delimiter=delimiter))
def write_tsv(path,rows,fields):
    path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if str(path).endswith(".gz"):
        with path.open("wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz, io.TextIOWrapper(gz,encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
    else:
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def guard():
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>RSS_LIMIT_KIB:raise RuntimeError("12 GiB RSS guard exceeded")
class U(pickle.Unpickler):
    def find_class(self,module,name):
        if name in ("Node","DiskVector") and module in ("__main__","phase_positive_control_rds_stream","rds_stream"):return getattr(rs,name)
        return super().find_class(module,name)
def load(p):
    with open(p,"rb") as f:return U(f).load()
def array(n):
    v=rs.value(n)
    if isinstance(v,rs.DiskVector):return np.memmap(v.path,dtype=v.dtype,mode="r",shape=(v.length,))
    if isinstance(n,rs.Node) and n.type==238:
        info,state=n.value;kind=rs.value(info[0][1]) if isinstance(info,list) and info else None
        if kind in ("compact_intseq","compact_realseq"):
            st=array(state);return np.arange(st[0],dtype=np.int64)*st[2]+st[1]
        if kind=="deferred_string":return array(state[0][1]).astype(str)
        if kind in ("wrap_integer","wrap_real","wrap_logical","wrap_string"):return array(state[0][1])
        raise ValueError(f"Unsupported ALTREP {kind}")
    if isinstance(v,dict) and "discarded_length" in v:raise ValueError("Required value was discarded")
    return np.asarray(v)
def decoded(n):
    a=rs.attrs(n);v=array(n)
    if "levels" in a:
        levels=array(a["levels"]);return np.array([levels[j-1] if 1<=j<=len(levels) else None for j in v],dtype=object)
    if isinstance(n,rs.Node) and n.type==10:return np.array([None if j==-2147483648 else bool(j) for j in v],dtype=object)
    return v
def strings(a):return [None if v is None else str(v) for v in a]

def main():
    global EXECUTION_SCRIPT_SHA;EXECUTION_SCRIPT_SHA=sha(pathlib.Path(__file__));started=time.monotonic()
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root",type=pathlib.Path,default=HERE.parent);ap.add_argument("--extract",required=True,type=pathlib.Path)
    ap.add_argument("--rds",type=pathlib.Path,default=None);ap.add_argument("--download-receipt",type=pathlib.Path,default=None)
    ap.add_argument("--ledger",default="manifests/phase_positive_control_discovery_candidates.tsv")
    ap.add_argument("--osf-dir",default="data/input_routes/osf/selected/04_Analysis/S3_COMPLETION_EXTENSION/GSE266265_DONOR_RECONSTRUCTION")
    ap.add_argument("--out",default="data/processed/phase_positive_control/discovery");ap.add_argument("--receipt",default="receipts/phase_positive_control_discovery_counts.json")
    a=ap.parse_args();root=a.root.resolve();out=root/a.out;out.mkdir(parents=True,exist_ok=True);osf=root/a.osf_dir
    lineage=json.loads((root/"config/discovery_lineage.json").read_text());included=set(lineage["include_exact_refined_labels"]);mural=set(lineage["exclude_exact_mural_labels"])
    ledger=tsv(root/a.ledger);donors=[r["donor_id"] for r in ledger]
    if len(donors)!=len(set(donors)):raise ValueError("Duplicate donors in ledger")
    v1_ledger_path=root/"manifests/final_discovery_clinical_selected.tsv"
    if sha(v1_ledger_path)!=V1_LEDGER_SHA:raise ValueError("Frozen v1 discovery ledger differs from its pin")
    v1_rows={r["donor_id"]:r for r in tsv(v1_ledger_path)};v1_donors=[d for d in donors if d in v1_rows]
    if [r["donor_id"] for r in ledger if r["counts_available"]=="yes"]!=v1_donors or len(v1_donors)!=22:raise ValueError("v1 donor set in ledger differs from the frozen ledger")
    v1_counts=root/"data/processed/discovery/counts_stable.tsv.gz";v1_counts_sha_before=sha(v1_counts)
    if v1_counts_sha_before!=V1_COUNTS_SHA:raise ValueError("Frozen v1 discovery matrix differs from its pin")
    # ---- structure -----------------------------------------------------------------------------
    report=json.loads((a.extract/"parse_report.json").read_text())
    if report["status"]!="PARSED":raise ValueError("Parser report is not PARSED")
    for q in report["payloads"]:
        if q["retained"] and sha(q["file"])!=q["sha256"]:raise ValueError("Retained payload hash changed: "+q["path"])
    rootnode=load(a.extract/"structure.pickle");ad=rs.attrs(rootnode)
    assay_names=strings(array(rs.attrs(ad["assays"])["names"]).tolist())
    if assay_names[0]!="RNA":raise ValueError(f"First assay is {assay_names[0]}, not RNA")
    c=rs.attrs(rs.attrs(rs.value(ad["assays"])[0])["counts"])
    if "dgCMatrix" not in strings(array(c["class"]).tolist()):raise ValueError("Counts are not a dgCMatrix")
    dim=array(c["Dim"]).astype(int).tolist();dn=rs.value(c["Dimnames"]);genes=strings(decoded(dn[0]));barcodes=strings(decoded(dn[1]))
    p=np.asarray(array(c["p"]),dtype=np.int64);i_dv=rs.value(c["i"]);x_dv=rs.value(c["x"])
    if not (isinstance(i_dv,rs.DiskVector) and isinstance(x_dv,rs.DiskVector)):raise ValueError("Expected on-disk i/x payloads")
    n_i=int(i_dv.length);n_x=int(x_dv.length)
    def read_slice(dv,lo,hi):  # explicit offset reads keep file pages out of the process RSS, unlike memmap slices
        return np.fromfile(dv.path,dtype=dv.dtype,count=hi-lo,offset=lo*np.dtype(dv.dtype).itemsize)
    ng,nc=dim
    if len(genes)!=ng or len(barcodes)!=nc or len(p)!=nc+1 or p[0]!=0 or p[-1]!=n_i or n_i==0 or n_i!=n_x or np.any(np.diff(p)<0):raise ValueError("CSC geometry invalid")
    if len(set(genes))!=ng or len(set(barcodes))!=nc:raise ValueError("Duplicate gene or barcode identifiers")
    md=ad["meta.data"];mnames=strings(array(rs.attrs(md)["names"]).tolist());mcols=rs.value(md);rn=strings(decoded(rs.attrs(md)["row.names"]))
    if rn!=barcodes:raise ValueError("Metadata row names differ from count barcodes")
    def col(name):
        if name not in mnames:return None
        return decoded(mcols[mnames.index(name)])
    sample=strings(col("sample"));refined=strings(col("AnnotationRefined"));main_=strings(col("AnnotationMain"))
    ncount=col("nCount_RNA");status=col("EndometriosisStatus");phase=col("MenstrualCyclePhase");batch=col("X10x.capturing.batch") if "X10x.capturing.batch" in mnames else col("batch")
    if sample is None or refined is None or main_ is None or ncount is None:raise ValueError("Required metadata columns missing")
    ncount=np.asarray(ncount,dtype=np.float64)
    structure={"assays":assay_names,"dim":dim,"nnz":n_x,"samples_in_metadata":len(set(sample)),"metadata_columns":mnames,"matches_author_record":{"features":ng==AUTHOR_RECORD["features"],"cells":nc==AUTHOR_RECORD["cells"],"samples":len(set(sample))==AUTHOR_RECORD["samples"],"nonzero":n_x==AUTHOR_RECORD["nonzero"]}}
    if not all(structure["matches_author_record"].values()):raise ValueError("Decoded structure differs from the author record: "+json.dumps(structure))
    # ---- donor and lineage selection -------------------------------------------------------------
    sample_arr=np.array(sample,dtype=object);donor_index={d:k for k,d in enumerate(donors)}
    # Lineage = source main label "stromal" AND refined label in the frozen 12-label list (config/discovery_lineage.json).
    # Cells whose refined label is in the list but whose main label is missing or different are excluded and counted,
    # never relabeled. v1 recorded no such conflict among its 22 donors, so their membership must be unchanged.
    col_donor=np.full(nc,-1,dtype=np.int64);member=np.zeros(nc,dtype=bool);conflict=collections.Counter();conflict_main=collections.Counter()
    for j in range(nc):
        d=donor_index.get(sample_arr[j]);
        if d is None:continue
        lab=refined[j]
        if lab is not None and lab in included:
            if lab in mural:raise ValueError("Mural label enters lineage")
            if main_[j]!="stromal":conflict[donors[d]]+=1;conflict_main[str(main_[j])]+=1;continue
            col_donor[j]=d;member[j]=True
    all_cells=collections.Counter(sample_arr[j] for j in range(nc) if sample_arr[j] in donor_index)
    cells=np.bincount(col_donor[member],minlength=len(donors))
    author_manifest={r["sample"]:int(r["cells"]) for r in tsv(osf/"DONOR_CELLTYPE_MANIFEST.csv",",") if r["population_id"]=="P077"}
    cross=[]
    for k,r in enumerate(ledger):
        d=r["donor_id"];idx=np.flatnonzero(sample_arr==d)
        nn=lambda arr:[v for v in strings(arr[idx]) if v is not None] if arr is not None else []
        st=collections.Counter(nn(status));ph=collections.Counter(nn(phase));bt=collections.Counter(nn(batch))
        disease_meta={"Endo":"ENDO","CTL":"non-ENDO"}.get(next(iter(st),None),next(iter(st),None)) if st else None
        meta_phase=next(iter(ph),None) if ph else None;meta_batch=next(iter(bt),None) if bt else None
        low=str(meta_phase).lower() if meta_phase is not None else ""
        coarse_meta=("secretory" if "secretory" in low else "proliferative" if "prolif" in low and "peri" not in low else low) if meta_phase is not None else None
        if r["counts_available"]=="yes" and conflict.get(d,0):raise ValueError("v1 donor has main/refined label conflicts, which v1 recorded as absent: "+d)
        row={"donor_id":d,"phase":r["phase"],"sub_phase":r["sub_phase"],"disease":r["disease"],"cells_all_labels":int(all_cells.get(d,0)),"primary_lineage_cells":int(cells[k]),"refined_in_lineage_but_main_not_stromal_excluded":int(conflict.get(d,0)),"ledger_primary_lineage_cells":r["primary_lineage_cells"],"author_P077_cells":author_manifest.get(d),"metadata_disease":disease_meta,"metadata_phase":meta_phase,"metadata_capture_batch":meta_batch,"ledger_capture_batch":r["capture_batch"],"single_status_phase_batch":len(st)<=1 and len(ph)<=1 and len(bt)<=1,"disease_agrees":disease_meta==r["disease"] if disease_meta is not None else None,"coarse_phase_agrees":coarse_meta==r["phase"] if coarse_meta is not None else None,"sub_phase_agrees":str(meta_phase).lower()==r["sub_phase"].lower() if meta_phase is not None else None,"batch_agrees":str(meta_batch).replace("batch","")==r["capture_batch"].replace("batch","") if meta_batch is not None else None,"meets_50_cells":int(cells[k])>=50}
        if r["counts_available"]=="yes":
            if int(cells[k])!=int(r["primary_lineage_cells"]) or int(cells[k])!=author_manifest.get(d,-1) or int(cells[k])!=int(v1_rows[d]["primary_lineage_cells"]):raise ValueError("v1 donor lineage cell count differs from the frozen ledger: "+d)
        if row["disease_agrees"] is False or row["coarse_phase_agrees"] is False or not row["single_status_phase_batch"]:raise ValueError("Metadata disagrees with the frozen ledger: "+json.dumps(row,default=str))
        cross.append(row)
    eligible_secretory=[r for r in cross if r["phase"]=="secretory" and r["meets_50_cells"]]
    # ---- one sequential pass: validation, per-cell totals, per-donor sums ----------------------------
    sums=np.zeros((len(donors),ng),dtype=np.float64);cell_tot=np.zeros(nc,dtype=np.float64);negative=nonfinite=fractional=0;total=0.0;maxv=0.0;entries=0
    for c0 in range(0,nc,CHUNK_COLUMNS):
        c1=min(nc,c0+CHUNK_COLUMNS);lo=int(p[c0]);hi=int(p[c1])
        if hi==lo:continue
        ii=np.asarray(read_slice(i_dv,lo,hi),dtype=np.int64);xv=np.asarray(read_slice(x_dv,lo,hi),dtype=np.float64);lens=np.diff(p[c0:c1+1])
        if len(ii)!=hi-lo or len(xv)!=hi-lo:raise ValueError("Short read from payload file")
        if np.any(ii<0) or np.any(ii>=ng):raise ValueError("Row index out of range")
        colid=np.repeat(np.arange(c0,c1,dtype=np.int64),lens)
        negative+=int((xv<0).sum());nonfinite+=int((~np.isfinite(xv)).sum());fractional+=int((xv!=np.floor(xv)).sum());total+=float(xv.sum());maxv=max(maxv,float(xv.max(initial=0)));entries+=len(xv)
        cell_tot[c0:c1]=np.bincount(colid-c0,weights=xv,minlength=c1-c0)
        dsel=col_donor[colid];use=dsel>=0
        if use.any():sums+=np.bincount(dsel[use]*ng+ii[use],weights=xv[use],minlength=len(donors)*ng).reshape(len(donors),ng)
        if (c0//CHUNK_COLUMNS)%20==0:guard();print(json.dumps({"columns_done":c1,"elapsed_s":round(time.monotonic()-started)}),flush=True)
    if entries!=n_x:raise ValueError("Entry count differs from nnz")
    if negative or nonfinite:raise ValueError("Counts fail the finite non-negative requirement")
    if fractional!=AUTHOR_RECORD["fractional"]:raise ValueError(f"Fractional entry count {fractional} differs from the author record")
    # conservation: per-cell totals against metadata nCount_RNA, all cells and selected cells
    diff=cell_tot-ncount;ok=np.isclose(cell_tot,ncount,rtol=1e-12,atol=1e-6)
    conservation={"cells":nc,"cells_matching_nCount_RNA":int(ok.sum()),"selected_cells":int(member.sum()),"selected_cells_matching_nCount_RNA":int(ok[member].sum()),"max_abs_difference_all":float(np.abs(diff).max()),"max_abs_difference_selected":float(np.abs(diff[member]).max()) if member.any() else 0.0}
    if conservation["selected_cells_matching_nCount_RNA"]!=conservation["selected_cells"]:raise ValueError("Selected-cell totals do not match metadata nCount_RNA")
    donor_total=sums.sum(axis=1);donor_cell_total=np.bincount(col_donor[member],weights=cell_tot[member],minlength=len(donors))
    if not np.allclose(donor_total,donor_cell_total,rtol=1e-12,atol=1e-6):raise ValueError("Per-donor gene sums differ from per-donor cell totals")
    # ---- v1 feature rules: irreversible groups excluded, pinned HGNC mapping -----------------------
    mapping_rows=table(osf/"SOURCE_GENE_MAPPING.csv",",")
    if [r["source_gene"] for r in mapping_rows]!=genes:raise ValueError("Source gene order differs from the author mapping table")
    groups=collections.defaultdict(list)
    for r in mapping_rows:groups[r["canonical_gene"]].append(r["source_gene"])
    many={g:v for g,v in groups.items() if len(v)>1};excluded_features={s for v in many.values() for s in v}
    singleton_idx=[k for k,g in enumerate(genes) if g not in excluded_features]
    if len(many)!=50 or len(singleton_idx)!=34550:raise ValueError(f"Expected 50 irreversible groups / 34550 singletons, found {len(many)} / {len(singleton_idx)}")
    annotation=json.loads((root/"config/gene_mapping_reference.json").read_text());annpath=root/annotation["path"]
    if sha(annpath)!=annotation["sha256"]:raise ValueError("Pinned HGNC annotation differs")
    mappings=map_features([genes[k] for k in singleton_idx],table(annpath))
    v1_map={r["source_feature_id"]:r for r in tsv(root/"manifests/discovery_features_mapping.tsv")}
    same_mapping=len(v1_map)==len(mappings) and all(v1_map.get(m["source_feature_id"],{}).get("decision")==m["decision"] and v1_map.get(m["source_feature_id"],{}).get("stable_gene_id")==m["stable_gene_id"] for m in mappings)
    include=[(singleton_idx[j],m["stable_gene_id"]) for j,m in enumerate(mappings) if m["decision"]=="include"]
    include.sort(key=lambda t:t[1])
    if len({g for _,g in include})!=len(include):raise ValueError("Stable gene mapping is nonunique")
    stable=np.stack([sums[:,k] for k,_ in include],axis=1) if include else np.zeros((len(donors),0))
    stable_ids=[g for _,g in include]
    singleton_total=sums[:,singleton_idx].sum(axis=1);mapped_total=stable.sum(axis=1);many_total=donor_total-singleton_total
    # ---- identity check against the frozen v1 matrix -------------------------------------------------
    with gzip.open(v1_counts,"rt") as f:
        reader=csv.DictReader(f,delimiter="\t");v1_tab=list(reader);v1_cols=[cc for cc in reader.fieldnames if cc!="gene_id"]
    same_genes=[r["gene_id"] for r in v1_tab]==stable_ids and v1_cols==v1_donors
    max_rel=None;same_values=False
    if same_genes:
        old=np.array([[float(r[d]) for d in v1_cols] for r in v1_tab]).T;new=np.stack([stable[donor_index[d],:] for d in v1_donors])
        rel=np.abs(new-old)/np.maximum(1.0,np.abs(old));max_rel=float(rel.max());same_values=bool(max_rel<=1e-6)
    v1_led=tsv(root/"data/processed/discovery/selected_donor_count_ledger.tsv") if (root/"data/processed/discovery/selected_donor_count_ledger.tsv").exists() else []
    v1_totals={r["donor_id"]:float(r["original_P077_estimated_UMIs"]) for r in v1_led} if v1_led else {}
    v1_total_agree=all(np.isclose(donor_total[donor_index[d]],t,rtol=1e-9,atol=1e-3) for d,t in v1_totals.items()) if v1_totals else None
    # ---- outputs --------------------------------------------------------------------------------------
    fmt=lambda v:repr(float(v))
    write_tsv(out/"source_counts_all_features.tsv.gz",[{"gene_id":g,**{d:fmt(sums[k,gi]) for k,d in enumerate(donors)}} for gi,g in enumerate(genes)],["gene_id"]+donors)
    write_tsv(out/"counts_stable.tsv.gz",[{"gene_id":g,**{d:fmt(stable[k,j]) for k,d in enumerate(donors)}} for j,g in enumerate(stable_ids)],["gene_id"]+donors)
    write_tsv(out/"discovery_features_mapping.tsv",mappings,list(mappings[0]))
    write_tsv(out/"irreversible_group_exclusions.tsv",[{"canonical_gene":g,"original_source_features":"|".join(v),"source_row_count":len(v)} for g,v in sorted(many.items())],["canonical_gene","original_source_features","source_row_count"])
    led=[]
    for k,r in enumerate(ledger):
        led.append({**r,"primary_lineage_cells":int(cells[k]),"cells_all_labels":cross[k]["cells_all_labels"],"meets_50_cells":cross[k]["meets_50_cells"],"metadata_phase":cross[k]["metadata_phase"],"sub_phase_agrees":cross[k]["sub_phase_agrees"],"metadata_capture_batch":cross[k]["metadata_capture_batch"],"batch_agrees":cross[k]["batch_agrees"],"all_feature_UMIs":fmt(donor_total[k]),"irreversible_group_excluded_UMIs":fmt(many_total[k]),"singleton_UMIs":fmt(singleton_total[k]),"mapped_retained_UMIs":fmt(mapped_total[k]),"mapping_excluded_UMIs":fmt(singleton_total[k]-mapped_total[k]),"overall_excluded_fraction":fmt((donor_total[k]-mapped_total[k])/donor_total[k]) if donor_total[k] else "NA","nCount_RNA_cells_matching":int(ok[member&(col_donor==k)].sum()),"counts_source":"GEO RDS stream; fractional estimates unrounded"})
    write_tsv(out/"selected_donor_count_ledger.tsv",led,list(led[0]))
    selcells=np.flatnonzero(member)
    write_tsv(out/"selected_cell_totals_vs_nCount_RNA.tsv.gz",[{"donor_id":donors[col_donor[j]],"barcode":barcodes[j],"refined_label":refined[j],"raw_total":fmt(cell_tot[j]),"nCount_RNA":fmt(ncount[j]),"difference":fmt(cell_tot[j]-ncount[j])} for j in selcells],["donor_id","barcode","refined_label","raw_total","nCount_RNA","difference"])
    download=json.loads(a.download_receipt.read_text()) if a.download_receipt and a.download_receipt.exists() else None
    rds_info={"file":a.rds.name if a.rds else None,"bytes":a.rds.stat().st_size if a.rds and a.rds.exists() else None,"download_receipt":download}
    def label(pp):  # receipt keys never carry host-specific absolute paths
        try:return str(pp.resolve().relative_to(root))
        except ValueError:return "extract/"+pp.name
    inputs=[root/a.ledger,v1_ledger_path,v1_counts,root/"config/discovery_lineage.json",root/"config/gene_mapping_reference.json",root/"config/gene_mapping_policy.json",root/"config/phase_positive_control.json",root/"manifests/discovery_features_mapping.tsv",osf/"SOURCE_GENE_MAPPING.csv",osf/"DONOR_CELLTYPE_MANIFEST.csv",a.extract/"parse_report.json",a.extract/"structure.pickle",HERE/"phase_positive_control_rds_stream.py",HERE/"prepare_discovery_counts.py",pathlib.Path(__file__)]
    outputs=[out/"source_counts_all_features.tsv.gz",out/"counts_stable.tsv.gz",out/"discovery_features_mapping.tsv",out/"irreversible_group_exclusions.tsv",out/"selected_donor_count_ledger.tsv",out/"selected_cell_totals_vs_nCount_RNA.tsv.gz"]
    phase_counts=collections.Counter(r["phase"] for r in ledger);pd_counts=collections.Counter((r["phase"],r["disease"]) for r in ledger)
    status_word="PHASE_POSITIVE_CONTROL_DISCOVERY_COUNTS_PREPARED_NO_DE" if (same_values and same_mapping) else "FAILED_V1_IDENTITY"
    rec={"status":status_word,"protocol":"docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md","source":"GEO GSE266265_EndoAtlas.rds RNA counts streamed by the ported author parser; author labels and metadata inherited","rds":rds_info,"structure":structure,"parser_report":{k:report.get(k) for k in ("header","nodes","uncompressed_bytes","elapsed_s","peak_RSS_bytes","rss_unit_multiplier","rss_ceiling_bytes","disk_reserve_bytes","archived_parser_sha256")},"selected_donors":donors,"phase_counts":dict(phase_counts),"phase_disease_counts":{f"{p_}|{d_}":n for (p_,d_),n in pd_counts.items()},"primary_lineage_cells_total":int(cells.sum()),"lineage_rule":"AnnotationMain == stromal AND AnnotationRefined in the 12 frozen labels; missing or conflicting main labels excluded and counted, not inferred","lineage_main_label_conflicts_excluded":{"cells":int(sum(conflict.values())),"donors":len(conflict),"by_main_label":dict(conflict_main)},"secretory_donors_meeting_50_cells":len(eligible_secretory),"secretory_donors_total":int(phase_counts.get("secretory",0)),"donor_cross_checks":cross,"count_scan":{"nonzero":entries,"negative":negative,"nonfinite":nonfinite,"fractional":fractional,"total":total,"maximum":maxv},"conservation_vs_nCount_RNA":conservation,"features":{"source":ng,"irreversible_groups":len(many),"features_in_irreversible_groups":len(excluded_features),"singletons":len(singleton_idx),"mapped_stable":len(stable_ids),"mapping_decisions":dict(collections.Counter(m["decision"] for m in mappings)),"hgnc_annotation_sha256":annotation["sha256"]},"v1_identity":{"v1_lineage_cell_counts_identical":True,"v1_feature_mapping_identical":same_mapping,"v1_gene_rows_identical":same_genes,"v1_donor_columns_within_1e-6_relative":same_values,"v1_max_relative_difference":max_rel,"v1_all_feature_totals_agree":v1_total_agree,"v1_counts_sha256_before":v1_counts_sha_before,"v1_counts_sha256_after":sha(v1_counts)},"normalized":False,"rounded":False,"disease_or_phase_effects_computed":False,"v1_artifacts_written":False,"runtime":{"host":"SRV01","python":platform.python_version(),"numpy":np.__version__,"peak_RSS_KiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"wall_seconds":round(time.monotonic()-started,3)},"execution_script_sha256_at_start":EXECUTION_SCRIPT_SHA,"input_hashes":{label(pp):sha(pp) for pp in inputs},"outputs":[{"path":str(pp.relative_to(root)),"bytes":pp.stat().st_size,"sha256":sha(pp)} for pp in outputs]}
    if sha(pathlib.Path(__file__))!=EXECUTION_SCRIPT_SHA:raise ValueError("On-disk script changed during execution; final receipt withheld")
    (out/"provenance.json").write_text(json.dumps(rec,indent=2,default=str)+"\n");rp=root/a.receipt;rp.parent.mkdir(parents=True,exist_ok=True);rp.write_text(json.dumps(rec,indent=2,default=str)+"\n")
    print(json.dumps({k:v for k,v in rec.items() if k not in ("donor_cross_checks","input_hashes","outputs")},indent=2,default=str),flush=True)
    if status_word!="PHASE_POSITIVE_CONTROL_DISCOVERY_COUNTS_PREPARED_NO_DE":raise ValueError("v1 identity check failed; see receipt")
if __name__=="__main__":main()
