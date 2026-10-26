#!/usr/bin/env python3
"""Check the ported RDS stream parser against the R-written synthetic fixture (NumPy only).

Runs the parser twice: once normally (values, metadata decoding, retained/discarded
payload policy are compared with FIXTURE_EXPECTED.json) and once with a 1 MiB memory
ceiling, which must abort, proving the Linux resident-memory unit fix is active.
Writes FIXTURE_CHECK.json into the fixture directory and exits nonzero on failure."""
import argparse,hashlib,json,pathlib,pickle,shutil,subprocess,sys
import numpy as np
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import phase_positive_control_rds_stream as rs

class U(pickle.Unpickler):
    def find_class(self,module,name):
        if name in ("Node","DiskVector") and module in ("__main__","phase_positive_control_rds_stream","rds_stream"):return getattr(rs,name)
        return super().find_class(module,name)
def load(p):
    with open(p,"rb") as f:return U(f).load()
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def array(n):  # same decoding rules as the archived inspect_extracted.array, including ALTREP forms
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

def main():
    ap=argparse.ArgumentParser();ap.add_argument("fixture_dir");ap.add_argument("--python",default=sys.executable);ap.add_argument("--disk-reserve-gib",default="75");a=ap.parse_args()
    d=pathlib.Path(a.fixture_dir).resolve();expected=json.loads((d/"FIXTURE_EXPECTED.json").read_text());parser=HERE/"phase_positive_control_rds_stream.py"
    out=d/"extract"
    if out.exists():shutil.rmtree(out)
    run=subprocess.run([a.python,"-I",str(parser),str(d/"fixture.rds"),str(out),"--disk-reserve-gib",a.disk_reserve_gib],capture_output=True,text=True)
    report=json.loads((out/"parse_report.json").read_text());checks={"parser_exit_zero":run.returncode==0,"status_parsed":report["status"]=="PARSED"}
    if report["status"]!="PARSED":
        print(json.dumps({"status":"FAIL","parser_error":report.get("error"),"stderr_tail":run.stderr[-1500:]},indent=2));sys.exit(1)
    root=load(out/"structure.pickle");ad=rs.attrs(root)
    names=[str(s) for s in array(rs.attrs(ad["assays"])["names"]).tolist()];checks["assay_names"]=names==expected["assay_names"]
    assays=rs.value(ad["assays"]);c=rs.attrs(rs.attrs(assays[0])["counts"])
    checks["class_dgCMatrix"]="dgCMatrix" in [str(s) for s in array(c["class"]).tolist()]
    i=array(c["i"]);p=array(c["p"]);x=array(c["x"]);dim=array(c["Dim"]).astype(int).tolist();dn=rs.value(c["Dimnames"])
    genes=decoded(dn[0]);cells=decoded(dn[1])
    checks["dim"]=dim==expected["counts_dim"];checks["nnz"]=len(x)==expected["counts_nnz"]==len(i)
    checks["sum"]=bool(np.isclose(float(np.sum(np.asarray(x,dtype=np.float64))),expected["counts_sum"],rtol=1e-12,atol=0))
    cols=np.array([float(np.sum(np.asarray(x[p[j]:p[j+1]],dtype=np.float64))) for j in range(len(p)-1)])
    rows=np.bincount(np.asarray(i,dtype=np.int64),weights=np.asarray(x,dtype=np.float64),minlength=dim[0])
    checks["colSums"]=bool(np.allclose(cols,expected["colSums"],rtol=1e-12,atol=1e-9));checks["rowSums"]=bool(np.allclose(rows,expected["rowSums"],rtol=1e-12,atol=1e-9))
    checks["fractional_entries"]=int(np.sum(np.asarray(x)!=np.floor(np.asarray(x))))==expected["fractional_entries"]
    checks["genes_head"]=[str(g) for g in genes[:3]]==expected["genes_head"];checks["cells_head"]=[str(g) for g in cells[:3]]==expected["cells_head"]
    checks["i_p_x_on_disk"]=all(isinstance(rs.value(c[k]),rs.DiskVector) for k in ("i","p","x")) or isinstance(rs.value(c["x"]),rs.DiskVector)
    pl=report["payloads"];dense=[q for q in pl if q["path"].endswith("/@/data")]
    checks["dense_discarded"]=len(dense)==2 and not any(q["retained"] for q in dense) and all(q["bytes"]==expected["dense_bytes_each"] for q in dense)
    checks["retained_only_counts_and_meta"]=all(q["path"].startswith("root/@/assays/[0]/@/counts") or q["path"].startswith("root/@/meta.data") for q in pl if q["retained"])
    checks["second_assay_counts_discarded"]=any(q["path"].startswith("root/@/assays/[1]/@/counts") and not q["retained"] for q in pl)
    md=ad["meta.data"];mnames=[str(s) for s in array(rs.attrs(md)["names"]).tolist()];cols_=rs.value(md);m={k:decoded(v) for k,v in zip(mnames,cols_)};rn=decoded(rs.attrs(md)["row.names"])
    e=expected["meta"];checks["meta_columns"]=mnames==e["columns"];checks["meta_rows"]=all(len(v)==e["rows"] for v in m.values())
    checks["sample_levels"]=[str(s) for s in array(rs.attrs(cols_[0])["levels"]).tolist()]==e["sample_levels"]
    checks["sample_codes"]=array(cols_[0])[:10].tolist()==e["sample_codes_head"]
    checks["intcol_sum"]=int(np.sum(array(cols_[4])))==e["intcol_sum"]
    checks["refined_na"]=sum(v is None for v in m["AnnotationRefined"])==e["refined_na"]
    checks["refined_head"]=[None if v is None else str(v) for v in m["AnnotationRefined"][:10]]==e["refined_head"]
    checks["flag_head"]=list(m["flag"][:10])==e["flag_head"]
    checks["rownames"]=[str(s) for s in rn[:3]]==e["rownames_head"] and len(rn)==len(cells) and all(str(r)==str(cc) for r,cc in zip(rn,cells))
    checks["nCount_RNA_retained_on_disk"]=isinstance(rs.value(cols_[1]),rs.DiskVector) and bool(np.allclose(np.asarray(array(cols_[1]),dtype=np.float64),expected["colSums"],rtol=1e-12,atol=1e-9))
    checks["peak_rss_plausible"]=10*1024**2<report["peak_RSS_bytes"]<4*1024**3
    checks["rss_unit_multiplier_linux"]=(report["rss_unit_multiplier"]==1024)==sys.platform.startswith("linux")
    out2=d/"extract_tiny_ceiling"
    if out2.exists():shutil.rmtree(out2)
    run2=subprocess.run([a.python,"-I",str(parser),str(d/"fixture.rds"),str(out2),"--rss-ceiling-gib","0.001"],capture_output=True,text=True)
    rep2=json.loads((out2/"parse_report.json").read_text())
    checks["tiny_ceiling_aborts"]=run2.returncode!=0 and rep2["status"]=="FAILED" and "MemoryError" in rep2["error"]
    status="PASS" if all(checks.values()) else "FAIL"
    result={"status":status,"checks":checks,"parser_sha256":sha(parser),"archived_parser_sha256":rs.ARCHIVED_SHA256,"fixture_sha256":sha(d/"fixture.rds"),"fixture_bytes":(d/"fixture.rds").stat().st_size,"fixture_r_version":expected.get("r_version"),"parse_report":{k:report.get(k) for k in ("header","nodes","uncompressed_bytes","elapsed_s","peak_RSS_bytes","rss_unit_multiplier","rss_ceiling_bytes","disk_reserve_bytes","platform","compression")},"tiny_ceiling_error":rep2.get("error"),"synthetic_only":True,"biological_validation":False}
    (d/"FIXTURE_CHECK.json").write_text(json.dumps(result,indent=2,default=str)+"\n");print(json.dumps(result,indent=2,default=str))
    if status!="PASS":print(run.stderr[-3000:],file=sys.stderr);sys.exit(1)
if __name__=="__main__":main()
