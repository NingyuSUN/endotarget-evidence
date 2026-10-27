#!/usr/bin/env python3
"""Phase positive control F1 freeze: verify and hash the discovery SEC_minus_PRO fit, the
B0/B1/G1 rankings and the marker-panel readout before any HECA phase effect is computed,
and publish compact copies of the frozen products into the repository tree.
Standard library only. Paths in the receipt are relative to the project root."""
import argparse,csv,datetime,gzip,hashlib,json,pathlib,shutil

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def rows(p):
    with open(p,newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def gz_copy(src,dst):  # deterministic gzip (no name, mtime 0) so the published copy hashes reproducibly
    with open(src,"rb") as f, open(dst,"wb") as raw, gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz:
        shutil.copyfileobj(f,gz)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1])
    ap.add_argument("--de-dir",required=True);ap.add_argument("--rank-dir",required=True)
    ap.add_argument("--markers",required=True);ap.add_argument("--marker-summary",required=True)
    ap.add_argument("--git-head",required=True,help="commit of the published checkout holding this code")
    ap.add_argument("--receipt",default="receipts/phase_positive_control_f1_discovery_freeze.json")
    ap.add_argument("--publish-dir",default="reports/phase_positive_control_f1");a=ap.parse_args()
    root=a.root.resolve();de=root/a.de_dir;rk=root/a.rank_dir;pub=root/a.publish_dir
    rel=lambda p:str(pathlib.Path(p).resolve().relative_to(root))
    dp=read(de/"provenance.json");rp=read(rk/"provenance.json");ms=read(root/a.marker_summary)
    if dp.get("contrast")!="SEC_minus_PRO" or dp.get("cohort")!="discovery" or dp.get("data_kind")!="public_cohort" or dp.get("heca_phase_effects_computed") is not False:
        raise ValueError("Discovery fit provenance is not the phase positive-control discovery fit")
    if rp.get("status")!="DISCOVERY_RANKING_ARTIFACT_NOT_COMPLETE_EXTERNAL_FREEZE" or rp.get("external_data_used") is not False:
        raise ValueError("Ranking provenance invalid")
    for original,digest in dp["output_sha256"].items():
        if sha(de/pathlib.Path(original).name)!=digest:raise ValueError("Changed DE output: "+original)
    for name,item in rp["outputs"].items():
        if sha(rk/name)!=item["sha256"]:raise ValueError("Changed ranking output: "+name)
    if rp["inputs"]["discovery_de"]["sha256"]!=sha(de/"de_fitted.tsv") or rp["inputs"]["discovery_de_all"]["sha256"]!=sha(de/"de_all.tsv"):
        raise ValueError("Ranking did not consume this fit's tables")
    if ms.get("cohort")!="discovery":raise ValueError("Marker summary is not the discovery readout")
    # inputs the fit declared, re-hashed now
    inputs={}
    for p,digest in dp["input_sha256"].items():
        if sha(root/p)!=digest:raise ValueError("Fit input changed since the fit: "+p)
        inputs[p]=digest
    for p in ["scripts/rank_targets.py","scripts/phase_positive_control_markers.py","scripts/phase_positive_control_freeze_f1.py","config/ranking_algorithm.json","data/reference/string/direct_evidence_gene_edges.tsv.gz","docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md"]:
        inputs[p]=sha(root/p)
    if inputs["data/reference/string/direct_evidence_gene_edges.tsv.gz"]!=rp["inputs"]["reference_edges"]["sha256"]:raise ValueError("Ranking used a different network")
    # publish compact copies
    if pub.exists() and any(pub.iterdir()):raise ValueError("Publish directory is not empty: "+str(pub))
    pub.mkdir(parents=True,exist_ok=True)
    published={}
    for name in ["B0_top20.tsv","B0_top50.tsv","B1_top20.tsv","B1_top50.tsv","G1_top20.tsv","G1_top50.tsv","universe.tsv","excluded_genes.tsv","provenance.json"]:
        shutil.copyfile(rk/name,pub/("ranking_"+name));published["ranking_"+name]=rel(pub/("ranking_"+name))
    for name in ["expression_filter.tsv","normalization.tsv","design.tsv","provenance.json"]:
        shutil.copyfile(de/name,pub/("de_"+name));published["de_"+name]=rel(pub/("de_"+name))
    gz_copy(de/"de_all.tsv",pub/"de_all.tsv.gz");published["de_all.tsv.gz"]=rel(pub/"de_all.tsv.gz")
    shutil.copyfile(root/a.markers,pub/"markers_discovery.tsv");shutil.copyfile(root/a.marker_summary,pub/"markers_discovery_summary.json")
    published["markers_discovery.tsv"]=rel(pub/"markers_discovery.tsv");published["markers_discovery_summary.json"]=rel(pub/"markers_discovery_summary.json")
    products={"discovery_de_all":rel(de/"de_all.tsv"),"discovery_de_fitted":rel(de/"de_fitted.tsv"),"expression_filter":rel(de/"expression_filter.tsv"),"normalization":rel(de/"normalization.tsv"),"design":rel(de/"design.tsv"),"de_provenance":rel(de/"provenance.json"),"universe":rel(rk/"universe.tsv"),"induced_edges":rel(rk/"induced_edges.tsv"),"excluded_genes":rel(rk/"excluded_genes.tsv"),"rank_provenance":rel(rk/"provenance.json"),"markers_discovery":a.markers,"markers_discovery_summary":a.marker_summary}
    for m in ("B0","B1","G1"):
        for suffix in ("","_top20","_top50"):products[f"{m}{suffix or '_full'}"]=rel(rk/f"{m}{suffix}.tsv")
    artifacts={**inputs,**{p:sha(root/p) for p in products.values()},**{v:sha(root/v) for v in published.values()}}
    top={m:[{"rank":int(r["rank"]),"gene_id":r["gene_id"],"discovery_sign":int(r["discovery_sign"]),"t":float(r["t"]),"adj.P.Val":float(r["adj.P.Val"])} for r in rows(rk/f"{m}_top20.tsv")] for m in ("B0","B1","G1")}
    de_all=rows(de/"de_all.tsv");q=[float(r["adj.P.Val"]) for r in de_all if r["adj.P.Val"]!="NA"]
    summary={"tested_family_size":dp["tested_family_size"],"estimable_genes":dp["estimable_genes"],"BH_q_lt_0_05":dp["significant_BH05"],"BH_q_lt_0_05_up":dp.get("significant_BH05_up"),"BH_q_lt_0_05_down":dp.get("significant_BH05_down"),"BH_q_lt_0_01":sum(1 for v in q if v<0.01),"design":dp["audit"],"filter_audit":dp["filter_audit"],"captured_warnings":dp.get("captured_warnings"),"ranking_counts":rp["counts"],"propagation":rp["propagation"],"top20":top,"markers_discovery":{k:ms[k] for k in ("panel_genes","in_expression_family","expected_sign_agrees","expected_sign_agrees_and_q_lt_0_05","threshold","pass")}}
    p1=read(root/"config/phase_positive_control.json")["endpoints"]["P1"]["pass"]["discovery_min"]
    receipt={"schema_version":1,"phase":"PHASE_POSITIVE_CONTROL_F1_DISCOVERY_FROZEN_PRE_HECA","created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"git_head":a.git_head,"protocol":"docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md","f0_freeze":"receipts/phase_positive_control_freeze.json","f0_artifact_check_on_analysis_host":dp.get("freeze_artifact_check"),"artifacts":artifacts,"products":products,"published_copies":published,"summary":summary,"P1_discovery":{"genes_BH_q_lt_0_05":dp["significant_BH05"],"threshold":p1,"discovery_component_met":dp["significant_BH05"]>=p1,"note":"P1 passes only if HECA also meets its threshold; HECA is not fitted"},"heca_phase_effects_computed":False,"lists_may_change_after_this_freeze":False,"software":{"fit":dp.get("software"),"ranking":rp.get("software")},"interpretation":"frozen discovery lists and directions for a pipeline positive control; not disease evidence, not targets"}
    out=root/a.receipt;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ("phase","git_head","summary","P1_discovery")},indent=2))
if __name__=="__main__":main()
