#!/usr/bin/env python3
"""Phase positive control, endpoint S1: read the frozen seven-gene marker panel from a full
DE table (de_all.tsv) for one cohort. Direction-only check with the prespecified signs;
genes outside the expression-pass family count as failures. Standard library only."""
import argparse,csv,gzip,json,pathlib

def tsv(p):
    opener=gzip.open if str(p).endswith(".gz") else open
    with opener(p,"rt",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--de-all",required=True);ap.add_argument("--mapping",required=True,help="features mapping TSV with approved_symbol, stable_gene_id, decision")
    ap.add_argument("--phase-config",default="config/phase_positive_control.json");ap.add_argument("--cohort",required=True,choices=["discovery","heca"])
    ap.add_argument("--out",required=True);ap.add_argument("--summary",required=True);a=ap.parse_args()
    pcfg=json.loads(pathlib.Path(a.phase_config).read_text());panel=pcfg["endpoints"]["S1"]["panel"];thresholds=pcfg["endpoints"]["S1"]["pass"]
    ids={}
    for r in tsv(a.mapping):
        if r["decision"]=="include":ids.setdefault(r["approved_symbol"],set()).add(r["stable_gene_id"])
    de={r["gene_id"]:r for r in tsv(a.de_all)}
    rows=[]
    for item in panel:
        sym=item["gene"];exp=1 if item["expected_sign"]=="+" else -1;cand=sorted(ids.get(sym,()))
        row={"gene":sym,"expected_sign":item["expected_sign"],"stable_gene_id":cand[0] if len(cand)==1 else ("|".join(cand) if cand else "unmapped"),"source":item["source"],"in_expression_family":False,"fit_estimable":None,"logFC":None,"t":None,"P.Value":None,"adj.P.Val":None,"CI_low":None,"CI_high":None,"observed_sign":None,"expected_sign_agrees":None,"q_lt_0_05":None}
        if len(cand)==1 and cand[0] in de:
            d=de[cand[0]];row["in_expression_family"]=True;row["fit_estimable"]=d["fit_estimable"]=="TRUE"
            for k in ("logFC","t","P.Value","adj.P.Val","CI_low","CI_high"):row[k]=d[k]
            if row["fit_estimable"]:
                lf=float(d["logFC"]);obs=(lf>0)-(lf<0);row["observed_sign"]=obs;row["expected_sign_agrees"]=obs==exp
                row["q_lt_0_05"]=d["adj.P.Val"]!="NA" and float(d["adj.P.Val"])<0.05
        rows.append(row)
    fields=list(rows[0]);pathlib.Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    with open(a.out,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader()
        for r in rows:w.writerow({k:("NA" if v is None else v) for k,v in r.items()})
    agree=sum(1 for r in rows if r["expected_sign_agrees"] is True);agree_q=sum(1 for r in rows if r["expected_sign_agrees"] is True and r["q_lt_0_05"] is True)
    if a.cohort=="discovery":need=thresholds["discovery_expected_sign_and_q_lt_0_05_min"];passed=agree_q>=need;criterion="expected sign and BH q<0.05"
    else:need=thresholds["heca_expected_sign_min"];passed=agree>=need;criterion="expected sign"
    summary={"cohort":a.cohort,"panel_genes":len(rows),"in_expression_family":sum(r["in_expression_family"] for r in rows),"expected_sign_agrees":agree,"expected_sign_agrees_and_q_lt_0_05":agree_q,"criterion":criterion,"threshold":need,"pass":passed,"filter_failures_count_as_failures":True,"genes":[{k:r[k] for k in ("gene","expected_sign","stable_gene_id","in_expression_family","logFC","adj.P.Val","expected_sign_agrees","q_lt_0_05")} for r in rows],"dependency_note":pcfg["endpoints"]["S1"].get("dependency_note")}
    pathlib.Path(a.summary).write_text(json.dumps(summary,indent=2)+"\n");print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
