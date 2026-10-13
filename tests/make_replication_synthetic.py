#!/usr/bin/env python3
"""Create artificial HECA inputs paired with completed artificial discovery outputs."""
import argparse, csv, datetime, gzip, hashlib, json, pathlib, uuid
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
CODE=["replication.nf","replication.config","scripts/fit_donor_de.R","scripts/donor_de_lib.R",
      "scripts/evaluate_replication.py","scripts/statistical_runtime_env.sh",
      "config/statistical_analysis.json","requirements/statistical-runtime.explicit.txt"]
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def rel(p):return str(pathlib.Path(p).resolve().relative_to(ROOT))
def table(p,columns,rows):
    with p.open("x",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(columns);w.writerows(rows)
def dump(p,o):
    with p.open("x") as f:json.dump(o,f,indent=2);f.write("\n")
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--discovery-results",type=pathlib.Path,required=True)
    a=ap.parse_args();dis=a.discovery_results.resolve()
    if not dis.is_relative_to(ROOT/"work/synthetic"):raise ValueError("artificial discovery required")
    dp=json.loads((dis/"de/provenance.json").read_text())
    if dp.get("data_kind")!="synthetic_test_only":raise ValueError("not synthetic")
    genes=[]
    with (dis/"de/expression_filter.tsv").open() as f:
        genes=[r["gene_id"] for r in csv.DictReader(f,delimiter="\t")]
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    folder=ROOT/"work/synthetic"/(stamp+"-replication-"+uuid.uuid4().hex[:8]);folder.mkdir()
    rng=np.random.default_rng(31415);donors=[f"TOYHECA{i:02d}" for i in range(22)]
    group=np.array([0]*7+[1]*15);age=rng.normal(0,5,22)
    # Independent artificial effects, missing original IDs, and low-count genes.
    genes=genes[:-20];base=rng.lognormal(5,.6,len(genes));effect=rng.normal(0,.2,len(genes))
    means=base[:,None]*np.exp(effect[:,None]*group+rng.normal(0,.15,(len(genes),22)))
    counts=rng.poisson(means);counts[-30:]=0
    with gzip.open(folder/"counts.tsv.gz","wt") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["gene_id",*donors]);w.writerows([[g,*row] for g,row in zip(genes,counts)])
    table(folder/"design.tsv",["donor_id","Intercept","ENDO_minus_non_ENDO","age_centered"],
          [[d,1,int(g),float(a)] for d,g,a in zip(donors,group,age)])
    table(folder/"clinical.tsv",["donor_id","disease","primary_nuclei","age","library_group_Lane","adenomyosis"],
          [[d,"ENDO" if g else "non-ENDO",100,35+float(a),"synthetic",False] for d,g,a in zip(donors,group,age)])
    inputs=[rel(folder/n) for n in ("counts.tsv.gz","design.tsv","clinical.tsv")]
    pre={"phase":"SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE","scientific_validation":False,
         "artifacts":{n:sha(ROOT/n) for n in CODE+inputs}}
    dump(folder/"pre_effect.json",pre)
    products={"all_U_directions":rel(dis/"nulls/all_U_directions.tsv"),
              "null_memberships":rel(dis/"nulls/null_memberships.tsv"),
              "null_status":rel(dis/"nulls/status.json"),
              "statistical_config":"config/statistical_analysis.json"}
    for m in ("B0","B1","G1"):
        for role,suffix in (("full",""),("top20","_top20"),("top50","_top50")):
            products[m+"_"+role]=rel(dis/"ranking"/(m+suffix+".tsv"))
    paths=set(CODE+inputs+[rel(folder/"pre_effect.json")]+list(products.values()))
    freeze={"phase":"SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE","scientific_validation":False,
            "external_disease_effects_computed":False,"products":products,
            "artifacts":{n:sha(ROOT/n) for n in sorted(paths)}}
    dump(folder/"discovery_freeze.json",freeze)
    paths.add(rel(folder/"discovery_freeze.json"))
    dump(folder/"evaluation_freeze.json",{"phase":"SYNTHETIC_REPLICATION_TEST_ONLY",
         "scientific_validation":False,"artifacts":{n:sha(ROOT/n) for n in sorted(paths)}})
    params=dict(synthetic_test=True,counts=inputs[0],design=inputs[1],clinical=inputs[2],
                pre_effect_freeze=rel(folder/"pre_effect.json"),discovery_freeze=rel(folder/"discovery_freeze.json"),
                evaluation_freeze=rel(folder/"evaluation_freeze.json"))
    dump(folder/"params.json",params)
    dump(folder/"fixture.json",{"data_kind":"synthetic_test_only","n_genes":len(genes),"n_donors":22,
         "ENDO":15,"non_ENDO":7,"seed":31415,"source_discovery":rel(dis),"scientific_validation":False})
    print(rel(folder))
if __name__=="__main__":main()
