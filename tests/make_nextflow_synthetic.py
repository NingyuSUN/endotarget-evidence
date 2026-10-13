#!/usr/bin/env python3
"""Create artificial data only under work/synthetic; never open real counts."""
import csv
import datetime
import gzip
import hashlib
import json
import pathlib
import uuid

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
CODE = [
    "main.nf", "nextflow.config", "scripts/fit_donor_de.R", "scripts/donor_de_lib.R",
    "scripts/rank_targets.py", "scripts/sample_matched_nulls.py",
    "scripts/statistical_runtime_env.sh", "config/statistical_analysis.json",
    "config/ranking_algorithm.json", "requirements/statistical-runtime.explicit.txt",
    "requirements/statistical-runtime.yml", "receipts/runtime_statistics_install.json",
    "receipts/runtime_statistics_R_smoke.json",
]
def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def write_table(path, header, rows):
    with path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)

def main():
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    folder = ROOT / "work/synthetic" / (stamp + "-" + uuid.uuid4().hex[:8])
    folder.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(42)
    n_genes, n_donors = 1500, 20
    genes = [f"ENSG{i:011d}" for i in range(1, n_genes + 1)]
    donors = [f"SYN{i:02d}" for i in range(n_donors)]
    group = np.array([0] * 10 + [1] * 10)
    age = rng.normal(0, 5, n_donors)
    base = rng.lognormal(5, .65, n_genes)
    effect = rng.normal(0, .35, n_genes)
    means = base[:, None] * np.exp(effect[:, None] * group + rng.normal(0, .12, (n_genes,n_donors)))
    counts = rng.poisson(means)
    with gzip.open(folder/"counts.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *donors])
        writer.writerows([[gene, *row.tolist()] for gene,row in zip(genes,counts)])
    write_table(folder/"design.tsv", ["donor_id","Intercept","ENDO_minus_non_ENDO","age_centered"],
                [[d,1,int(g),float(a)] for d,g,a in zip(donors,group,age)])
    write_table(folder/"clinical.tsv", ["donor_id","disease","primary_lineage_cells","age","capture_batch","adenomyosis"],
                [[d,"ENDO" if g else "non-ENDO",100,35+float(a),"synthetic",False] for d,g,a in zip(donors,group,age)])
    with gzip.open(folder/"edges.tsv.gz", "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene1","gene2","weight"])
        writer.writerows([[genes[i],genes[(i+1)%n_genes],.7+.05*(i%7)] for i in range(n_genes)])
    inputs = {key:str((folder/key).relative_to(ROOT)) for key in ["counts.tsv.gz","design.tsv","clinical.tsv","edges.tsv.gz"]}
    artifact_paths = CODE + list(inputs.values())
    frozen = {"schema_version":1,"phase":"SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE",
              "data_kind":"synthetic","scientific_validation":False,
              "created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "artifacts":{name:sha(ROOT/name) for name in artifact_paths}}
    (folder/"freeze.json").write_text(json.dumps(frozen,indent=2)+"\n")
    params = {"synthetic_test":True,"pre_effect_freeze":str((folder/"freeze.json").relative_to(ROOT)),
              "counts":inputs["counts.tsv.gz"],"design":inputs["design.tsv"],
              "clinical":inputs["clinical.tsv"],"edges":inputs["edges.tsv.gz"]}
    (folder/"params.json").write_text(json.dumps(params,indent=2)+"\n")
    (folder/"fixture.json").write_text(json.dumps({"data_kind":"synthetic","seed":42,
        "genes":n_genes,"donors":n_donors,"non_ENDO":10,"ENDO":10,
        "scope":"flow_mechanics_only_not_real_disease_or_target_validation"},indent=2)+"\n")
    print(folder)
if __name__ == "__main__":
    main()
