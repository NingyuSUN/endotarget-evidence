#!/usr/bin/env python3
"""Evaluate immutable discovery predictions against complete external donor-level fits."""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import math
import pathlib
import platform
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
METHODS = ("B0", "B1", "G1")
REASONS = ("absent_id", "expression_filter_fail", "invalid_fit", "discovery_zero_sign", "external_zero_sign")
DETAIL_FIELDS = ["method","k","rank","gene_id","discovery_sign","evaluable","original_evaluable",
                 "missing_reason","external_sign","external_logFC","external_t","external_p","external_q",
                 "external_CI_low","external_CI_high","agree","signed_external_t"]
DE_COLUMNS = ("gene_id","logFC","t","P.Value","adj.P.Val","CI_low","CI_high","fit_estimable")

def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4*1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def project_path(path):
    p = pathlib.Path(path)
    if not p.is_absolute():
        p = ROOT / p
    p = p.resolve()
    p.relative_to(ROOT)
    return p

def read_json(path):
    return json.loads(pathlib.Path(path).read_text())

def read_tsv(path, required=()):
    with pathlib.Path(path).open() as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames or not set(required).issubset(r.fieldnames):
            raise ValueError("Missing TSV columns: " + str(path))
        if len(set(r.fieldnames)) != len(r.fieldnames):
            raise ValueError("Duplicate TSV columns: " + str(path))
        ans = list(r)
        if any(None in row or any(v is None for v in row.values()) for row in ans):
            raise ValueError("Malformed TSV row: " + str(path))
        return ans

def index_unique(rows, key="gene_id"):
    out = {}
    for row in rows:
        if not row[key] or row[key] in out:
            raise ValueError("Empty or duplicate " + key)
        out[row[key]] = row
    return out

def number(value):
    if value is None or str(value).strip() in ("", "NA", "NaN", "nan", "None"):
        return None
    try:
        n = float(value)
    except (ValueError, TypeError):
        return None
    return n if math.isfinite(n) else None

def boolean(value):
    if value is True or str(value).upper() == "TRUE":
        return True
    if value is False or str(value).upper() == "FALSE":
        return False
    raise ValueError("Invalid logical field: " + str(value))

def sign(value):
    n = number(value)
    return None if n is None else (1 if n > 0 else -1 if n < 0 else 0)

def valid_fit(row):
    return (boolean(row["fit_estimable"]) and
            all(number(row[k]) is not None for k in ("logFC","t","P.Value")) and
            0 <= number(row["P.Value"]) <= 1)

def bh_complete_family(rows):
    # n includes all expression-pass genes, including failed fits.
    valid = [(row["gene_id"], number(row["P.Value"])) for row in rows if valid_fit(row)]
    valid.sort(key=lambda item: (item[1], item[0]))
    adjusted = {}
    upper = 1.
    for i in range(len(valid)-1, -1, -1):
        gene, p = valid[i]
        upper = min(upper, p*len(rows)/(i+1))
        adjusted[gene] = upper
    return {row["gene_id"]: adjusted.get(row["gene_id"]) for row in rows}

def validate_bh(rows):
    expected = bh_complete_family(rows)
    for row in rows:
        q = number(row["adj.P.Val"])
        if expected[row["gene_id"]] is None:
            if q is not None:
                raise ValueError("A failed fit must have NA q, never zero")
        elif q is None or not math.isclose(q, expected[row["gene_id"]], rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("BH q does not match complete testing family")
    return expected

def validate_family(filters, family, fitted):
    by_gene = index_unique(family)
    expected = {gene for gene, keep in filters.items() if keep}
    if set(by_gene) != expected:
        raise ValueError("DE does not equal full expression-pass family")
    valid = {gene: row for gene, row in by_gene.items() if valid_fit(row)}
    fitted_index = index_unique(fitted)
    if set(fitted_index) != set(valid) or any(fitted_index[g] != valid[g] for g in valid):
        raise ValueError("de_fitted does not equal exact valid-fit subset of de_all")
    for row in family:
        if boolean(row["fit_estimable"]) != valid_fit(row):
            raise ValueError("fit_estimable conflicts with finite effect/t/p")
    validate_bh(family)
    return by_gene

def classify(gene, discovery_sign, filters, external):
    if gene not in filters:
        return "absent_id"
    if not filters[gene]:
        return "expression_filter_fail"
    if gene not in external or not valid_fit(external[gene]):
        return "invalid_fit"
    if discovery_sign == 0:
        return "discovery_zero_sign"
    if sign(external[gene]["logFC"]) == 0:
        return "external_zero_sign"
    return ""

def evaluate_set(selected, filters, external, k):
    if len(selected) != k or len({r["gene_id"] for r in selected}) != k:
        raise ValueError("Selected set has wrong size or duplicate IDs")
    details = []
    for rank, selected_row in enumerate(selected, 1):
        gene = selected_row["gene_id"]
        ds = int(selected_row["discovery_sign"])
        if ds not in (-1,0,1):
            raise ValueError("Invalid frozen discovery sign")
        reason = classify(gene, ds, filters, external)
        er = external.get(gene, {})
        fields = {dest: number(er.get(src)) for dest,src in (
            ("external_logFC","logFC"),("external_t","t"),("external_p","P.Value"),
            ("external_q","adj.P.Val"),("external_CI_low","CI_low"),("external_CI_high","CI_high"))}
        es = sign(er.get("logFC"))
        details.append({"rank": int(selected_row.get("rank",rank)), "gene_id": gene,
                        "discovery_sign": ds,"evaluable": not reason,"original_evaluable": not reason,
                        "missing_reason": reason,"external_sign": es,
                        "agree": ds == es if not reason else None,
                        "signed_external_t": fields["external_t"]*ds if er and valid_fit(er) and ds != 0 else None, **fields})
    evaluable = sum(r["evaluable"] for r in details)
    agree = sum(r["agree"] is True for r in details)
    missing = k-evaluable
    return details, {"selected":k, "evaluable":evaluable,"agree":agree,
                     "opposite":evaluable-agree,"unevaluable":missing,
                     "concordance":agree/evaluable if evaluable else None,
                     "coverage":evaluable/k, "coverage_adequate":evaluable >= math.ceil(.8*k),
                     "missing_bounds":[agree/k,(agree+missing)/k],
                     "missing_reasons":{reason:sum(r["missing_reason"]==reason for r in details) for reason in REASONS},
                     "significant_same_direction":sum(r["agree"] is True and r["external_q"] is not None and
                                                     r["external_q"]<.05 for r in details)}

def primary_contrast(b0, g1):
    lo = g1["missing_bounds"][0]-b0["missing_bounds"][1]
    hi = g1["missing_bounds"][1]-b0["missing_bounds"][0]
    adequate = min(b0["evaluable"],g1["evaluable"]) >= 16
    diff = None if b0["concordance"] is None or g1["concordance"] is None else g1["concordance"]-b0["concordance"]
    return {"contrast":"G1_minus_B0_top20","concordance_difference":diff,"missing_bounds":[lo,hi],
            "B0_evaluable":b0["evaluable"],"G1_evaluable":g1["evaluable"],
            "coverage_adequate_both":adequate,"missing_bounds_include_zero":lo <= 0 <= hi,
            "eligible_for_bootstrap_superiority_assessment": adequate and lo>0,
            "superiority_claim":False,
            "interpretation":"fixed-list transcriptomic comparison; paired bootstrap still required, not therapeutic validation"}

def null_distribution(summaries, observed):
    adequate = [s["concordance"] for s in summaries if s["coverage_adequate"] and s["concordance"] is not None]
    return {"total_sets":len(summaries),"adequate_sets":len(adequate),
            "inadequate_sets":len(summaries)-len(adequate),
            "adequate_concordances":adequate,
            "observed_descriptive_percentile":100*sum(x<=observed for x in adequate)/len(adequate)
                if adequate and observed is not None else None,
            "percentile_rule":"100 times fraction of coverage-adequate sets with concordance <= observed; ties included",
            "p_value":None,"interpretation":"descriptive expression/degree matching comparison, not a permutation test"}

def dossier_category(row):
    if row["evaluable"] and row["agree"] is True and row["external_q"] is not None and row["external_q"]<.05:
        return "supported_same_direction"
    low,high = row["external_CI_low"],row["external_CI_high"]
    if row["evaluable"] and row["agree"] is False and low is not None and high is not None and (low>0 or high<0):
        return "strongly_conflicting"
    return "inconclusive"

def select_dossiers(pool):
    ordered = sorted(pool,key=lambda r:(r["B0_rank"],r["gene_id"]))
    output = []
    for category,limit in (("supported_same_direction",3),("strongly_conflicting",1),("inconclusive",1)):
        eligible = [r for r in ordered if dossier_category(r)==category]
        output.extend({**r,"dossier_category":category} for r in eligible[:limit])
    return output

def verify_discovery(path, synthetic=False):
    manifest = read_json(path)
    phase = "SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE" if synthetic else "DISCOVERY_FROZEN_PRE_EXTERNAL"
    if manifest.get("phase") != phase or manifest.get("external_disease_effects_computed") is not False:
        raise ValueError("Invalid pre-external discovery freeze")
    if not manifest.get("artifacts"):
        raise ValueError("Empty discovery freeze")
    for name,digest in manifest["artifacts"].items():
        if sha(project_path(name)) != digest:
            raise ValueError("Changed frozen artifact: " + name)
    products = manifest["products"]
    for role in ["all_U_directions","null_memberships","null_status","statistical_config"] + [
            m+"_"+suffix for m in METHODS for suffix in ("full","top20","top50")]:
        if products.get(role) not in manifest["artifacts"]:
            raise ValueError("Missing frozen mandatory role: " + role)
    all_directions = index_unique(read_tsv(project_path(products["all_U_directions"]),
                                         ("gene_id","logFC","discovery_sign")))
    for gene,row in all_directions.items():
        if int(row["discovery_sign"]) != sign(row["logFC"]):
            raise ValueError("Saved discovery sign/logFC mismatch")
    lists,full_ranks = {},{}
    for method in METHODS:
        full = read_tsv(project_path(products[method+"_full"]),("rank","gene_id","method","discovery_sign","logFC"))
        if set(index_unique(full)) != set(all_directions):
            raise ValueError("Full rankings differ from fixed all-U directions")
        for rank,row in enumerate(full,1):
            ds = all_directions[row["gene_id"]]
            if int(row["rank"]) != rank or row["method"] != method or int(row["discovery_sign"]) != int(ds["discovery_sign"]):
                raise ValueError("Invalid full frozen ranking or direction")
            if not math.isclose(float(row["logFC"]),float(ds["logFC"]),rel_tol=1e-12,abs_tol=1e-14):
                raise ValueError("Frozen ranking coefficient differs from all-U prediction")
        full_ranks[method] = {row["gene_id"]:int(row["rank"]) for row in full}
        for k in (20,50):
            top = read_tsv(project_path(products[f"{method}_top{k}"]))
            if top != full[:k] or len(top) != k:
                raise ValueError("Frozen list is not exact full-rank prefix")
            lists[(method,k)] = top
    return manifest,lists,full_ranks,all_directions

def validate_external(directory, discovery_path, manifest, synthetic=False):
    prov = read_json(directory/"provenance.json")
    if prov.get("cohort") != "heca" or prov.get("data_kind") != ("synthetic_test_only" if synthetic else "public_cohort"):
        raise ValueError("Wrong external cohort/data kind")
    if prov.get("scientific_effect_run") is not (not synthetic):
        raise ValueError("Inconsistent scientific effect-run provenance")
    counts_path = "data/processed/heca/counts_stable.tsv.gz"
    required = {counts_path,"manifests/final_heca_clinical_selected.tsv",
                "reports/final_clinical_design_heca_primary.tsv","config/statistical_analysis.json",
                "scripts/fit_donor_de.R","scripts/donor_de_lib.R"}
    if synthetic:
        required = set(prov.get("input_sha256",{}))
    inputs = {}
    for original,digest in prov.get("input_sha256",{}).items():
        path = project_path(original)
        relative = str(path.relative_to(ROOT))
        expected = sha(discovery_path) if path==discovery_path else manifest["artifacts"].get(relative)
        if expected is None or digest != expected or sha(path) != digest:
            raise ValueError("External provenance differs from original frozen input: " + original)
        inputs[relative] = digest
    if not required.issubset(inputs) or str(discovery_path.relative_to(ROOT)) not in inputs:
        raise ValueError("External provenance lacks locked inputs/discovery-freeze signature")
    if not synthetic and "receipts/statistical_pre_effect_freeze.json" not in inputs:
        raise ValueError("External provenance lacks original implementation freeze")
    outputs = {}
    for original,digest in prov.get("output_sha256",{}).items():
        original_path = project_path(original)
        name = original_path.name
        if name in outputs or sha(directory/name) != digest:
            raise ValueError("External output hash mismatch or ambiguous basename")
        outputs[name] = digest
    if not {"de_all.tsv","de_fitted.tsv","expression_filter.tsv"}.issubset(outputs):
        raise ValueError("Missing signed external output artifacts")
    filter_rows = read_tsv(directory/"expression_filter.tsv",("gene_id","expression_pass"))
    filters = {g:boolean(r["expression_pass"]) for g,r in index_unique(filter_rows).items()}
    # A full all-input-gene filter ledger distinguishes absent IDs from filtered IDs.
    if not synthetic:
        with gzip.open(project_path(counts_path),"rt") as f:
            counts_reader = csv.reader(f,delimiter="\t")
            header = next(counts_reader)
            if header[0] != "gene_id":
                raise ValueError("Counts gene ID column mismatch")
            counts_ids = [r[0] for r in counts_reader]
        if counts_ids != [r["gene_id"] for r in filter_rows]:
            raise ValueError("Expression filter ledger does not cover exact original count IDs/order")
        if len(counts_ids) != prov["audit"]["input_genes"] or len(header)-1 != prov["audit"]["n_donors"]:
            raise ValueError("Counts dimensions disagree with external provenance")
    family = read_tsv(directory/"de_all.tsv",DE_COLUMNS)
    fitted = read_tsv(directory/"de_fitted.tsv",DE_COLUMNS)
    external = validate_family(filters,family,fitted)
    if len(family) != prov.get("tested_family_size") or len(fitted) != prov.get("estimable_genes"):
        raise ValueError("Testing family / fit counts disagree with external provenance")
    return filters,external,prov,{"inputs":inputs,"outputs":outputs}

def evaluate_nulls(manifest, directions, lists, filters, external, observed):
    products = manifest["products"]
    status = read_json(project_path(products["null_status"]))
    rows = read_tsv(project_path(products["null_memberships"]),("method","draw","position","gene_id","discovery_sign","logFC"))
    grouped = {}
    for row in rows:
        method = row["method"]
        if method not in METHODS:
            raise ValueError("Unknown null method")
        key=(method,int(row["draw"]))
        grouped.setdefault(key,[]).append(row)
    results,summaries = [],{}
    for method in METHODS:
        state=status["methods"][method]
        draws={draw:members for (m,draw),members in grouped.items() if m==method}
        if state["status"] == "NOT_ESTIMABLE":
            if draws or state["sampled_sets"] != 0 or not state["reason"]:
                raise ValueError("Invalid frozen non-estimable null experiment")
            summaries[method]={"status":"NOT_ESTIMABLE","reason":state["reason"],**null_distribution([],observed[method])}
            continue
        if state["status"] != "ESTIMABLE" or state["sampled_sets"] != 200 or set(draws) != set(range(1,201)):
            raise ValueError("Expected exactly 200 original frozen null sets per estimable method")
        local=[]
        selected_ids={r["gene_id"] for r in lists[(method,20)]}
        for draw,members in sorted(draws.items()):
            ordered=sorted(members,key=lambda r:int(r["position"]))
            if [int(r["position"]) for r in ordered] != list(range(1,21)):
                raise ValueError("Frozen null draw positions are incomplete")
            for row in ordered:
                if row["gene_id"] not in directions or row["gene_id"] in selected_ids:
                    raise ValueError("Frozen null gene outside U or overlaps own top20")
                ds=directions[row["gene_id"]]
                if int(row["discovery_sign"]) != int(ds["discovery_sign"]) or not math.isclose(float(row["logFC"]),float(ds["logFC"]),rel_tol=1e-12,abs_tol=1e-14):
                    raise ValueError("Null member does not use its own frozen discovery prediction")
            _,summary=evaluate_set(ordered,filters,external,20)
            summary.update(method=method,draw=draw)
            local.append(summary)
            results.append(summary)
        summaries[method]={"status":"ESTIMABLE",**null_distribution(local,observed[method])}
    return results,summaries

def serial(value):
    if value is None:
        return "NA"
    if value is True:
        return "TRUE"
    if value is False:
        return "FALSE"
    return value

def write_tsv(path,rows,fields):
    with path.open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key:serial(row.get(key)) for key in fields})

def write_json(path,obj):
    with path.open("x") as f:
        json.dump(obj,f,indent=2,allow_nan=False)
        f.write("\n")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery-freeze",required=True,type=pathlib.Path)
    parser.add_argument("--external-dir",required=True,type=pathlib.Path)
    parser.add_argument("--output-dir",required=True,type=pathlib.Path)
    parser.add_argument("--synthetic-test",action="store_true",help="only fixtures wholly under work/synthetic/")
    args=parser.parse_args()
    freeze,external_dir,out=map(project_path,(args.discovery_freeze,args.external_dir,args.output_dir))
    if args.synthetic_test:
        for path in (freeze,external_dir,out):
            path.relative_to(ROOT/"work/synthetic")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be empty; previous evaluation preserved")
    manifest,lists,ranks,directions=verify_discovery(freeze,args.synthetic_test)
    # All provenance/family validations run before writing any evaluation.
    filters,external,prov,validated=validate_external(external_dir,freeze,manifest,args.synthetic_test)
    details,summaries=[],{}
    for method in METHODS:
        summaries[method]={}
        for k in (20,50):
            records,summary=evaluate_set(lists[(method,k)],filters,external,k)
            details.extend({**r,"method":method,"k":k} for r in records)
            summaries[method][str(k)]=summary
    observed={m:summaries[m]["20"]["concordance"] for m in METHODS}
    null_rows,null_summary=evaluate_nulls(manifest,directions,lists,filters,external,observed)
    pool_ids={r["gene_id"] for m in ("B0","G1") for r in lists[(m,20)]}
    predictions=[directions[g] for g in sorted(pool_ids)]
    pool,_=evaluate_set(predictions,filters,external,len(predictions))
    for r in pool:
        gene=r["gene_id"]
        r.update(B0_rank=ranks["B0"][gene],G1_rank=ranks["G1"][gene],
                 in_B0_top20=gene in {x["gene_id"] for x in lists[("B0",20)]},
                 in_G1_top20=gene in {x["gene_id"] for x in lists[("G1",20)]},
                 dossier_category=dossier_category(r))
    pool.sort(key=lambda r:(r["B0_rank"],r["gene_id"]))
    chosen=select_dossiers(pool)
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/"selected_genes.tsv",details,DETAIL_FIELDS)
    write_tsv(out/"bootstrap_interface.tsv",details,DETAIL_FIELDS)
    null_fields=["method","draw","selected","evaluable","agree","opposite","unevaluable","concordance","coverage","coverage_adequate","missing_lower","missing_upper",*REASONS]
    flattened=[{**r,"missing_lower":r["missing_bounds"][0],"missing_upper":r["missing_bounds"][1],**r["missing_reasons"]} for r in null_rows]
    write_tsv(out/"null_sets.tsv",flattened,null_fields)
    write_json(out/"null_summary.json",null_summary)
    pool_fields=["gene_id","B0_rank","G1_rank","in_B0_top20","in_G1_top20","dossier_category",*DETAIL_FIELDS[4:]]
    write_tsv(out/"dossier_pool.tsv",pool,pool_fields)
    write_tsv(out/"dossier_selection.tsv",chosen,pool_fields)
    summary={"schema_version":1,"cohort":"heca","data_kind":prov["data_kind"],
             "scientific_scope":"source-separated donor-level transcriptomic associations, not causal or therapeutic validation",
             "methods":summaries,"primary":primary_contrast(summaries["B0"]["20"],summaries["G1"]["20"]),
             "top50_role":"secondary","dossier_pool_size":len(pool),"dossier_selection_count":len(chosen),
             "dossier_conflict_rule":"opposite nonzero direction and unadjusted 95% CI excludes zero; BH q reported separately",
             "dossier_selection_rule":"fixed union of B0/G1 top20; ascending full B0 rank then stable ID; up to 3 supported, 1 conflict, 1 inconclusive; no forced category",
             "null_scope":"200 original memberships per estimable method, own saved discovery signs, descriptive only"}
    write_json(out/"summary.json",summary)
    provenance={"schema_version":1,"completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "python":platform.python_version(),"git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                "implementation_sha256":sha(pathlib.Path(__file__)),
                "discovery_freeze_sha256":sha(freeze),"external_provenance_sha256":sha(external_dir/"provenance.json"),
                "validated_external":validated,"data_kind":prov["data_kind"],
                "output_sha256":{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
    write_json(out/"provenance.json",provenance)
    print(json.dumps({"status":"REPLICATION_EVALUATION_COMPLETE","data_kind":prov["data_kind"],
                      "output_dir":str(out),"primary":summary["primary"]}))

if __name__ == "__main__":
    main()
