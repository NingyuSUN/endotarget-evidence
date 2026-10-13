#!/usr/bin/env python3
"""Report all prespecified clinical subsets and untouched original ascertainment."""
import csv
import datetime
import json
import pathlib
import sys
from evaluate_replication import (read_tsv,index_unique,number,boolean,sign,valid_fit,
    validate_family,evaluate_set,sha,DE_COLUMNS)
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/"work/statistical/clinical_sensitivities"
METHODS=("B0","B1","G1")
MISSING={"","NA","nan","NaN","None","unknown","Unknown"}

def write_tsv(path,rows):
    if not rows:raise ValueError("Empty report table")
    with path.open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows({k:"NA" if v is None else v for k,v in r.items()} for r in rows)

def ascertainment_row(cohort,clinical,candidate):
    if cohort=="discovery":
        if candidate is None:raise ValueError("Discovery original clinical join missing")
        age=candidate["Patient.age"];day=candidate["Menstrual cycle day (Base 1)"]
        progesterone=candidate["Progesterone in serum (nmol/L)"];histology=candidate["pathologist HE evaluation"]
        expression="";author=candidate["MenstrualCyclePhase"]
        strict=boolean(clinical["strict_histology_sensitivity"]);free=boolean(clinical["adenomyosis_free_sensitivity"])
        hormone=candidate["Hormonal treatment/IUD"];h5=""
        expression_status="separate_expression_phase_field_not_recorded_in_cached_discovery_table"
    elif cohort=="heca":
        age=clinical["age"];day=clinical["day_original"];progesterone=""
        histology=clinical["histology_phase_original"];expression=clinical["expression_phase_original"]
        author=clinical["phase_original"];strict=not boolean(clinical["histology_conflict"]);free=float(clinical["adenomyosis"])==0
        hormone=clinical["hormone_clinical_original"];h5=clinical["hormone_h5_original"]
        expression_status="source_recorded_field"
    else:raise ValueError("Unknown cohort")
    return {"cohort":cohort,"donor_id":clinical["donor_id"],"sample_id":clinical["sample_id"],"disease":clinical["disease"],
        "age_original":age,"age_frozen":clinical["age"],"age_numeric":number(age),"cycle_day_original":day,"cycle_day_numeric":number(day),
        "progesterone_original_nmol_L":progesterone,"progesterone_numeric_nmol_L":number(progesterone),
        "histology_original":histology,"frozen_histology_phase_original":clinical["histology_phase_original"],
        "expression_phase_original":expression,"expression_phase_field_status":expression_status,"author_phase_original":author,
        "strict_histology_membership":strict,"adenomyosis_free_membership":free,"hormone_clinical_original":hormone,"hormone_h5_original":h5,
        "hormone_eligibility_basis":clinical["hormone_eligibility_basis"]}

def ascertainment_summaries(rows):
    result=[]
    fields=(("age","age_original","age_numeric"),("cycle_day","cycle_day_original","cycle_day_numeric"),
        ("progesterone_nmol_L","progesterone_original_nmol_L","progesterone_numeric_nmol_L"),
        ("histology","histology_original",None),("expression_phase","expression_phase_original",None),
        ("per_person_hormone_clinical","hormone_clinical_original",None))
    for cohort in ("discovery","heca"):
        for group in ("ALL","ENDO","non-ENDO"):
            subset=[r for r in rows if r["cohort"]==cohort and (group=="ALL" or r["disease"]==group)]
            for field,raw,numeric in fields:
                present=sum(r[raw] not in MISSING for r in subset)
                values=[r[numeric] for r in subset if r[numeric] is not None] if numeric else []
                result.append({"cohort":cohort,"disease":group,"field":field,"denominator_people":len(subset),"recorded_original_value":present,
                    "missing_or_unknown":len(subset)-present,"numeric_value_count":len(values) if numeric else None,"numeric_min":min(values) if values else None,"numeric_max":max(values) if values else None,
                    "note":"not separately recorded" if field=="expression_phase" and cohort=="discovery" else "not recorded in selected cached HECA table" if field=="progesterone_nmol_L" and cohort=="heca" else "source original; no imputation or relabeling"})
    return result

def paired_gene_details(primary,filters,external):
    details=[]
    for gene,p in primary.items():
        e=external.get(gene,{})
        if not valid_fit(p):reason="primary_invalid_fit"
        elif gene not in filters:reason="sensitivity_absent_id"
        elif not filters[gene]:reason="sensitivity_expression_filter_fail"
        elif not e or not valid_fit(e):reason="sensitivity_invalid_fit"
        elif sign(p["logFC"])==0:reason="primary_zero_sign"
        elif sign(e["logFC"])==0:reason="sensitivity_zero_sign"
        else:reason=""
        ps=sign(p["logFC"]);es=sign(e.get("logFC"));low=number(e.get("CI_low"));high=number(e.get("CI_high"))
        row={"gene_id":gene,"primary_logFC":number(p["logFC"]),"sensitivity_logFC":number(e.get("logFC")),"primary_sign":ps,"sensitivity_sign":es,
             "direction_evaluable":not reason,"reason":reason,"same_direction":ps==es if not reason else None,
             "opposite_direction_unadjusted_CI_excludes_zero":bool(not reason and ps!=es and low is not None and high is not None and (low>0 or high<0))}
        for prefix,source in (("primary",p),("sensitivity",e)):
            for key in ("t","P.Value","adj.P.Val","CI_low","CI_high"):row[prefix+"_"+key]=number(source.get(key))
        details.append(row)
    return details

def summarize_pairs(rows,label):
    evaluable=[r for r in rows if r["direction_evaluable"]]
    same=sum(r["same_direction"] for r in evaluable)
    return {"comparison_scope":label,"primary_genes":len(rows),"direction_evaluable":len(evaluable),"same_direction":same,"opposite_direction":len(evaluable)-same,
        "unevaluable":len(rows)-len(evaluable),"same_direction_rate":same/len(evaluable) if evaluable else None,
        "opposite_direction_unadjusted_CI_excludes_zero":sum(r["opposite_direction_unadjusted_CI_excludes_zero"] for r in evaluable),
        "interpretation":"descriptive paired effect stability; missing genes are not opposite effects; CI flag not multiplicity-adjusted significance"}

def read_model(directory):
    filters={r["gene_id"]:boolean(r["expression_pass"]) for r in read_tsv(directory/"expression_filter.tsv",("gene_id","expression_pass"))}
    if len(filters)!=len(read_tsv(directory/"expression_filter.tsv")):raise ValueError("Duplicate filter IDs")
    family=read_tsv(directory/"de_all.tsv",DE_COLUMNS);fitted=read_tsv(directory/"de_fitted.tsv",DE_COLUMNS)
    return filters,validate_family(filters,family,fitted)

def main():
    os_paths=[]
    snapshot=json.loads((OUT/"input_snapshot.json").read_text())
    for p,h in snapshot["input_sha256"].items():
        if sha(ROOT/p)!=h:raise ValueError("Sensitivity source changed: "+p)
        os_paths.append(ROOT/p)
    summary=json.loads((OUT/"summary.json").read_text());process=json.loads((OUT/"process.json").read_text())
    if summary["models_fit"]!=3 or process["status"]!="COMPLETED" or process["exit_status"]!=0:raise ValueError("Sensitivity models not complete")
    freeze=json.loads((ROOT/"receipts/discovery_rank_freeze.json").read_text())
    for p,h in freeze["artifacts"].items():
        if sha(ROOT/p)!=h:raise ValueError("Discovery frozen artifact changed")
    lists={(m,k):read_tsv(ROOT/freeze["products"][m+"_top"+str(k)],("gene_id","discovery_sign","rank")) for m in METHODS for k in (20,50)}
    for m,k in lists:
        if len(lists[m,k])!=k:raise ValueError("Frozen list size mismatch")
    reports=ROOT/"reports";generated=[]
    scenario_rows=[];robust=[];external_rows=[];external_details=[]
    discovery_filters,primary_discovery=read_model(ROOT/"work/statistical/discovery_primary/de")
    for name,status in summary["scenarios"].items():
        audit=status["audit"]
        scenario_rows.append({"cohort":status["cohort"],"scenario":status["scenario"],"status":status["status"],"fit_executed":status["fit_executed"],
            "donors":audit.get("n_donors"),"ENDO":audit["ENDO"],"non_ENDO":audit["non_ENDO"],"rank":audit["rank"],"residual_df":audit["residual_df"],
            "expression_pass_genes":status.get("expression_pass_family"),"estimable_genes":status.get("estimable_genes"),"BH_q_lt05":status.get("BH05"),
            "minimum_BH_q":status.get("minimum_BH_q"),"reason":audit.get("reason","")})
        if not status["fit_executed"]:
            if (OUT/name/"model.rds").exists() or (OUT/name/"de_all.tsv").exists():raise ValueError("Failed minimum was fitted")
            continue
        provenance=json.loads((OUT/name/"provenance.json").read_text())
        for p,h in provenance["output_sha256"].items():
            if sha(OUT/name/p)!=h:raise ValueError("Sensitivity output integrity failed")
        filters,external=read_model(OUT/name)
        if status["cohort"]=="discovery":
            details=paired_gene_details(primary_discovery,filters,external)
            fullpath=OUT/name/"paired_primary_gene_effects.tsv";write_tsv(fullpath,details);generated.append(fullpath)
            by_gene=index_unique(details)
            robust.append({"scenario":status["scenario"],**summarize_pairs(details,"all_primary_expression_pass_genes")})
            for m,k in lists:
                selected=[by_gene[r["gene_id"]] for r in lists[m,k]]
                robust.append({"scenario":status["scenario"],**summarize_pairs(selected,m+"_top"+str(k))})
        else:
            for m,k in lists:
                details,result=evaluate_set(lists[m,k],filters,external,k)
                external_rows.append({"scenario":status["scenario"],"method":m,"k":k,"selected":result["selected"],"sensitivity_evaluable":result["evaluable"],"agree":result["agree"],
                    "opposite":result["opposite"],"unevaluable":result["unevaluable"],"concordance":result["concordance"],"coverage":result["coverage"],"coverage_adequate":result["coverage_adequate"],
                    "missing_bound_low":result["missing_bounds"][0],"missing_bound_high":result["missing_bounds"][1],"significant_same_direction_BH05":result["significant_same_direction"],"lists_replaced":False})
                for r in details:
                    r.pop("original_evaluable");r["sensitivity_evaluable"]=r.pop("evaluable")
                    external_details.append({"scenario":status["scenario"],"method":m,"k":k,**r})
    discovery=read_tsv(ROOT/"manifests/final_discovery_clinical_selected.tsv");heca=read_tsv(ROOT/"manifests/final_heca_clinical_selected.tsv")
    candidates=index_unique(read_tsv(ROOT/"manifests/discovery_clinical_candidates.tsv"),"sample_id")
    ascertain=[ascertainment_row("discovery",r,candidates.get(r["sample_id"])) for r in discovery]+[ascertainment_row("heca",r,None) for r in heca]
    if len(ascertain)!=44 or len({(r["cohort"],r["donor_id"]) for r in ascertain})!=44:raise ValueError("Exact44ascertainment people required")
    if any(r["age_numeric"] is None for r in ascertain):raise ValueError("Complete-age audit changed; do not invent unadjusted model")
    tables={"clinical_sensitivity_scenarios.tsv":scenario_rows,"clinical_sensitivity_discovery_robustness.tsv":robust,
        "clinical_sensitivity_heca_fixed_list_summary.tsv":external_rows,"clinical_sensitivity_heca_fixed_list_details.tsv":external_details,
        "clinical_sensitivity_ascertainment_donors.tsv":ascertain,"clinical_sensitivity_ascertainment_summary.tsv":ascertainment_summaries(ascertain)}
    for name,rows in tables.items():p=reports/name;write_tsv(p,rows);generated.append(p)
    lines=["# 临床敏感性与分期证据报告","","已执行三个预先定义的可估计子集；两个人数不足的组合没有拟合。所有分析仍是观察性供者表达关联，不能证明病因或治疗作用。","",
        "本轮包装脚本实现及执行发生在主发现和外部结果已观察之后。子集定义、资格标记和设计文件此前已存在，并在本轮首次拟合前保存哈希；这是可审计的预设敏感性执行，不是预注册或盲法验证。","",
        "|队列/场景|人数（ENDO/非ENDO）|秩/残余自由度|通过表达过滤基因|BH q<0.05|状态|","|---|---:|---:|---:|---:|---|"]
    for r in scenario_rows:
        lines.append("|{cohort}/{scenario}|{donors}（{ENDO}/{non_ENDO}）|{rank}/{residual_df}|{expression_pass_genes}|{BH_q_lt05}|{status}|".format(**{k:"NA" if v is None else v for k,v in r.items()}))
    lines += ["","过滤后的基因数会随原规则重新计算：严格发现子集的较小组由10人变为9人，表达支持人数也按同一规则重算，因此供者更少但通过过滤的基因可以略多。这不是根据显著性放宽规则。","","主分析与发现子集逐基因比较，分母是能在两次拟合中得到有效非零方向的基因；过滤失败或缺失另计，不算反方向。完整逐基因效果表保存在各发现子集的 paired_primary_gene_effects.tsv，原发现名单及方向没有更换。","",
        "|发现子集|原主分析基因|方向可比较|同方向|反方向|不可比较|","|---|---:|---:|---:|---:|---:|"]
    for r in robust:
        if r["comparison_scope"]=="all_primary_expression_pass_genes":lines.append("|{scenario}|{primary_genes}|{direction_evaluable}|{same_direction}|{opposite_direction}|{unevaluable}|".format(**r))
    lines += ["","无腺肌症 HECA 子集仍只评估冻结原名单，每个方法的可评估分母单独报告。新的子集结果不会改写主外部 bootstrap 的原固定分母。","",
        "|固定名单|同方向/可评估/选择数|同方向率|缺失上下界|BH q<0.05 同方向|","|---|---:|---:|---:|---:|"]
    for r in external_rows:lines.append("|{} top{}|{}/{}/{}|{:.4f}|[{:.4f}, {:.4f}]|{}|".format(r["method"],r["k"],r["agree"],r["sensitivity_evaluable"],r["selected"],r["concordance"],r["missing_bound_low"],r["missing_bound_high"],r["significant_same_direction_BH05"]))
    lines += ["","44 人年龄均有原始记录。发现组周期日22/22可读为数值；HECA周期日16/22，六人未知。原始日28的对照被如实保留，不做新的周期日排除。","",
        "发现组组织学12/22缺失，两个显式分期矛盾在预设严格子集中排除；缺失组织学并没有因此获得独立确认。HECA组织学11人Proliferative、8人Secretory、1人混合、2人缺失；表达辅助分期12/22有记录、10/22为空。发现缓存没有单独的逐人表达辅助分期字段，不从作者最终标签补造。","",
        "HECA个人激素记录NA/nan保留为缺失；自然周期资格依据来自研究队列定义，不是把每人缺失值改写为未用激素。全部原值和每组分母见 ascertainment 表。","",
        "30/100细胞阈值与主队列选择相同，复用主结果；这不是两次独立验证。所有年龄完整，因此未新做去掉年龄模型；也没有事后制定周期日界限。子集减少供者，组织学/作者表达辅助分期依赖、病例其他病理及残余混杂仍限制结论。","",
        "YELLOW — 理解尚未评估：应能区分组织学矛盾排除与独立组织学确认；供者子集稳健性与因果证明；完整年龄记录与年龄混杂完全消除。"]
    report=reports/"clinical_sensitivity_report.md";report.write_text("\n".join(lines)+"\n");generated.append(report)
    input_paths=[OUT/"input_snapshot.json",OUT/"summary.json",ROOT/"manifests/discovery_clinical_candidates.tsv",ROOT/"receipts/discovery_rank_freeze.json",pathlib.Path(__file__),ROOT/"scripts/evaluate_replication.py"]
    receipt={"status":"COMPLETE_PRESPECIFIED_CLINICAL_SENSITIVITIES","completed_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"models_fit":3,"scenarios_not_estimable":2,"scenario_summary":scenario_rows,"process_receipt":process,
        "primary_external_effects_previously_observed":True,"preregistered":False,"original_discovery_lists_replaced":False,"original_bootstrap_masks_replaced":False,"clinical_ascertainment_people":44,
        "new_age_unadjusted_model":False,"new_cycle_day_cutoff":False,"threshold30_100":"identicalprimaryselection_reuse_no_refits",
        "input_sha256":{str(p.relative_to(ROOT)):sha(p) for p in input_paths},"output_sha256":{str(p.relative_to(ROOT)):sha(p) for p in generated},"scientific_scope":"observationalclinicalsensitivity_notcausalvalidation","scientific_validation":False}
    target=ROOT/"receipts/clinical_sensitivity_execution.json"
    if target.exists():raise ValueError("Existing sensitivity execution receipt preserved")
    target.write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"models_fit":3,"not_estimable":2,"ascertainment_people":44,"reports":len(generated)},indent=2))
if __name__=="__main__":
    try:main()
    except (ValueError,OSError,KeyError) as exc:sys.exit("Sensitivity summary failed: "+str(exc))
