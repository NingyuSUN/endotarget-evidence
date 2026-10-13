#!/usr/bin/env python3
"""Read-only report assembly from frozen/evaluated artifacts. Never fits a model."""
import argparse, csv, hashlib, json, math, os, platform, re, statistics, subprocess, tempfile
from datetime import datetime, timezone
from pathlib import Path

METHODS = ("B0", "B1", "G1")
REASONS = ("absent_id", "expression_filter_fail", "invalid_fit", "discovery_zero_sign", "external_zero_sign")
def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1048576), b""): h.update(block)
    return h.hexdigest()
def load(path): return json.loads(Path(path).read_text())
def rows(path):
    with open(path, newline="") as f: return list(csv.DictReader(f, delimiter="\t"))
def truth(v): return str(v).lower() == "true"
def num(v):
    try:
        n = float(v)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError): return None
def fmt(v):
    n = num(v)
    return "NA" if n is None else f"{n:.4g}"
def pct(v): return "NA" if num(v) is None else f"{100*float(v):.1f}%"
def require(ok, message):
    if not ok: raise ValueError(message)
def unique(rs, key="gene_id"):
    out = {r[key]: r for r in rs}
    require(len(out) == len(rs), f"Duplicate {key}")
    return out
def write_tsv(path, rs, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rs)
def verify_outputs(directory, provenance, key):
    for name, entry in provenance[key].items():
        want = entry["sha256"] if isinstance(entry, dict) else entry
        path = Path(directory) / Path(name).name
        require(path.is_file() and sha(path) == want, f"Changed or missing artifact: {path}")

def bootstrap_summary(directory, synthetic):
    p = Path(directory)
    if not (p / "summary.json").is_file(): return None, []
    s = load(p / "summary.json"); rs = rows(p / "paired_differences.tsv")
    require(bool(s.get("synthetic", False)) == synthetic, "Bootstrap data-kind mismatch")
    require(s["total_required"] == 500 and s["completed"] == 500 and len(rs) == 500, "Incomplete bootstrap")
    require({int(r["replicate"]) for r in rs} == set(range(1, 501)), "Missing/duplicate bootstrap replicate")
    good = [r for r in rs if truth(r["valid"])]
    require(all(num(r["G1_minus_B0_top20"]) is not None for r in good), "Nonfinite bootstrap difference")
    require(s["valid"] == len(good) and s["invalid"] == 500-len(good), "Bootstrap counts mismatch")
    suppress = len(good) < 450
    require(bool(s["interval_suppressed"]) == suppress, "Bootstrap suppression gate mismatch")
    if suppress: require(s["interval"] is None, "Suppressed interval must be absent")
    else:
        require(len(s["interval"]) == 2 and all(num(x) is not None for x in s["interval"]), "Invalid interval")
        values = sorted(float(r["G1_minus_B0_top20"]) for r in good)
        def q(prob):
            pos = (len(values)-1)*prob; lo = int(pos); hi = math.ceil(pos)
            return values[lo] + (values[hi]-values[lo])*(pos-lo)
        require(all(abs(float(a)-b) < 1e-10 for a,b in zip(s["interval"], (q(.025),q(.975)))), "Interval not valid-draw type7 quantiles")
    return s, rs

def lodo_summary(directory, synthetic):
    p = Path(directory)
    if not (p / "summary.json").is_file(): return None, []
    s = load(p / "summary.json"); rs = rows(p / "fold_metrics.tsv")
    require(bool(s.get("synthetic", False)) == synthetic, "LODO data-kind mismatch")
    n = s["total_folds"]
    require(n == 22 and len(rs) == n*3, "Incomplete discovery LODO: all22 donors and3 methods required")
    seen = {(r["fold"], r["method"]) for r in rs}
    require(len(seen) == len(rs) and all(sum(r["fold"] == str(i) for r in rs) == 3 for i in range(1,n+1)), "Missing/duplicate fold")
    failed = set()
    for r in rs:
        require(r["method"] in METHODS, "Unknown LODO method")
        if not truth(r["valid"]): failed.add(r["fold"]); continue
        overlap = int(r["intersection"])
        require(0 <= overlap <= 20 and abs(float(r["retention"])-overlap/20) < 1e-10, "LODO retention denominator changed")
        require(abs(float(r["jaccard"])-overlap/(40-overlap)) < 1e-10, "LODO Jaccard changed")
    require(s["failed_folds"] == len(failed) and s["valid_folds"] == n-len(failed), "LODO failed-fold mismatch")
    require(not s.get("failed_folds_replaced", False), "Replaced LODO failures")
    return s, rs

def dossier_choices(pool):
    ordered = sorted(pool, key=lambda r: (int(r["B0_rank"]), r["gene_id"]))
    selected = []
    for cat, limit in (("supported_same_direction",3), ("strongly_conflicting",1), ("inconclusive",1)):
        selected.extend([r for r in ordered if r["dossier_category"] == cat][:limit])
    return selected

def validate_counts(summary, selected):
    require(len(selected) == 210, "Selected table must retain all210 method/k rows")
    for method in METHODS:
        for k in (20,50):
            rs = [r for r in selected if r["method"] == method and int(r["k"]) == k]
            require(len(rs) == k and len({r["gene_id"] for r in rs}) == k, "Changed selected-list membership")
            s = summary["methods"][method][str(k)]
            evaluable = sum(truth(r["evaluable"]) for r in rs)
            agree = sum(truth(r["evaluable"]) and truth(r["agree"]) for r in rs)
            require((s["selected"],s["evaluable"],s["agree"],s["unevaluable"]) == (k,evaluable,agree,k-evaluable), "Summary count mismatch")
            expected = (agree/k,(agree+k-evaluable)/k)
            require(all(abs(float(a)-b)<1e-10 for a,b in zip(s["missing_bounds"],expected)), "Missing bounds mismatch")
            for reason in REASONS:
                require(s["missing_reasons"][reason] == sum(r["missing_reason"] == reason for r in rs), "Missing reason count mismatch")
            sig = sum(truth(r["evaluable"]) and truth(r["agree"]) and num(r["external_q"]) is not None and float(r["external_q"])<.05 for r in rs)
            require(s["significant_same_direction"] == sig, "Significant same-direction count mismatch")

def validate_resampling_binding(bootstrap, lodo, bs, ls, de, ev, ranking, freeze_path, track):
    if bs:
        p=Path(bootstrap); df=load(track(p/"draw_freeze.json"))
        require(bs["binding_sha256"] == sha(track(p/"binding.rds")), "Bootstrap binding hash changed")
        require(df["ledger_sha256"] == sha(track(p/"draw_ledger.tsv")) and df["ledger_saved_before_any_fit"], "Bootstrap draw ledger changed")
        binding=df["binding"]["inputs"]
        for name in (de/"model.rds", de/"provenance.json", ev/"bootstrap_interface.tsv", Path(freeze_path)):
            require(binding.get(str(name)) == sha(name), "Bootstrap belongs to a different evaluation/model/freeze")
        for name,digest in binding.items(): require(sha(track(name))==digest,"Bootstrap source changed")
    if ls:
        p=Path(lodo); bp=track(p/"binding.json"); binding=load(bp)
        require(ls["binding_sha256"] == sha(bp), "LODO binding hash changed")
        for name in [Path(freeze_path)]+[ranking/(m+"_top20.tsv") for m in METHODS]:
            require(binding["sha256"].get(str(name)) == sha(name), "LODO belongs to a different discovery baseline")
        for name,digest in binding["sha256"].items(): require(sha(track(name))==digest,"LODO source changed")

def validate_prefix(prefix):
    require(isinstance(prefix,str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,63}",prefix) is not None, "Unsafe output prefix: use1-64 letters/digits/underscore/hyphen, starting with a letter")
    return prefix

def build(args):
    root = Path.cwd()
    ev = Path(args.evaluation); de = Path(args.external_de); ranking = Path(args.ranking)
    output = Path(args.output_dir); receipt = Path(args.receipt)
    prefix = validate_prefix(args.prefix or ("replication_draft" if args.partial else "replication"))
    report_name = prefix + "_results.md"
    inputs = {}
    def track(p):
        p = Path(p); require(p.is_file(), f"Missing source: {p}")
        inputs[str(p)] = sha(p); return p
    if args.synthetic_test:
        for p in (ev,de,ranking,output,receipt):
            require(p.resolve().is_relative_to((root/"work/synthetic").resolve()), "Synthetic paths must remain under work/synthetic")
    summary = load(track(ev/"summary.json")); ep = load(track(ev/"provenance.json"))
    kind = "synthetic_test_only" if args.synthetic_test else "public_cohort"
    require(summary["data_kind"] == kind and ep["data_kind"] == kind, "Evaluation data-kind mismatch")
    verify_outputs(ev,ep,"output_sha256")
    external = load(track(de/"provenance.json"))
    require(external["data_kind"] == kind and bool(external["scientific_effect_run"]) != args.synthetic_test, "External model provenance mismatch")
    require(sha(de/"provenance.json") == ep["external_provenance_sha256"], "Evaluator/model provenance mismatch")
    verify_outputs(de,external,"output_sha256")
    selected = rows(track(ev/"selected_genes.tsv")); pool = rows(track(ev/"dossier_pool.tsv")); choices = rows(track(ev/"dossier_selection.tsv"))
    validate_counts(summary,selected); unique(pool)
    require([r["gene_id"] for r in choices] == [r["gene_id"] for r in dossier_choices(pool)], "Dossier choices violate deterministic rule")
    require(len(pool) == summary["dossier_pool_size"], "Pool count mismatch")
    freeze_path = next((p for p in external["input_sha256"] if Path(p).name == "discovery_freeze.json" or Path(p).name == "discovery_rank_freeze.json"),None)
    require(freeze_path is not None, "Missing discovery checkpoint link")
    freeze = load(track(freeze_path))
    require(sha(freeze_path) == ep["discovery_freeze_sha256"], "Evaluation discovery freeze mismatch")
    ranks = {}
    for method in METHODS:
        fullpath = track(ranking/(method+".tsv"))
        require(freeze["artifacts"].get(str(fullpath)) == sha(fullpath), "Ranking differs from discovery freeze")
        ranks[method] = rows(fullpath)
        for k in (20,50):
            pred = {r["gene_id"]:r for r in ranks[method][:k]}
            rs = [r for r in selected if r["method"] == method and int(r["k"]) == k]
            require(set(pred) == {r["gene_id"] for r in rs}, "Selected genes differ from frozen prefix")
            require(all(float(r["discovery_sign"]) == float(pred[r["gene_id"]]["discovery_sign"]) for r in rs), "Saved prediction sign changed")
    require({r["gene_id"] for r in pool} == {r["gene_id"] for m in ("B0","G1") for r in ranks[m][:20]}, "Changed fixed pool")
    if not args.synthetic_test: require(len(pool) == 40, "Real fixed pool must retain40 genes")
    genetic = {}
    if args.genetic:
        genetic = unique(rows(track(args.genetic)))
        require(set(genetic) == {r["gene_id"] for r in pool}, "Genetic table must cover exact fixed pool")
    elif not args.synthetic_test: raise ValueError("Real report requires complete genetic table")
    if args.genetic_context:
        contexts = rows(track(args.genetic_context))
    else: contexts = []
    for r in pool:
        g = genetic.get(r["gene_id"], {})
        for key,value in g.items():
            if key != "gene_id": r["genetic_"+key] = value
        r["gene_symbol"] = g.get("approved_symbol",r["gene_id"])
    joined = {r["gene_id"]:r for r in pool}; choices = [joined[r["gene_id"]] for r in choices]
    bs,br = bootstrap_summary(args.bootstrap,args.synthetic_test)
    ls,lr = lodo_summary(args.lodo,args.synthetic_test)
    pending = []
    for name,p,s,table in (("bootstrap",Path(args.bootstrap),bs,"paired_differences.tsv"),("discovery LODO",Path(args.lodo),ls,"fold_metrics.tsv")):
        if s is None: pending.append(name)
        else:
            track(p/"summary.json"); track(p/table)
    require(args.partial or not pending, "Resampling pending; use --partial for an explicitly unfinished draft")
    if not args.synthetic_test:
        validate_resampling_binding(args.bootstrap,args.lodo,bs,ls,de,ev,ranking,freeze_path,track)
        track("docs/CURRENT_STATE.md"); track("docs/REMAINING_SCIENTIFIC_WORK.md"); track("reports/discovery_primary_results.md")
    external_rows = rows(track(de/"de_all.tsv")); external_index = unique(external_rows)
    significant = sum(num(r["adj.P.Val"]) is not None and float(r["adj.P.Val"])<.05 for r in external_rows)
    require(len(external_rows) == external["tested_family_size"] and significant == external["significant_BH05"], "External full-family size/significance mismatch")
    primary = summary["primary"]; b0 = summary["methods"]["B0"]["20"]; g1 = summary["methods"]["G1"]["20"]
    missing_bounds = [g1["missing_bounds"][0]-b0["missing_bounds"][1],g1["missing_bounds"][1]-b0["missing_bounds"][0]]
    require(all(abs(a-b)<1e-10 for a,b in zip(missing_bounds,primary["missing_bounds"])), "Primary missing bounds mismatch")
    superiority = bool(b0["evaluable"]>=16 and g1["evaluable"]>=16 and missing_bounds[0]>0 and bs and not bs["interval_suppressed"] and bs["interval"][0]>0)
    lines = ["# 外部复核与研究候选档案 / External replication and research dossiers","",
             "**人工数据测试：以下数字没有生物学含义。 / SYNTHETIC SOFTWARE TEST ONLY.**" if args.synthetic_test else f"**主要结果：B0前20有{b0['agree']}/{b0['evaluable']}/{b0['selected']}同向/可评价/选入，G1有{g1['agree']}/{g1['evaluable']}/{g1['selected']}；外部{significant}/{len(external_rows)}基因达到BH q<0.05。没有证据据此确认治疗靶点或干预方向。**",
             "", f"状态 / Status: {'PARTIAL — pending: '+', '.join(pending) if pending else 'COMPLETE FOR REPORTED ANALYSES'}.",
             "", "## 本轮回答的问题 / Question",
             "", "在增殖期在位子宫内膜基质细胞范围内，事先冻结的表达排序B0、连接度排序B1和网络传播排序G1，能否在另一队列保留表达变化方向？每位供者先汇总计数，以人作为统计单位；不按细胞数增加独立样本。排名、方向、前20/50与随机集合均未利用外部疾病效应重新选择。",
             "", f"外部输入：{external['audit']['ENDO']}位病例、{external['audit']['non_ENDO']}位非病例，共{external['audit']['n_donors']}位供者；{len(external_rows)}个基因进入完整检验家族，{significant}/{len(external_rows)}达到BH q<0.05。BH q同时考虑全部受检基因，不能只在前列名单内重新校正。",
             "", "模型分别在各队列进行库大小/组成归一化和供者级加权拟合，调整冻结的年龄、批次及腺肌症。log2表达差异表示病例相对对照的表达变化；它可能包含细胞状态组成与未测混杂，不能解释为因果效应。对照可能有其他妇科疾病。发现使用估计小数计数、外部使用原始整数核计数，scRNA/snRNA测量差异与相位/激素判定限制仍有效。",
             "", "## 冻结名单复核 / Fixed-list results","",
             "| 方法 / Method | k | 同向 / agree | 可评价 / evaluable | 原始选入 / selected | 同向比例 | 缺失结果范围 | 同向且完整家族q<0.05 |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    plot_counts=[]
    for method in METHODS:
        for k in (20,50):
            s=summary["methods"][method][str(k)]
            lines.append(f"| {method} | {k} | {s['agree']} | {s['evaluable']} | {s['selected']} | {pct(s['concordance'])} | {pct(s['missing_bounds'][0])}–{pct(s['missing_bounds'][1])} | {s['significant_same_direction']} |")
            plot_counts.append(dict(method=method,k=k,agree=s["agree"],opposite=s["opposite"],unevaluable=s["unevaluable"]))
    lines += ["", "外部完整家族最小未校正p="+fmt(min((float(r["P.Value"]) for r in external_rows if num(r["P.Value"]) is not None),default=None))+"，最小BH q="+fmt(min((float(r["adj.P.Val"]) for r in external_rows if num(r["adj.P.Val"]) is not None),default=None))+"。没有达到q阈值并不证明没有表达变化；方向一致与显著证据回答不同问题。"]
    lines += ["","同向比例只回答可评价基因的符号是否一致；微小、不精确的效应也可能同向。缺失范围把未评价基因分别视为全部不同向或全部同向，不补位、不把缺失当反向。前20是主要描述终点，前50是次要终点。",
              "",f"G1−B0前20比例差：{fmt(primary['concordance_difference'])}；缺失结果范围[{fmt(missing_bounds[0])}, {fmt(missing_bounds[1])}]。两方法均须至少16/20可评价；若范围包含零或可逆转，则不能主张优势。当前"+("在覆盖、缺失范围及条件重采样区间下有描述性正向差异；仍不代表目标准确率、临床效用或治疗优越性。" if superiority else "没有满足全部预设条件的G1优于B0证据。"),
              "","| 方法 | absent ID | 表达过滤失败 | 拟合失败 | 发现零方向 | 外部零方向 |","|---|---:|---:|---:|---:|---:|"]
    for m in METHODS: lines.append("| "+m+" | "+" | ".join(str(summary["methods"][m]["20"]["missing_reasons"][r]) for r in REASONS)+" |")
    overlap=len({r["gene_id"] for r in ranks["B1"][:20]} & {r["gene_id"] for r in ranks["G1"][:20]})
    lines += ["",f"B1与G1前20重合{overlap}/20。网络传播可能偏向高连接基因；B1用于识别这种倾向。STRING边表示功能关联，不能视为因果调控或物理结合证据。"]
    null=load(track(ev/"null_summary.json"))
    lines += ["","## 匹配随机集合 / Matched sets",""]
    for m in METHODS:
        s=null[m]
        lines.append(f"- {m}: 状态{s['status']}；充分覆盖{s.get('adequate_sets','NA')}/{s.get('total_sets','NA')}组；观察百分位{fmt(s.get('observed_descriptive_percentile'))}；原观察名单覆盖{summary['methods'][m]['20']['evaluable']}/20。")
    lines += ["","每种可估方法使用冻结的200组表达量/连接度匹配集合，方向来自各基因自己的发现效应。百分位是描述性位置；这些集合共享基因与供者，不是200个独立生物学重复，也不是校准的p值。",
              "","## 供者敏感性 / Donor sensitivity",""]
    if bs:
        lines += [f"外部配对重采样完成{bs['completed']}/500次；有效{bs['valid']}/500，失败{bs['invalid']}/500。所有方法共享相同供者抽样，固定发现名单/方向和原始可评价分母。失败不替换。",
                  "失败超过50/500，按预设规则抑制区间；不从成功子集宣称稳定性。" if bs["interval_suppressed"] else f"G1−B0前20条件区间[{fmt(bs['interval'][0])}, {fmt(bs['interval'][1])}]（有效抽样的type7百分位）；"+("存在失败，区间以成功抽样为条件。" if bs["invalid"] else "没有失败抽样。"),
                  "该区间只描述外部供者抽样变化，条件于已经选定的发现名单；不包含发现选择、不确定注释或未来实验的不确定性。",
                  f"失败原因 / Failure reasons: {json.dumps(bs['failure_reasons'],ensure_ascii=False)}"]
    else: lines.append("PENDING：外部500次配对供者重采样尚未完成，不提供区间。")
    if ls:
        lines += ["",f"发现逐人留出：有效{ls['valid_folds']}/{ls['total_folds']}轮，失败{ls['failed_folds']}/{ls['total_folds']}轮；每次移除一人后重新过滤、拟合和排名，失败轮保留。","",
                  "| 方法 | 有效轮 / 全部轮 | 前20保留比例中位数 [范围] | Jaccard中位数 [范围] |","|---|---:|---:|---:|"]
        for m in METHODS:
            good=[r for r in lr if r["method"]==m and truth(r["valid"])]
            def stat(key):
                vals=[float(r[key]) for r in good]
                return f"{fmt(statistics.median(vals))} [{fmt(min(vals))}, {fmt(max(vals))}]" if vals else "NA"
            lines.append(f"| {m} | {len(good)}/{ls['total_folds']} | {stat('retention')} | {stat('jaccard')} |")
        lines += ["","保留比例=交集/原始20；Jaccard=交集/两名单并集。表中的中位数只描述有效轮，不能遮盖失败或当作独立重复置信区间。"]
    else: lines.append("PENDING：发现22位供者逐人留出尚未完成；不提供稳定性结论。")
    lines += ["","## 完整候选池与确定性档案 / Pool and deterministic dossiers","",
              f"冻结的B0/G1前20并集共{len(pool)}个基因，全部保留在附表。最多3个同向且完整家族q<0.05、1个精确反向、1个缺失/不确定案例；每类按完整B0名次和稳定ID选择。精确反向指未经多重校正的95%区间不跨零，q值另列。没有合格阳性时不强行产生阳性名单。",
              "",f"本次确定性选择{len(choices)}个档案，其中同向且q<0.05的档案{sum(r['dossier_category']=='supported_same_direction' for r in choices)}个。",
              "","遗传证据来自单一预设2023论文及其公开映射表，同一论文多个注释不是独立确认。缺少映射不等于阴性遗传检验；病例表达变化方向与遗传QTL方向也不能直接等同。没有证据建立独立个体非重叠、增殖期基质细胞特异性或干预方向。"]
    for r in choices:
        g=genetic.get(r["gene_id"],{})
        category=r["dossier_category"]; rec={"supported_same_direction":"follow up / 研究随访","strongly_conflicting":"deprioritize / 在本表达假设下暂缓优先","inconclusive":"hold / 保留待证"}[category]
        lines += ["",f"### {r['gene_symbol']} ({r['gene_id']})",
                  f"- 表达 / Expression: {category}；发现log2差异{fmt(unique(ranks['B0'])[r['gene_id']]['logFC'])}；外部log2差异{fmt(r['external_logFC'])}，95%未校正区间[{fmt(r['external_CI_low'])}, {fmt(r['external_CI_high'])}]，完整家族q={fmt(r['external_q'])}；缺失原因={r['missing_reason'] or 'none'}。",
                  f"- 遗传 / Genetics: {g.get('genetic_status','unavailable')}；映射类型={g.get('mapping_categories','unavailable') or 'unavailable'}；底层发表数={g.get('publication_count','unavailable')}。"+
                  ("SKIDA1的公开证据限于位置注释和共享疾病/疼痛分析中的bulk内膜QTL；不构成共定位、孟德尔随机化或因果治疗证据。" if r['gene_symbol']=="SKIDA1" else ""),
                  "- 特异性、矛盾与干预 / Specificity, contradictions and intervention: 单个基因的细胞状态特异性未独立建立；"+("精确反向与原发现方向冲突，未经多重校正区间不能等同显著q。" if category=="strongly_conflicting" else "当前证据范围内未建立遗传机制矛盾；这不是证明不存在矛盾。")+"干预方向未知，关联不能指导增加或抑制基因。",
                  f"- 研究建议 / Recommendation: {rec}；"+("可作为组织/状态适用性与扰动实验的候选，尚未完成这些实验。" if category=="supported_same_direction" else "先解决表达不一致、精度或缺失证据；不作药物或临床建议。")]
    lines += ["","## 仍未解决与执行边界 / Limits and execution","",
              "临床敏感性、细胞状态组成及相位/周期日来源说明属于另外的进行中交付，见[剩余科学工作审计](../docs/REMAINING_SCIENTIFIC_WORK.md)；本报告只完成已列出的主分析复核和供者敏感性，不宣称全研究完成。30/50/100细胞阈值保留相同的22+22位供者、相同库及设计，因此复用主结果，不进行三个重复拟合，也不将阈值等价称为三份独立验证。44/44位供者年龄完整，没有因年龄缺失排除人；原方案的缺失年龄条件不触发额外无年龄调整拟合，但年龄完整不证明不存在其他混杂。来源可分离不等于已经证明所有历史供者彼此独立。遗传注释只覆盖预设论文；不能称为2026全面证据综述。",
              "","源数据限制：发现复用OSF作者二级汇总，保留Alevin估计小数计数；完整源RDS没有重算。HECA52,744个核的原始计数与保存H5总计数差9,607/360,765,431（0.00266%），其中1,359可由全局缺失特征解释、8,248尚未直接归因，未修正以强行吻合。继承的相位/细胞标签和标记表达支持群体合理性，不能证明纯度或独立重注释。来源分离不等于基因型证实个体非重叠。",
              "","本步骤只读取已有模型/评价结果并生成图文，没有重新拟合、重排或下载。实际环境为SRV01原生运行；Docker与AWS尚未执行，不能称为容器或云验证。测试、哈希和生成图形成功只证明工程可复现，不证明生物学、治疗或临床有效。",
              "","YELLOW — 理解尚未评估 / Understanding not yet assessed：需能解释供者分母、缺失范围、完整家族q、固定名单重采样与发现逐人留出的不同、不显著与无效应的区别、遗传映射与干预因果的区别。",
              "","本报告构建无需用户作新决定。仍待步骤保持待完成，不扩大已授权分析范围。",
              "","## 可追溯输入与图形 / Traceable inputs and figures",""]
    for p in inputs: lines.append(f"- `{p}`")
    for name in ("counts","scatter","lodo") + (("bootstrap",) if bs and not bs["interval_suppressed"] else ()):
        if name=="lodo" and not ls: continue
        lines += [f"- [{name} PDF](figures/{prefix}_{name}.pdf) / [PNG](figures/{prefix}_{name}.png)"]
    report="\n".join(lines)+"\n"
    destinations=[output/report_name,output/(prefix+"_dossier_pool.tsv"),output/(prefix+"_dossiers.tsv"),receipt]
    fig_names=["counts","scatter"]+(["lodo"] if ls else [])+(["bootstrap"] if bs and not bs["interval_suppressed"] else [])
    destinations += [output/"figures"/(prefix+"_"+n+"."+ext) for n in fig_names for ext in ("pdf","png")]
    require(all(not p.exists() for p in destinations), "Refuse overwrite: report/receipt/figure already exists")
    require(all(sha(p)==digest for p,digest in inputs.items()), "Input changed while report was assembled")
    output.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".replication-report-",dir=output) as tmp:
        stage=Path(tmp); (stage/"figures").mkdir()
        (stage/report_name).write_text(report)
        fields=list(pool[0].keys())
        write_tsv(stage/(prefix+"_dossier_pool.tsv"),pool,fields); write_tsv(stage/(prefix+"_dossiers.tsv"),choices,fields)
        write_tsv(stage/"counts.tsv",plot_counts,("method","k","agree","opposite","unevaluable"))
        discovery=unique(ranks["B0"]); memberships={r["gene_id"]:{m:False for m in METHODS} for r in pool}
        for m in METHODS:
            for r in ranks[m][:20]: memberships.setdefault(r["gene_id"],{m:False for m in METHODS})[m]=True
        scatter=[]
        for gid,r in discovery.items():
            e=external_index.get(gid)
            if e and num(e["logFC"]) is not None and truth(e["fit_estimable"]):
                scatter.append(dict(gene_id=gid,discovery_logFC=r["logFC"],external_logFC=e["logFC"],selected=gid in memberships))
        write_tsv(stage/"scatter.tsv",scatter,("gene_id","discovery_logFC","external_logFC","selected"))
        if ls: write_tsv(stage/"lodo.tsv",lr,lr[0].keys())
        if bs and not bs["interval_suppressed"]:
            write_tsv(stage/"bootstrap.tsv",br,br[0].keys())
            write_tsv(stage/"bootstrap_interval.tsv",[dict(low=bs["interval"][0],high=bs["interval"][1],valid=bs["valid"],failed=bs["invalid"])],("low","high","valid","failed"))
        plot_script=track(args.plot_script)
        run=subprocess.run([args.rscript,str(plot_script),str(stage),str(stage/"figures"),prefix,"SYNTHETIC TEST ONLY" if args.synthetic_test else "HECA donor-level replication"],text=True,capture_output=True)
        require(run.returncode==0,f"R figures failed: {run.stderr}")
        generated=[p for p in stage.rglob("*") if p.is_file() and (p.parent==stage/"figures" or p.name in {report_name,prefix+"_dossier_pool.tsv",prefix+"_dossiers.tsv"})]
        require(all((stage/"figures"/(prefix+"_"+n+"."+e)).stat().st_size>1000 for n in fig_names for e in ("pdf","png")), "Missing/empty figures")
        rendering_revision=None
        if args.supersedes_receipt:
            previous=load(track(args.supersedes_receipt))
            require(previous["data_kind"]==kind and previous["bootstrap"]==bs and previous["lodo"]==ls, "Prior rendering has different scientific results")
            checked=0
            for name,digest in previous["input_sha256"].items():
                if Path(name).resolve()==Path(args.plot_script).resolve(): continue
                require(inputs.get(name)==digest, "Prior rendering has different scientific input: "+name)
                checked+=1
            require(checked>0, "Prior receipt has no comparable inputs")
            rendering_revision=dict(supersedes_receipt=str(args.supersedes_receipt),supersedes_receipt_sha256=sha(args.supersedes_receipt),reason="layout-only title/margin correction; v1 preserved",scientific_inputs_identical=True,bootstrap_and_lodo_identical=True,shared_input_count=checked)
        receipt_data=dict(schema_version=1,rendering_revision=rendering_revision,completed_utc=datetime.now(timezone.utc).isoformat(),data_kind=kind,partial=bool(pending),pending=pending,scientific_validation=False,therapeutic_validation=False,container_execution=False,aws_execution=False,execution_host=platform.node(),python=platform.python_version(),input_sha256=inputs,implementation_sha256=sha(__file__),bootstrap=bs,lodo=ls,superiority_descriptive_guards_pass=superiority,plot_stdout=run.stdout,plot_stderr=run.stderr,output_sha256={str(output/p.relative_to(stage)):sha(p) for p in generated})
        for p in generated:
            target=output/p.relative_to(stage); target.parent.mkdir(parents=True,exist_ok=True)
            require(not target.exists(),"Concurrent output exists"); os.link(p,target)
        receipt.parent.mkdir(parents=True,exist_ok=True)
        with open(receipt,"x") as f: json.dump(receipt_data,f,indent=2,ensure_ascii=False); f.write("\n")
    return dict(report=str(output/report_name),receipt=str(receipt),pending=pending,figures=len(fig_names)*2)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evaluation",default="work/statistical/heca_primary/evaluation")
    p.add_argument("--external-de",default="work/statistical/heca_primary/de")
    p.add_argument("--ranking",default="work/statistical/discovery_primary/ranking")
    p.add_argument("--genetic",default="reports/genetic_evidence_summary.tsv")
    p.add_argument("--genetic-context",default="reports/genetic_evidence_context.tsv")
    p.add_argument("--bootstrap",default="work/statistical/bootstrap")
    p.add_argument("--lodo",default="work/statistical/discovery_lodo")
    p.add_argument("--output-dir",default="reports")
    p.add_argument("--receipt",default=None)
    p.add_argument("--prefix",default=None)
    p.add_argument("--supersedes-receipt",default=None,help="Previous immutable rendering receipt; require identical scientific inputs")
    p.add_argument("--plot-script",default="scripts/plot_replication.R")
    p.add_argument("--rscript",default="Rscript")
    p.add_argument("--partial",action="store_true")
    p.add_argument("--synthetic-test",action="store_true")
    args=p.parse_args()
    if args.prefix is None:
        args.prefix="replication_draft" if args.partial else "replication"
    try: validate_prefix(args.prefix)
    except ValueError as e: p.error(str(e))
    if args.receipt is None:
        args.receipt="receipts/"+args.prefix+"_report.json"
    try: print(json.dumps(build(args),ensure_ascii=False))
    except (ValueError,KeyError,FileNotFoundError) as e: p.exit(1,str(e)+"\n")
if __name__=="__main__": main()

