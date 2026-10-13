# 外部复核与研究候选档案 / External replication and research dossiers

**主要结果：B0前20有14/20/20同向/可评价/选入，G1有6/19/20；外部0/17643基因达到BH q<0.05。没有证据据此确认治疗靶点或干预方向。**

状态 / Status: COMPLETE FOR REPORTED ANALYSES.

## 本轮回答的问题 / Question

在增殖期在位子宫内膜基质细胞范围内，事先冻结的表达排序B0、连接度排序B1和网络传播排序G1，能否在另一队列保留表达变化方向？每位供者先汇总计数，以人作为统计单位；不按细胞数增加独立样本。排名、方向、前20/50与随机集合均未利用外部疾病效应重新选择。

外部输入：15位病例、7位非病例，共22位供者；17643个基因进入完整检验家族，0/17643达到BH q<0.05。BH q同时考虑全部受检基因，不能只在前列名单内重新校正。

模型分别在各队列进行库大小/组成归一化和供者级加权拟合，调整冻结的年龄、批次及腺肌症。log2表达差异表示病例相对对照的表达变化；它可能包含细胞状态组成与未测混杂，不能解释为因果效应。对照可能有其他妇科疾病。发现使用估计小数计数、外部使用原始整数核计数，scRNA/snRNA测量差异与相位/激素判定限制仍有效。

## 冻结名单复核 / Fixed-list results

| 方法 / Method | k | 同向 / agree | 可评价 / evaluable | 原始选入 / selected | 同向比例 | 缺失结果范围 | 同向且完整家族q<0.05 |
|---|---:|---:|---:|---:|---:|---:|---:|
| B0 | 20 | 14 | 20 | 20 | 70.0% | 70.0%–70.0% | 0 |
| B0 | 50 | 36 | 49 | 50 | 73.5% | 72.0%–74.0% | 0 |
| B1 | 20 | 3 | 13 | 20 | 23.1% | 15.0%–50.0% | 0 |
| B1 | 50 | 16 | 35 | 50 | 45.7% | 32.0%–62.0% | 0 |
| G1 | 20 | 6 | 19 | 20 | 31.6% | 30.0%–35.0% | 0 |
| G1 | 50 | 26 | 47 | 50 | 55.3% | 52.0%–58.0% | 0 |

外部完整家族最小未校正p=6.993e-05，最小BH q=0.9699。没有达到q阈值并不证明没有表达变化；方向一致与显著证据回答不同问题。

同向比例只回答可评价基因的符号是否一致；微小、不精确的效应也可能同向。缺失范围把未评价基因分别视为全部不同向或全部同向，不补位、不把缺失当反向。前20是主要描述终点，前50是次要终点。

G1−B0前20比例差：-0.3842；缺失结果范围[-0.4, -0.35]。两方法均须至少16/20可评价；若范围包含零或可逆转，则不能主张优势。当前没有满足全部预设条件的G1优于B0证据。

| 方法 | absent ID | 表达过滤失败 | 拟合失败 | 发现零方向 | 外部零方向 |
|---|---:|---:|---:|---:|---:|
| B0 | 0 | 0 | 0 | 0 | 0 |
| B1 | 0 | 7 | 0 | 0 | 0 |
| G1 | 0 | 1 | 0 | 0 | 0 |

B1与G1前20重合9/20。网络传播可能偏向高连接基因；B1用于识别这种倾向。STRING边表示功能关联，不能视为因果调控或物理结合证据。

## 匹配随机集合 / Matched sets

- B0: 状态ESTIMABLE；充分覆盖200/200组；观察百分位96.5；原观察名单覆盖20/20。
- B1: 状态ESTIMABLE；充分覆盖200/200组；观察百分位0.5；原观察名单覆盖13/20。
- G1: 状态ESTIMABLE；充分覆盖200/200组；观察百分位0.5；原观察名单覆盖19/20。

每种可估方法使用冻结的200组表达量/连接度匹配集合，方向来自各基因自己的发现效应。百分位是描述性位置；这些集合共享基因与供者，不是200个独立生物学重复，也不是校准的p值。

## 供者敏感性 / Donor sensitivity

外部配对重采样完成500/500次；有效487/500，失败13/500。所有方法共享相同供者抽样，固定发现名单/方向和原始可评价分母。失败不替换。
G1−B0前20条件区间[-0.3895, 0.08684]（有效抽样的type7百分位）；存在失败，区间以成功抽样为条件。
该区间只描述外部供者抽样变化，条件于已经选定的发现名单；不包含发现选择、不确定注释或未来实验的不确定性。
失败原因 / Failure reasons: {"design not full rank; nonconstant nuisance cannot be dropped": 13}

发现逐人留出：有效22/22轮，失败0/22轮；每次移除一人后重新过滤、拟合和排名，失败轮保留。

| 方法 | 有效轮 / 全部轮 | 前20保留比例中位数 [范围] | Jaccard中位数 [范围] |
|---|---:|---:|---:|
| B0 | 22/22 | 0.725 [0.3, 1] | 0.5692 [0.1765, 1] |
| B1 | 22/22 | 1 [1, 1] | 1 [1, 1] |
| G1 | 22/22 | 0.975 [0.8, 1] | 0.9524 [0.6667, 1] |

保留比例=交集/原始20；Jaccard=交集/两名单并集。表中的中位数只描述有效轮，不能遮盖失败或当作独立重复置信区间。

## 完整候选池与确定性档案 / Pool and deterministic dossiers

冻结的B0/G1前20并集共40个基因，全部保留在附表。最多3个同向且完整家族q<0.05、1个精确反向、1个缺失/不确定案例；每类按完整B0名次和稳定ID选择。精确反向指未经多重校正的95%区间不跨零，q值另列。没有合格阳性时不强行产生阳性名单。

本次确定性选择1个档案，其中同向且q<0.05的档案0个。

遗传证据来自单一预设2023论文及其公开映射表，同一论文多个注释不是独立确认。缺少映射不等于阴性遗传检验；病例表达变化方向与遗传QTL方向也不能直接等同。没有证据建立独立个体非重叠、增殖期基质细胞特异性或干预方向。

### PARP9 (ENSG00000138496)
- 表达 / Expression: inconclusive；发现log2差异0.4659；外部log2差异-0.1415，95%未校正区间[-0.4126, 0.1296]，完整家族q=0.9699；缺失原因=none。
- 遗传 / Genetics: no_available_mapping_in_audited_source；映射类型=unavailable；底层发表数=0。
- 特异性、矛盾与干预 / Specificity, contradictions and intervention: 单个基因的细胞状态特异性未独立建立；当前证据范围内未建立遗传机制矛盾；这不是证明不存在矛盾。干预方向未知，关联不能指导增加或抑制基因。
- 研究建议 / Recommendation: hold / 保留待证；先解决表达不一致、精度或缺失证据；不作药物或临床建议。

## 仍未解决与执行边界 / Limits and execution

临床敏感性、细胞状态组成及相位/周期日来源说明属于另外的进行中交付，见[剩余科学工作审计](../docs/REMAINING_SCIENTIFIC_WORK.md)；本报告只完成已列出的主分析复核和供者敏感性，不宣称全研究完成。30/50/100细胞阈值保留相同的22+22位供者、相同库及设计，因此复用主结果，不进行三个重复拟合，也不将阈值等价称为三份独立验证。44/44位供者年龄完整，没有因年龄缺失排除人；原方案的缺失年龄条件不触发额外无年龄调整拟合，但年龄完整不证明不存在其他混杂。来源可分离不等于已经证明所有历史供者彼此独立。遗传注释只覆盖预设论文；不能称为2026全面证据综述。

源数据限制：发现复用OSF作者二级汇总，保留Alevin估计小数计数；完整源RDS没有重算。HECA52,744个核的原始计数与保存H5总计数差9,607/360,765,431（0.00266%），其中1,359可由全局缺失特征解释、8,248尚未直接归因，未修正以强行吻合。继承的相位/细胞标签和标记表达支持群体合理性，不能证明纯度或独立重注释。来源分离不等于基因型证实个体非重叠。

本步骤只读取已有模型/评价结果并生成图文，没有重新拟合、重排或下载。实际环境为SRV01原生运行；Docker与AWS尚未执行，不能称为容器或云验证。测试、哈希和生成图形成功只证明工程可复现，不证明生物学、治疗或临床有效。

YELLOW — 理解尚未评估 / Understanding not yet assessed：需能解释供者分母、缺失范围、完整家族q、固定名单重采样与发现逐人留出的不同、不显著与无效应的区别、遗传映射与干预因果的区别。

本报告构建无需用户作新决定。仍待步骤保持待完成，不扩大已授权分析范围。

## 可追溯输入与图形 / Traceable inputs and figures

- `work/statistical/heca_primary/evaluation/summary.json`
- `work/statistical/heca_primary/evaluation/provenance.json`
- `work/statistical/heca_primary/de/provenance.json`
- `work/statistical/heca_primary/evaluation/selected_genes.tsv`
- `work/statistical/heca_primary/evaluation/dossier_pool.tsv`
- `work/statistical/heca_primary/evaluation/dossier_selection.tsv`
- `receipts/discovery_rank_freeze.json`
- `work/statistical/discovery_primary/ranking/B0.tsv`
- `work/statistical/discovery_primary/ranking/B1.tsv`
- `work/statistical/discovery_primary/ranking/G1.tsv`
- `reports/genetic_evidence_summary.tsv`
- `reports/genetic_evidence_context.tsv`
- `work/statistical/bootstrap/summary.json`
- `work/statistical/bootstrap/paired_differences.tsv`
- `work/statistical/discovery_lodo/summary.json`
- `work/statistical/discovery_lodo/fold_metrics.tsv`
- `work/statistical/bootstrap/draw_freeze.json`
- `work/statistical/bootstrap/binding.rds`
- `work/statistical/bootstrap/draw_ledger.tsv`
- `data/processed/heca/counts_stable.tsv.gz`
- `manifests/final_heca_clinical_selected.tsv`
- `config/statistical_analysis.json`
- `scripts/resample_donor_models.R`
- `scripts/donor_de_lib.R`
- `scripts/statistical_runtime_env.sh`
- `requirements/statistical-runtime.explicit.txt`
- `receipts/statistical_pre_effect_freeze.json`
- `work/statistical/heca_primary/de/model.rds`
- `work/statistical/heca_primary/evaluation/bootstrap_interface.tsv`
- `work/statistical/discovery_lodo/binding.json`
- `data/processed/discovery/counts_stable.tsv.gz`
- `manifests/final_discovery_clinical_selected.tsv`
- `data/reference/string/direct_evidence_gene_edges.tsv.gz`
- `config/ranking_algorithm.json`
- `/home/analyst/te/scripts/run_discovery_lodo.py`
- `/home/analyst/te/scripts/resample_donor_models.R`
- `/home/analyst/te/scripts/donor_de_lib.R`
- `/home/analyst/te/scripts/rank_targets.py`
- `/home/analyst/te/scripts/statistical_runtime_env.sh`
- `/home/analyst/te/requirements/statistical-runtime.explicit.txt`
- `work/statistical/discovery_primary/ranking/B0_top20.tsv`
- `work/statistical/discovery_primary/ranking/B1_top20.tsv`
- `work/statistical/discovery_primary/ranking/G1_top20.tsv`
- `docs/CURRENT_STATE.md`
- `docs/REMAINING_SCIENTIFIC_WORK.md`
- `reports/discovery_primary_results.md`
- `work/statistical/heca_primary/de/de_all.tsv`
- `work/statistical/heca_primary/evaluation/null_summary.json`
- [counts PDF](figures/replication_counts.pdf) / [PNG](figures/replication_counts.png)
- [scatter PDF](figures/replication_scatter.pdf) / [PNG](figures/replication_scatter.png)
- [lodo PDF](figures/replication_lodo.pdf) / [PNG](figures/replication_lodo.png)
- [bootstrap PDF](figures/replication_bootstrap.pdf) / [PNG](figures/replication_bootstrap.png)
