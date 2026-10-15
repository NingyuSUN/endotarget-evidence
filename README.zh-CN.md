# EndoTarget Evidence

**以供者级转录组和网络分析评估子宫内膜异位症候选靶点的跨研究证据。**

[English](README.md) · [简体中文](README.zh-CN.md)

[![数据：公开研究](https://img.shields.io/badge/data-public_studies-466A83)](docs/INPUT_PROVENANCE.md)
[![Nextflow：25.10.2](https://img.shields.io/badge/Nextflow-25.10.2-287C72)](docs/LOCAL_EXECUTION.md)
[![Docker：本地复现已核验](https://img.shields.io/badge/Docker-local_replay_verified-287C72)](receipts/local_primary_container_replay.json)

生物网络能否帮助筛出在另一项研究中仍保留疾病表达关联的基因？本项目比较表达排序与两种网络基线，在另一子宫内膜队列中评价冻结的候选名单，再单独整理遗传证据，形成可追溯的研究判断。

**主要结论：未形成得到统计显著跨研究支持的靶点名单，也未建立网络传播优于表达排序的证据。研究建议是保留候选池，等待更强证据。**

| 项目概览 | 范围 |
|---|---|
| 生物学对象 | 增殖期在位子宫内膜基质成纤维细胞 |
| 数据 | 发现队列 22 名供者、150,788 个细胞；外部队列 22 名供者、52,744 个细胞核 |
| 方法比较 | 在共同的 11,749 个基因上比较表达、网络连接度和固定网络传播 |
| 评价 | 外部表达方向一致性、供者重采样、临床敏感性分析 |
| 工程复现 | Docker 中重新执行 5 个主要 Nextflow 步骤，34 张科学表和 8 份结果元数据一致 |

**结果入口：** [两页英文研究简报](reports/executive_brief_en.pdf) · [方法卡](docs/METHOD_CARD.md) · [完整结果与图表](reports/replication_v2_results.md)

## 研究设计

先按供者汇总计数，再进行统计建模；同一个人的大量细胞不会被当成大量独立患者。两个队列分别使用 edgeR 和 limma-voom 建模，调整年龄、技术批次和腺肌症。

```mermaid
flowchart LR
    A[公开计数与临床信息] --> B[供者筛选与计数汇总]
    B --> C[发现队列效应与三种排序]
    C --> D[冻结名单、方向与参数]
    D --> E[外部评价与敏感性分析]
    E --> F[证据报告与研究判断]
    G[单独整理的公开遗传证据] --> F
```

| 队列 | 病例 / 对照 | 保留的细胞或核 | 受检基因 | BH q < 0.05 的基因 |
|---|---:|---:|---:|---:|
| 发现：GSE266265 | 12 / 10 | 150,788 个细胞 | 17,733 | 0 |
| 外部：HECA，E-MTAB-14039 | 15 / 7 | 52,744 个细胞核 | 17,643 | 0 |

统计样本量是**每个队列 22 名供者**。对照为非内异症手术供者，可能有其他疾病。两研究招募来源分离，但没有基因型数据确认所有个体完全不重叠。[数据来源与纳入规则](docs/INPUT_PROVENANCE.md)。

## 基线比较结果

三种方法使用相同的 11,749 个基因和 123,313 条 STRING 功能关联边。名单、发现队列表达方向和参数均在计算外部疾病效应前冻结。

| 方法 | 排序依据 | 同向 / 可评价 / 原定入选 |
|---|---|---:|
| **表达排序 B0** | 调整后疾病统计量的绝对值 | **14 / 20 / 20** |
| **连接度排序 B1** | 加权网络连接度，作为拓扑对照 | **3 / 13 / 20** |
| **网络传播 G1** | 以发现表达证据为起点的带重启随机游走 | **6 / 19 / 20** |

G1 是固定图算法，没有模型训练或参数搜索。无法在外部评价的基因保留为缺失，不从后续名次补位。B1 未达到预定的至少 16/20 可评价门槛。两种网络名单的同向比例都低于各自 200 组表达量/连接度匹配参考集合（第 0.5 百分位，B0 为第 96.5 百分位）：这是描述性的“低于参考”模式，而不只是没有改进；该百分位不是校准的 p 值。

![固定前20与前50名单在外部队列中的同向、反向和无法评价基因数](reports/figures/replication_v2_counts.png)

**如何理解：** 同向只是表达变化正负号一致，不能作为靶点准确率。没有入选基因同时满足同向且外部 BH q < 0.05；BH 校正使用完整受检基因集合，控制多重比较造成的假阳性。不显著也不证明没有疾病效应。

固定样本中，G1 减 B0 的同向比例相差 **−38.4 个百分点**。500 次配对外部供者重采样中，487 次有效、13 次未通过模型设计检查；95% 条件区间为 **−38.9 至 +8.7 个百分点**，包含零，不能确认一般性方法优势。重采样噪声会把两种方法的同向比例都拉向 50%，因此该区间并不以观察差值为中心（见[英文结果摘要的解读说明](reports/results_summary_en.md#notes-for-interpretation)）。发现队列逐一排除供者的 22/22 轮分析均已完成。[不确定性与敏感性结果](reports/replication_v2_results.md)。

遗传层保留全部 40 个冻结候选基因：所选来源中，SKIDA1 有有限映射，另 39 个缺少可用映射。档案规则产生 1 个 PARP9 不确定案例，没有得到支持的靶点档案。[遗传证据报告](reports/genetic_evidence_report.md)。

## 方法与工程实现

| 组成 | 实现 | 核验入口 |
|---|---|---|
| 供者级分析 | 临床纳入审核、计数来源、协变量调整与明确排除记录 | [方法卡](docs/METHOD_CARD.md)、[供者筛选图](reports/figures/donor_flow_v1.svg) |
| 受控图方法比较 | 共同基因集合、固定参数、外部评价与缺失结果范围 | [排序代码](scripts/rank_targets.py)、[算法参数](config/ranking_algorithm.json) |
| 不确定性分析 | 配对供者重采样、逐供者排除和描述性匹配参考集合 | [重采样方法](docs/RESAMPLING.md)、[运行记录](receipts/real_resampling_execution.json) |
| 可复现流程 | Nextflow DSL2、锁定 R/Python 依赖、Docker、文件哈希和缓存检查 | [发现流程](main.nf)、[外部流程](replication.nf)、[复现说明](docs/LOCAL_EXECUTION.md) |

本地容器复现完成全部 5 个主要步骤，离散选择完全一致，连续数值最大差 **1.10 × 10⁻¹³**，小于事先规定的绝对误差 10⁻⁸。合成数据测试还完成 5 个首次执行步骤和 5 个缓存复用步骤。这些是对相同数据的工程核验，不能替代生物学复核。[主要结果复现记录](receipts/local_primary_container_replay.json) · [缓存检查](receipts/local_container_nextflow.json)。

## 查看与复现

报告、图表、代码和小型证据文件可以直接在 GitHub 阅读，也可以克隆到本地：

```bash
git clone https://github.com/NingyuSUN/endotarget-evidence.git
cd endotarget-evidence
```

**仅克隆仓库并不具备完整重跑条件。** 大型输入、中间结果、运行环境和缓存不在 Git 中；完整复现还需恢复记录引用的文件，并匹配对应环境和路径。

| 目的 | 入口 |
|---|---|
| 阅读研究结论 | [英文简报](reports/executive_brief_en.pdf)、[英文结果摘要](reports/results_summary_en.md)、[中文生物学报告](reports/biology_report_zh.md) |
| 核对方法与图表 | [方法卡](docs/METHOD_CARD.md) |
| 复现指定的主要结果 | [前提、命令和检查标准](docs/LOCAL_EXECUTION.md) |
| 查看云端部署方案 | [AWS 说明：未执行](docs/AWS_DEPLOYMENT.md) |
| 查找方法、结果和来源 | [文档目录](docs/README.md) |

## 适用范围与限制

- 这是回顾性公开数据分析，外部效应计算前有冻结记录，但不是预注册或盲法研究。作者标签、分期判定、残余混杂、细胞构成及细胞与细胞核测量差异仍限制解释。
- 发现数据采用作者二级汇总的估计计数，没有重算完整原 RDS。HECA 原始计数与 H5 差异中仍有 8,248 个计数未逐特征归因，详见[方法卡](docs/METHOD_CARD.md)。
- 遗传证据限定于一项 2023 年研究。缺少映射不等于阴性证据，表达关联不等于因果关系，也不能确定激活或抑制策略。
- 本项目用于早期靶点证据评估，不生成药物分子或抗体，没有训练 GNN 或基础模型，也没有验证治疗效果。AWS 部署尚未执行。

## 许可

代码采用 [MIT 许可证](LICENSE)。本项目撰写的文档、报告和图表采用 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.zh-hans)。第三方数据、参考资源和出版物仍适用其各自条款，见[来源登记](docs/DATA_SOURCES.md)和[输入来源说明](docs/INPUT_PROVENANCE.md)。

维护者：[Ningyu Sun](https://github.com/NingyuSUN)。实现与审查使用了 AI 辅助。原始研究及来源归属见[输入来源说明](docs/INPUT_PROVENANCE.md)和[遗传证据方法](docs/GENETIC_EVIDENCE_METHOD.md)。

历史运行编号与当前公开提交的对应关系见[来源追溯说明](docs/PUBLICATION_HISTORY.md)。
