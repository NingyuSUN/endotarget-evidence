# EndoTarget Evidence

这个项目用公开的子宫内膜数据，比较表达排序、连接度排序和网络传播能否在另一研究中保留基因表达变化方向，并独立整理遗传证据。原始科学分析在SRV01完成；现在按用户要求转到本地 G:\projects\te，进行容器验证和限定的主要结果复现。服务器原件保留；GitHub项目仓库：[NingyuSUN/endotarget-evidence](https://github.com/NingyuSUN/endotarget-evidence)。

**已经得到完整的主要科学结果：没有产生满足预设显著跨队列支持的靶点名单，也没有证据证明网络传播优于表达排序。** 失败和不确定结果均保留。

快速阅读：[两页英文简报](reports/executive_brief_en.pdf)、[方法卡与五类图](docs/METHOD_CARD.md)、[当前本地复现说明](docs/LOCAL_EXECUTION.md)。

从[面向生物学家的综合报告](reports/biology_report_zh.md)开始，或查看[带图表的详细复核报告](reports/replication_v2_results.md)。

| 项目 | 实际结果 |
|---|---|
| 发现供者 |12名病例＋10名对照；150,788个细胞 |
| 外部供者 |15名病例＋7名对照；52,744个细胞核 |
| 表达排序前20 |14同向／20可评价／20原定入选 |
| 网络传播前20 |6同向／19可评价／20原定入选 |
| 外部全基因校正 |0/17,643达到BH q<0.05 |
| 配对供者重采样 |500次中487有效、13失败；方法差值95%条件区间−38.9至+8.7个百分点 |
| 发现逐供者排除 |22/22轮完成；较稳定的网络名单未建立更好外部复现 |
| 临床敏感性 |3个可估计子集完成；2个因人数不足未拟合 |
| 遗传与档案 |40/40基因保留；SKIDA1有限映射；按预定规则产生1个PARP9不确定档案 |

同向比例不是靶点准确率，不显著也不证明没有效应。按供者汇总能避免把同一个人的细胞当作许多独立患者；它不能消除所有混杂。

- [临床敏感性和原始分期说明](reports/clinical_sensitivity_report.md)
- [供者级细胞亚型构成](reports/stromal_composition.md)
- [独立遗传证据](reports/genetic_evidence_report.md)
- 保留的SRV01科学检查点
- [AWS部署说明：未执行](docs/AWS_DEPLOYMENT.md)

Nextflow与独立统计环境已真实运行，输入、参数、名单和输出有冻结及校验记录。本地镜像构建、真实manifest身份和容器流程均已核验：人工测试5个步骤首次运行及5个缓存续跑通过；真实数据5个主要步骤重新执行后，34张科学表和8份结果元数据一致，离散选择完全相同，连续数值最大差约1.10e-13，小于预定绝对误差1e-8。AWS没有创建或运行资源。Docker固定环境、Nextflow管理步骤和缓存；二者不提高生物学结论的有效性。

本项目保留来源分离、作者标签、二级计数汇总和小量未完全归因计数差异的限制，不提供治疗或基因激活/抑制建议。大型输入和计算产物已复制到本地，SRV01原件保留；这些文件仍被Git忽略，只克隆源代码不足以恢复完整检查点。旧版中英文首页的原始字节保存在docs/archive，科学冻结文件保持不变。

核验入口：[原里程碑与交付对应表](docs/LOCAL_DELIVERY_ACCEPTANCE.md)、[真实主要结果复现凭据](receipts/local_primary_container_replay.json)、[最终交付完整性判定](receipts/local_delivery_verification.json)。
