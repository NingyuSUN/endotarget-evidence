# Discovery counts 已准备；M0 尚未全部通过

1. **这一轮做了什么**
从公开 OSF P077 fibroblast 原始分数 Alevin UMI 汇总输入准备可复用矩阵。原始 compact 来源覆盖 23 供者；按已采用临床决定排除 sample73 疾病字段矛盾后，主矩阵为 22（12 ENDO、10 non-ENDO）。五个临床旗标病例/对照本来已不在 compact23 中，不重复计入排除数。严谨组织学阶段敏感性另排除 samples76/123，保留 20（11/9）元数据。与临床代理的最终22供者账本逐字段比对 ID、疾病、年龄、捕获批次、腺肌症和敏感性标记，全部一致。
本轮没有重新下载/解码完整 RDS，也没有差异表达、排序或网络分析。历史 author aggregation 与本轮输入工程核验分开。

2. **怎么做的**
首先冻结自己的 lineage 边界：确切 source fibroblast/decidual stromal labels 纳入，Prv/VSMC mural 标签、其他主要群和缺失标签排除。来源是 Duempelmann 2026 Cell-type identification：Symphony v0.1.1 使用 Tan reference，低于0.5的标签丢弃，随后 menstrual-stage Louvain refinements，原作者用 DotPlot 确认 markers。C7/CXCL14 等标签内 marker 名和 STEAP4/MYH11/CCL19 mural 名来自作者标签，不是本轮重新证明的分类。`config/discovery_lineage.json` 和 `config/discovery_lineage_crosswalk.tsv` 保存方法/marker来源及明确局限。自己的边界逐细胞重算与 P077 成员计数完全相符，primary22有150,788个qualifying fibroblast，最少每供者2,273，22/22均≥50。此为标签成员审核，不是独立 marker expression 或生物身份重新验证。
其次排除全部50不可逆多源合并组（对应177原始features），从34,600 canonical行恢复34,550一对一原始feature IDs，包括9,349带版本ENSG和25,201原始符号/其他ID。没有拆分合并值，也没有把作者canonical别名当作本项目基因ID。
最后固定官方HGNC2026-10-06批准基因注释，只接受唯一官方稳定Ensembl目标。符号先精确匹配approved名称，才尝试唯一official previous symbol；不使用generic aliases、不猜大写名称。原ENSG仅剥离末端数字版本用于官方核对，原ID/版本始终保留。多个源feature对应同一target时整个冲突组剔除，不相加、不择最大值。

3. **结果意味着什么**
保存 source_counts.tsv.gz 为34,550原始features×22供者，counts_stable.tsv.gz 为26,109唯一stableENSG genes×22供者。接受24,676 exact approved symbols、1,069原ENSG和364唯一previous symbols。另排除8,309未映射、129重复target源行（63组）、2个无stabletarget、1个多义符号。实际多义例C18orf21关联HGNC:28802/RMP24和HGNC:57042/RMP24P1，没有任选一个。
每供者因为后续注释排除损失1.7169%–2.2385%的singleton估计UMI；连同50不可逆组，总损失为原P077的1.8432%–2.3796%。这说明丢弃计数规模，不说明被排除基因生物学不重要。原始source counts及完整排除理由保留，后续gene-level DE若采用这一入口，应明确限定为无歧义的HGNC批准named-gene范围；它也不等于之后的STRING网络范围。

4. **验证与问题处理**
6个新增初始边界测试先因功能未实现而失败，完成代码后新增实际固定HGNC C18orf21检查；全部32个项目unittest通过。逐值检查非负/有限，保留原输入十进制字符串；source→singleton→mapped/excluded每个供者的Decimal计数总量严格守恒，写出后逐值逐ID回读一致。两次最终重放输出hash一致。真实服务器准备25.73秒、峰值391,860KiB（约383MiB）；不需要完整41.7GB RDS或大内存反序列化。官方注释16,973,125bytes，固定object generation、官方MD5、SHA256、CC0许可证均保存。程序及测试使用服务器Python3.10标准库，未改变其他环境。
系统先前没有rg，文件定位改用pathlib；没有影响分析内容。没有隐藏失败或把旧历史parser PASS算作本轮测试。

5. **仍未解决的范围**
未独立从全原始RNA sparse slots重算作者donor sums；细胞标签与QC继承作者，不能当作本轮独立生物注释验证。HGNC是批准named-gene注释而非完整Ensembl非编码feature全集；原始counts参考版本未重新确认，剥离gene version不能保证不同reference release中基因模型完全同一。官方previous symbol关系也不是原始mapping坐标证据，局限公开保留。
HECA只做label-blind部分source-ID核对：已缓存22/51来源列中，21,395个discovery stable IDs能直接对应；其余4,714尚无该部分探针的一致对应。剩余29列、完整counts层和共同基因最终冻结尚未完成，不按此探针结果筛选疾病目标。最终设计、外部counts、表达过滤、network/ranking freeze仍由父任务推进；本轮不能宣布M0或科学/临床有效性通过。

6. **接下来的授权工作及理解状态**
本轮没有需要用户决定的事项。父任务可使用标准接口推进final design及外部count gate：稳定矩阵first column gene_id，后续列为exactoriginal sample IDs；donor账本含donor_id、disease、numeric age、capture batch及adenomyosis0/1，另保留20供者敏感性标记。需要双方完整source-ID一致性审核后再冻结共同范围，不能把缺失基因当作阴性。
YELLOW：理解尚未评估。应能解释细胞不是独立患者样本、为何分数UMI此步不取整、为何多源合并不可逆、为什么approved命名范围的丢失与外部coverage限制不等于疾病无效应。

## 文件与复现

- 真实矩阵及provenance：data/processed/discovery/source_counts.tsv.gz、counts_stable.tsv.gz、selected_donor_count_ledger.tsv、provenance.json。
- Feature/排除账本：manifests/discovery_features_mapping.tsv、discovery_features_irreversible_exclusions.tsv、discovery_features_donor_count_ledger.tsv。
- 冻结policy与lineage：config/gene_mapping_reference.json、gene_mapping_policy.json、discovery_lineage.json、discovery_lineage_crosswalk.tsv。
- 原始数据与下载不提交。只提交自己的脚本、测试、上述小型config/manifest/report/receipts；临床代理账本只读取不修改。
- 命令：python3 scripts/fetch_gene_annotation.py；python3 scripts/prepare_discovery_counts.py；python3 -m unittest discover -s tests -v。资源测量完整命令见 receipts/discovery_counts_validation.json。
