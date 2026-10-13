OSF ingestion 路线核查（2026-10-13）

VERIFIED：存在真正可取得的计数输入；仅覆盖作者选定的 23 个增殖期供者，不能当作完整 60 供者单细胞 atlas。主 counts 保留分数 Alevin UMI；没有采用作者疾病效应排名、检验阈值或最低组样本数。本轮未运行疾病模型。

1. 输入与方法
公开主入口 https://osf.io/yjm7h/ ，完整包 https://osf.io/download/feukh/ 共 601,984,856 bytes。OSF 声明 archive SHA-256 为 63de832fdb5da21b6f88624423774e4a9b8b1bf8a3706e5cf182e92605832941；未下载完整包，故完整包 SHA-256 未独立核验。采用 256 KiB 尾段读取 ZIP 目录，再按成员偏移提取小文件。每个保存成员已核对 ZIP CRC 和公开清单 SHA-256。所有实际下载和执行仅发生在 SRV01:/home/analyst/te。作者原创材料 CC BY 4.0 已从节点描述与许可证关系定义核实，第三方数据仍受原始权利约束；没有把外部代码接入项目流水线。

2. 真实 counts 与细胞映射
已保存 P001 endothelial、P005 Main::stromal、P077 Derived::fibroblasts 三个全基因 pseudobulk 输入，均为 34,600 canonical genes × 23 original sample IDs。P005/P077 包括所有 23 供者，病例 12、临床非内异症对照 11；没有把其他 37 atlas 供者纳入。用 PRIMARY_CELL_ANNOTATION_MANIFEST.csv.gz 独立逐细胞计数：177,351 unique barcodes，sample 前缀一致，没有重复。Stromal 160,191 cells、每供者最少 2,467；fibroblasts 156,171 cells、每供者最少 2,273。两者的 23/23 供者都超过 50 细胞。Derived fibroblasts 使用保存的明确标签划分，排除 Prv/VSMC mural 标签；此划分仍需项目独立认可，不能把 broad stromal 等同于纯 fibroblast。
每供者 counts 总和与对应细胞 nCount_RNA 元数据总和相符（rtol 1e-12、atol 1e-6）。两矩阵无负值或非有限值，分别含 109,514 / 109,305 个分数值。23 供者主输入包括 14,939 个作者 DEG cell flag FALSE 的细胞，因此不是 flag-TRUE 子集。此检查验证结构、映射和可得到的总量一致性；没有从原始 RNA sparse slots 独立重算全部 counts。

3. 原始来源与稳定 ID 限制
作者保存的原始来源是 GSE266265_EndoAtlas.rds 的 RNA counts，源文件 41,729,517,596 bytes，SHA-256 44a832cdbb7b43ebd3b4a2b6e375c52b08a8059229cc60228ab20a29249cf94f。作者记录 34,727 source features × 466,371 cells、60 samples、1,188,298,946 nonzero estimates，其中 104,078,188 fractional entries。这些全源核验是作者既有记录，本轮未重新下载源 RDS 核实。源 RNA assay 不是 SCT/integrated；作者报告 chrMT features 缺失，不应用其不存在来解释为细胞无线粒体表达。
SOURCE_GENE_MAPPING.csv 保存全部 34,727 原始 feature IDs，但输出已映射/合并为 34,600 行，共 1,794 行改名、50 个多源合并组、净减少 127 行。其中 METAZOA_SRP 合并 53 个原始 features、Y_RNA 合并 28 个。映射表保留来源，但不能把已经相加的 counts 逆向拆回原始 features。需原始 stable feature counts 时，OSF compact exports 不能直接满足，应使用原始 RDS 流式 route 或另一明确原始 counts 来源。

4. 流式代码与本轮测试
已取得并静态审阅 rds_stream.py、inspect_extracted.py、aggregate_author_counts.py、validate_pseudobulk.py、test_aggregation.py 等。解析器按 8 MiB payload chunks 写磁盘，保留第一个 assay 的 counts 与 meta.data，inspect_extracted 检查首 assay 为 RNA；不执行 R 代码。作者既有记录报告 stream extraction peak RSS 2,216,542,208 bytes 和 154.6 s，不是本轮服务器运行结果。
重要移植问题：Linux ru_maxrss 单位为 KiB，而代码按 macOS bytes 处理，因此现有 guard/report 不能直接用于 SRV01；35 GiB ceiling 也超过服务器 31 GiB RAM。还要求 75 GiB disk reserve。使用前必须修复单位、收紧内存上限，并在小型 R serialization oracle/fixture 核对后才执行真实源文件。完整 RDS 下载及解析本轮未做。
本轮用既有 Miniconda（numpy 1.26.4、pandas 2.2.1、scipy 1.15.3）执行作者已审阅 synthetic aggregation fixture：23 synthetic donors、1,035 synthetic cells、60 source genes → 59 canonical genes，fractional estimates 保留，独立 sparse 汇总和 CSV roundtrip 通过，已完成 checkpoint 未重新写入。该 fixture 没有真实病人或疾病模型，也没有测试真实 RDS parser；工程测试不能当作生物或外部验证。Synthetic mitochondrial missing-value warnings 不影响这些计数断言。

5. 临床补充副本与问题
已取得 archived publication MOESM3 XLSX，17,332 bytes、SHA-256 5749f05690ffc1c7b6ff7492eba047b5870610c3f6387a47d7342e67538bce69，CRC 与 OSF file manifest 相符。Publisher 原始字节等同性 NOT VERIFIED；已交临床映射代理。MOESM4 在 tighter steering 到达前已下载，未打开、未用于判断或提交。直接重用临时签名 redirect 地址出现 HTTP400；稳定 OSF URL 的逐请求 Content-Range 检查可成功。Optional 最后文档成员提取仍未完成，不影响已核实公开 license、code 与 counts；不把失败标为成功。

6. 最小下一步与理解状态
本轮没有需要用户决定的常规执行事项。父任务可继续用已取得 MOESM3 副本核查 23 sample IDs/age/phase/batch/疾病条件，并先独立审查 source→canonical 映射。若项目允许审查后的 gene-level aggregate，P077 已是小型真实全基因 counts 候选；若要求未合并原始 feature IDs，则将合法流式代码作有限 Linux 修正、运行 RDS 小型 oracle，再在明确资源预算下下载源 RDS。不得因 author ≥4/group 或 author DE thresholds 替代本项目门槛。
YELLOW：理解尚未评估。需要能独立解释 raw fractional UMI / normalized expression 的区别，按供者汇总与细胞级数据的区别，以及 many-to-one gene mapping 为什么不可逆。未推断 GREEN，标签不代表软件或科学验证状态。

文件与凭证：data/input_routes/osf/independent_donor_coverage.csv；receipts/osf_route_metadata.json、osf_route_ranges.json、osf_route_independent_audit.json、osf_route_fixture_test.json、osf_route_execution.json。原始下载均保留在被 Git 忽略的 data/input_routes/osf/。本轮仅提交自己的报告与 receipts，其他代理改动和原始下载不纳入。

PROPOSED 一对一保守入口（父任务审查后才能采用）：排除 50 个 multi-source canonical groups，可保留 34,550 个可逆 source feature IDs，不需要把合并值拆开。其中 9,349 个是原始带版本 ENSG stable ID；25,201 个原始为 gene symbol/other IDs，不能统称为 Ensembl stable IDs。1,376 个 singleton ENSG 已映射符号；267 个非 ENSG 单行名称变化均仅大写化。全部原始 9,499 ENSG 都带版本，源 ID 没有 pipe。静态检查聚合及最终 reexport 为 untouched per-donor float64 sums → 按 canonical label np.add.at → CSV；没有标准化、取整、表达筛选。实际映射仍需审查命名/版本规则；本轮只是计算候选排除影响，没有改动或生成任何新科学 counts 输入。逐供者影响在 data/input_routes/osf/proposed_singleton_retention_impact.csv，映射审核凭证 receipts/osf_route_gene_mapping_audit.json。
