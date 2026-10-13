# STRING 蛋白到稳定基因背景图

**本轮已生成159,312条基因边、14,554个稳定ENSG节点；原蛋白图保留，未下载、未计算疾病效应或筛选表达基因。**

1. 做了什么及输入
使用缓存 STRING12.5 protein.info、已固定direct-evidence protein edges和HGNC2026-10-06 approved set。完整input hashes在receipts/string_gene_mapping.json。config/network_gene_mapping.json于2026-10-13T08:47:38.235607+00:00保存后才进行映射，冻结case-sensitive approved-symbol优先、唯一official previous symbol后备、无generic aliases、无case guess。请求模型配置记为gpt-6.1-sol/high。

2. 怎么做及为什么
每个source STRING蛋白ID保留source preferred_name，用精确HGNC批准名称，或唯一official previous symbol，找到唯一Approved HGNC记录及一个unversioned stable ENSG target。没有通过序列、基因组坐标或当前匿名别名推断。多个蛋白允许映射到一个基因；这与RNA feature碰撞时不能再分解其counts的情况不同：此处只转换网络edge endpoints。基因边重复时取原edge最大权重，映射后形成gene selfloops就删除；不把重复source evidence相加。未映射/多义端点的边排除，原因/边权单独保存。

3. 实际结果和含义
19,699 source protein.info IDs中19,301 included proteins，对应19,297 unique stable genes。18,820以approved symbols匹配，481以唯一official previous symbols匹配。382无名称映射、15无合法stable target、1多义。4个基因各对应2个蛋白，保留该many-to-one关系；不是人为选择一个protein。
原161,548 protein edges中排除2,236（1.3841%），最终159,312 gene edges/14,554 nodes。排除原因分别unmapped2178edges、ambiguous46edges、no-stable-target12edges。原操作性edge weights总和139,296.265，endpoint排除1,852.205（1.3297%），最终137,444.060。真实数据新增gene loops和duplicates均为0，不能据此说处理代码没有必要；手算测试覆盖二者。权重之和用于工程守恒，不是独立证据条数、概率或生物效应大小。
Gene图SHA256为`84a4313292da6d6d69b76a4f721893aa89775beebeed9eabcb263fe93e797033`。原protein图仍为`803cb3f406f33a81857c0759878485c41f018a5442b493d14f598b5b71808e09`。

4. 验证与问题
6个手算测试先因未实现而失败，再全部通过：approved/previousname映射、multiple proteins同gene、unknown/ambiguous/caseguess排除、current approved优先、多条protein edges合并max、resultinggeneloop排除、未知endpoint/非法weight失败、duplicateproteinID不覆盖。手算5edges=2endpoint排除+1loop+1duplicatecollapse+1输出；weight4.000=1.600+0.800+0.700+0.900。
真实run严格Decimal edge-count/weight conservation通过，输出geneIDs/weights逐行回读一致，所有原始input hashes运行前后相同。实测9.48秒、peakRSS374185984bytes（约356.9MiB）。提交前Git检查发现successful crosswalk最后一列为空导致trailing-tab报告；已经将该列明确写为not_excluded，保留此前副本并重跑6tests通过，gene图SHA保持相同。按systematic-debugging先定位根因再作仅metadata格式修正；未改mapping decisions或edgeweights。本轮只做targeted tests，未修改其他代理输入或环境。

5. 尚未解决的科学范围
名称映射链是 ENSP → source preferred_name → HGNC，不能证明原protein和基因属于同一genome release、相同序列或转录本/isoform。HGNC snapshot和STRINGprotein.info来源时间不同，仍有gene-model/name变化与剩余mapping错误的可能；排除一条edge不是该gene没有生物关联。Functionalassociation不能当作physicalbinding或causalregulation。
ID-only coverage记录于receipts/string_gene_id_coverage.json：discovery26,109accepted IDs中17476能在name-mapped protein.info对应，13701出现在当前背景图节点中。这里只数ID，不读取expression值、不按expression过滤reference。最终discovery-expressed induced U、其中zero-degree移除、DE和排名仍未做。

6. 下一步和理解状态
本轮没有需要用户决定的事项。父任务可在已冻结表达规则和counts gates通过后，利用该背景图构建discovery-expressed gene induced subgraph，再在该子图中重新判定zero-degree；不得沿用完整背景图degree作为最终筛选。
YELLOW：理解尚未评估。应能解释名字映射不等于序列等价，多蛋白变成同gene为何会形成duplicate/loop，以及edge weight之和只验证工程账本。

交付：data/reference/string/direct_evidence_gene_edges.tsv.gz；manifests/string_protein_gene_crosswalk.tsv；manifests/string_gene_mapping_statistics.json；receipts/string_gene_mapping.json、string_gene_validation.json、string_gene_id_coverage.json。配置、代码、这些小型manifests/receipts/report与测试提交；原data参考和生成图不提交。
复现：python3 scripts/map_string_genes.py；python3 -m unittest discover -s tests -p test_string_gene_mapping.py -v。代码仅缓存读入，已有输出若字节不同就失败并保留old/partial，而非静默覆盖。
