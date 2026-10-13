# STRING 固定网络过滤审核

**本轮规则测试和缓存重放通过；网络边表未改变，未使用任何疾病输入。**

1. 这一轮做了什么
输入只读既有 STRING12.5 human protein.info 和 protein.links.full 缓存、config/network_reference.json、receipts/string_source_probe.json。原下载、ETag、length、SHA256、gzip CRC 以及首次77.799秒运行凭证完整保留在receipts/string_reference_download.json。本轮没有重新下载sourcefiles，也没有读入expression/disease/DE/ranking表。

2. 方法和修改理由
把原始main内过滤代码提取成纯streaming filter_edges函数，输入为一行一行的网络、允许protein集合和冻结policy。原末尾常量assert不能证明真实路径，已经移除。当前实现严格检查species9606/version12.5、direct experiments/database、1000分制、700阈值、max直接证据/1000、undirected最大重复权重和排除selfloops。配置若要求其他渠道、scale/threshold、weight或duplicate/selfloop语义，会失败而不是忽略。Transferred、textmining、coexpression和combined_score均不能提高筛选分数。达到直接证据阈值的protein必须同时在人类protein.info白名单及9606前缀内，否则失败。
主入口仅接受已有固定hash缓存；missingcache立即失败，不触发networkrequest。生成一个新临时replay副本，只有和指定SHA完全相同才通过；原边表已相同则保留原边表，删除本轮自己的相同临时副本。原下载/runtime JSON不覆盖；重放写receipts/string_reference_cached_replay.json。

3. 实际结果及含义
10,826,054 source directed rows中，10,502,958低于directscore700；323,096 qualifying directed rows合并为161,548 undirected edges、14,658 protein nodes。原source当前未出现qualifying selfloops，但手算测试明确覆盖selfloop删除。边表SHA256仍为`803cb3f406f33a81857c0759878485c41f018a5442b493d14f598b5b71808e09`，所有summary计数与首次运行一致。权重是max(experiments,database)/1000这一操作定义，不是经过校准的interaction概率。

4. 实际测试与出现的问题
新增7个手算边界测试先全部因缺少纯过滤函数而失败，重构后全部通过：高combined/transferred/text-only不能救回低direct；699失败而700在任一directchannel通过；AB/BA及同向重复取max；selfloops删除；unknown/nonhuman qualifyingprotein失败；非法配置策略失败；无效score/malformedrow失败。测试使用人为边，不能当作科学证据。
完整套件先在系统Python遇到HECA NumPy缺失，切换既有m0-h5环境，没有安装依赖。该快照47tests中8HECA测试错误，由HECA代理确认是其并行TDD故意失败基线，implementation仍在进行；本轮7个STRING及其他已实现测试没有失败。这些8错误的名称/日志和环境限制保存在新cached-replay receipt，未隐藏，也未修改其他代理的实现。父任务随后要求只做targeted STRING验证，后续整套重跑由父任务负责。

5. 仍未完成什么
该产物是固定protein功能关联参考边，不是最终discovery-expressed gene universe、physicalbinding或causalregulation。Gene-protein映射、诱导子图、零度节点排除、疾病模型与排名仍未执行。本轮仅证明过滤代码执行固定规则及字节结果不变，不能宣布科学或临床验证成功。STRING来源原许可记录为CC BY4.0及provider条款，见[官方许可页](https://string-db.org/cgi/access?footer_active_subpage=licensing)；本轮复用既有许可来源，不另作再授权。

6. 下一步及理解状态
本轮没有需要用户决定的事项。父任务可继续把固定蛋白边映射到无歧义gene集合，再按研究协议构建discovery表达gene induced subgraph；不能因为原参考图中有边就认为其在最终gene universe里非孤立。
YELLOW：理解尚未评估。需要能解释direct与transferred/combined分数的区别、为什么阈值700边界必须固定、AB/BA合并为何使用最大允许证据，以及functionalassociation不证明因果调控。

复现命令：python3 -m unittest discover -s tests -p test_prepare_string_reference.py -v；python3 scripts/prepare_string_reference.py。后一命令仅重放缓存，不覆盖原下载/runtime receipt。实际缓存重放elapsed/RSS见receipts/string_reference_cached_replay.json。全部任务执行及下载留在SRV01:/home/analyst/te。
