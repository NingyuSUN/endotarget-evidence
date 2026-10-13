# 临床补充表恢复：未成功取得工作簿

VERIFIED: 本轮只执行官方临床补充文件定位与有界下载，所有运行、临时接收、日志和文件均在 SRV01 `/home/analyst/te`。未下载表达对象，未分析疾病差异表达，未读取 HECA Tables5/6。

## 输入和范围

已保存官方主论文 XML 明确 discovery 的 PMC13385355 / DOI10.1038/s41467-026-73020-4，MOESM3 为 `41467_2026_73020_MOESM3_ESM.xlsx`（Supplementary Data1）；HECA 的 PMC11387200 / DOI10.1038/s41588-024-01873-w，MOESM3 为 `41588_2024_1873_MOESM3_ESM.xlsx`（Supplementary Tables1–9，Table1 为 harmonized sample metadata）。XML来源哈希和原文caption见 summary receipt。这些是论文文件标识，不是已恢复的临床数据。

已核对上一轮 Springer CDN DNS失败及 EuropePMC25MB上限/超时记录。本轮先用单独PMC媒体路径，然后每份官方EuropePMC补充包只做一次提高上限的完整读取，以及一次新策略顺序ZIP读取：目标文件完整时立即停止，而不是下载余下包。每次128MiB上限，整体180秒，连接/读超时40秒。未绕过CAPTCHA，未使用第三方转载，未补造值。

## 实际结果

| 目标 | 方法 | 本次接收字节 | 秒 | 结果 |
|---|---|---:|---:|---|
| 41467_2026_73020 | direct | 0 | 0.231 | HTTPError: HTTP Error 404: Not Found |
| 41467_2026_73020 | full package | 47185920 | 180.0 | TimeoutError: overall deadline 180s exceeded |
| 41588_2024_1873 | direct | 0 | 0.192 | HTTPError: HTTP Error 404: Not Found |
| 41588_2024_1873 | full package | 4718592 | 142.69 | TimeoutError: The read operation timed out |
| 41588_2024_1873 | sequential ZIP local records, stop at target | 0 | 41.059 | TimeoutError: The read operation timed out |
| 41467_2026_73020 | sequential ZIP local records, stop at target | 25624576 | 180.0 | OverallDeadline: overall deadline 180s exceeded |

VERIFIED: 单独PMC媒体两次404；旧域名OA元数据两次404，迁移域名HECA OA元数据也404，详细URL在OA receipts。Discovery流式包遇到的成员顺序为 Fig2.jpg、Fig2.gif、MOESM1.pdf；180秒结束前未抵达MOESM3。这表明本次网络速率下，前置补充PDF阻塞目标恢复；不能据此认定目标文件不存在。HECA最后流式尝试在收到任何字节前读超时。

VERIFIED: `data/clinical/` 仍为空；成功工作簿0/2。本轮总接收字节合计 77529088，这些是未完整的官方包流，未写为XLSX，也没有声称文件哈希有效。下载receipts保留URL、实际接收字节、时间、错误；不会用未完整ZIP代替原工作簿。

## 验证及结论边界

VERIFIED: `python3 -m py_compile scripts/fetch_clinical_supplements.py` 成功。实际下载验证失败：没有可执行XLSX CRC、结构、sheet schema检查的工作簿。因此未列工作表名、元数据表头或供体ID，未解析HECA Table1。脚本中成功路径包含XLSX组件与内部CRC检查，但这条成功路径没有真实工作簿验证，不能称为已通过。

NOT VERIFIED: 样本年龄数值、周期阶段的临床/组织学/转录推断来源、SourceName与individual对应关系、激素停药与重复活检的逐供体状态仍未知。不能据论文研究层面描述认定某个人符合临床入组条件。

## 可续作状态

下载和流式进程均已退出。本轮没有新的同路线重复尝试。下一步需要在相同官方单文件CDN网络可用时取得精确两份工作簿，或由公开官方存档提供直接文件；然后只检查metadata sheets，HECA只读Table1。不能因为本轮恢复失败而放宽临床门槛或以表达结果替代临床状态。

复现命令（会产生新网络请求，不建议未经新路线证据直接重跑）：

```bash
cd /home/analyst/te
python3 scripts/fetch_clinical_supplements.py
python3 scripts/fetch_clinical_supplements.py --stream discovery
python3 scripts/fetch_clinical_supplements.py --stream heca
python3 -m py_compile scripts/fetch_clinical_supplements.py
```

YELLOW — understanding not yet assessed：需要能区分“官方包收到部分字节”、“完整工作簿CRC通过”和“逐供体临床判定有证据”，三者不能替代。本轮没有需要用户决定的事项；恢复失败作为临床入组判定的明确阻塞证据交回主任务。
