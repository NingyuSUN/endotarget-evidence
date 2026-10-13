# 固定 HGNC 人类基因注释与共同 ID 的边界

采用官方 HGNC complete approved set 的2026-10-06 quarterly snapshot，共45,233批准记录，GRCh38 reference与alternative loci。下载object generation1791293906485185，文件16,973,125bytes，SHA256 `8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d`，官方Google object MD5 `r0P9VqZK1XwPzHLW72XlFg==`匹配。来源、完整generation URL、许可证本地HTML及其hash在config/gene_mapping_reference.json、receipts/gene_annotation_hgnc_download.json。重复下载脚本只接受固定大小、MD5、SHA与approved schema，已有文件损坏时拒绝复用。

HGNC archive定义、批准符号/previous symbol和Ensembl gene ID字段见[官方归档文档](https://hgnc.genenames.org/download/archive/)；[官方许可证](https://hgnc.genenames.org/about/license/)说明数据为CC0。不是匿名仓库，也没有套用OSF作者org.Hs.eg.db/legacy aliases。当前snapshot晚于原count reference，previous-name/版本模型变化属于公开限制，不因stableID相同断言raw-read quantification equivalence。

明确policy：sourceENSG保留完整版本、lookup时剥离末端版本；symbol严格精确approved先，然后唯一official previous symbol；无Ensembl target或歧义排除；重复target全组剔除。original symbol/features原数据保留，未因不是ENSG而整类丢弃。HGNC named-gene范围并非全Ensembl feature范围。26,109 accepted genes包括24,676批准符号、1,069原ENSG和364previous-name记录，8441排除行见完整mapping。

HECA source raw Ensembl gene IDs优先于当前名称猜测。目前只对缓存22/51 gene_ids来源列作共同ID检查，21,395个accepted discovery IDs能对应，无疾病效应信息；报告reports/gene_mapping_heca_shared_id_probe.json记录输入路径hash和未完成列。不冻结外部共同范围、不改候选列表、不用missing coverage补选其他基因。

计数损失按供者及排除原因公开：reports/gene_mapping_mass_summary.json和manifests/discovery_features_donor_count_ledger.tsv。后续inferential gene family与network comparison universe须另行固定，不把unmapped叫作无效应，也不把程序通过叫作科学验证。
