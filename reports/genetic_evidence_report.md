# Genetic evidence for the frozen 40-gene pool

**One of 40 genes, SKIDA1 (B0 rank 2), has limited published mapping support in the prespecified 2023 source. The other 39 have no available mapping in the audited gene-mapping result columns. This is not 39 negative genetic tests, and it does not establish that SKIDA1 is causal or therapeutically actionable.** No expression ranking or model was changed.

The pool is the already frozen B0 top-20 plus G1 top-20 union, without overlap. All 40 stable IDs map uniquely to an approved symbol in the locked HGNC snapshot. Genetic annotation began after creating the pool receipt; no HECA disease effects were read. `genetic_evidence_summary.tsv` retains every member.

## Actual evidence and its meaning

| Gene | Actual source coordinate | Finding | Interpretation |
|---|---|---|---|
| SKIDA1 | Supp11.CrediblesetSNPannot!J76, with A76/C76/G76 | Index rs10828249; credible-set variant rs946711; author `Nearest gene` column is SKIDA1 and consequence is `5_prime_UTR_variant` | Positional annotation at an endometriosis locus. A 99% credible set is a set-level fine-mapping object; it is not a 99% probability that this gene or this variant is causal. |
| SKIDA1 | Supp31.EndometriumeQTLs-shared!K11, with A11/H11/I11/J11 | rs11816415; endometrial eQTL beta 0.213136, SE 0.0544234, nominal P=8.99351e-5; A1=T/A2=C as published | Regulatory mapping at a locus shared in the endometriosis/pain MTAG analysis. Bulk UK/Australian endometrium and shared-trait selection limit applicability. No new multiplicity correction, allele harmonization, colocalization, MR or intervention-direction inference was performed. |

The main article Table 2 also lists SKIDA1 under endometrium eQTL at MLLT10/10p12.31. This repeats the same underlying evidence and is **not** counted as another publication. Neither of these two annotations demonstrates that increasing or inhibiting SKIDA1 would help patients. The current germline regulatory direction cannot be equated with a case-control expression difference, which may reflect disease consequence or cellular composition.

No pool member was found in the audited published significant SMR gene columns, closest-gene-to-mQTL columns or chromatin interaction columns. Absence from a selectively reported significant-results table does not show that a gene was adequately tested or exclude an effect. A conflicting disease mechanism was not established for any member; `not_established` is not a finding of no conflict.

## What was actually checked

The publisher XLSX containing Supplementary Tables 1–35 was downloaded on SRV01, 120,578,101 bytes, and hash-checked. Exact headers/cell values were inspected and validated before parsing gene columns. Tables 3, 10–18, 31–32 and 34–35 were audited. Tables 34–35 contain regional UK cis-eQTL test results, not disease-gene colocalization. Their restricted availability, tested-row counts, minimum nominal P and author top-variant rows are retained separately in `genetic_evidence_uk_qtl_context.tsv`; these do not promote a gene to disease support.

The two UK worksheets contain 1,431,756 row records including their header rows. None of the 40 pool IDs occurs in these tables; a separate check against the workbook shared-string inventory also found none of those stable IDs. This does not contradict the SKIDA1 symbol entry in Table 31: that table reports the combined UK/Australian QTL mapping, while Tables 34–35 expose a restricted UK regional subset. It does not show a failed SKIDA1 eQTL replication.

A mechanical whole-sheet string search would incorrectly call ACTB and EGFR hits: ACTB occurs in ARL14EP's gene-function prose at Supp31!L12; EGFR occurs in EPS15's prose at Supp32!M4. The actual gene-ID cells name different genes. Both were explicitly excluded and preserved in `genetic_evidence_rejected_prose_hits.tsv`. This is why matching named gene columns matters.

## Source independence and limits

The [primary 2023 publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/) establishes the source identity and accession. GCST90205183 excludes 23andMe; the published full meta-analysis includes it. The overall GWAS is predominantly European (98%, with 2% East Asian); European subsets underlie fine mapping and functional analyses. Generalization across ancestries is limited. No full GWAS summary file was downloaded.

The published genetic data are separate from the present expression ranking. The article describes removing Oxford/Melbourne GWAS cohorts in its endometrial QTL SMR analysis to reduce its own overlap. No reuse of GSE266265/E-MTAB-14039 is documented in the retrieved methods, but individual donor non-overlap across historical collections is not verified. Portal mirrors and repeated analyses reuse one underlying publication and must not be counted as independent confirmations.

The article reports 368 endometrial RNA-seq samples in Results and 369 (=206+163) in Methods; that source discrepancy remains visible. Bulk tissue, immortalized/cancer cell-line chromatin annotations and blood/uterus QTLs do not establish primary proliferative stromal specificity. No therapeutic claims or drug recommendations are added.

## Problems, checks and remaining work

PMC's supplement link returned an HTML download-preparation page, and the original GWAS Catalog REST endpoint returned HTTP 410. The HTML was rejected, retained, and replaced with the official publisher XLSX route; the accession was verified from the actual primary article. Source URLs, hashes, retrieval times, failed response status and byte counts are in `manifests/genetic_sources.tsv`. The paper HTML timestamp is taken from its saved-file modification time; the supplement download time was explicitly logged.

The full audit took 4 min 55.68 s, with maximum resident memory 288,172 KiB (~281 MiB) and exit status 0 on SRV01. The validation receipt checks exact frozen pool membership, 40/40 coverage, literal source column names, source hashes, two genuine SKIDA1 cells and two rejected prose mentions. These are engineering/provenance checks, not biological or clinical validation. The source scope is one prespecified publication and its published tables, not an exhaustive 2026 literature review.

The next authorized use is to join this separate evidence table to expression dossiers without reranking and to retain missingness and uncertainty. No user decision is needed for this annotation layer. Learning label **YELLOW: understanding not yet assessed**; the key distinction to defend is that association of a variant, assignment to a gene, and the effect of intervening on that gene are different claims.

Reproduce on SRV01 (reuse source files; no duplicate download):

```bash
envs/m0-h5/bin/python scripts/genetic_evidence_fetch.py
envs/m0-h5/bin/python scripts/genetic_evidence_audit.py
envs/m0-h5/bin/python scripts/genetic_evidence_finalize.py
```
