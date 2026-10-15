# Separate genetic evidence for the frozen dossier pool

This layer answers whether the already frozen B0/G1 top-20 union has traceable published genetic mapping evidence. It cannot select, rerank, seed, or fit an expression model. It reads neither HECA disease effects nor external concordance. Source-level gene annotation is performed before any dossier judgment.

## Fixed objects and source

`reports/genetic_evidence_frozen_pool.tsv` holds all 40 unique stable gene IDs, frozen list membership and approved symbols. `receipts/genetic_evidence_pool_freeze.json` records the original pre-annotation creation time and source hashes. Every B0/G1 input is checked against the committed discovery freeze. Approved HGNC symbols and official previous symbols are permitted only if uniquely assigned in the existing 2026-10-06 snapshot; generic aliases are excluded. Ambiguity stays unresolved.

The prespecified underlying publication is Rahmioglu et al., Nature Genetics 2023, DOI 10.1038/s41588-023-01323-z, PMID 36914876, PMCID PMC10042257. The primary article verifies GCST90205183 as the **excluding-23andMe** summary-statistic release; the published meta-analysis and its gene tables include 23andMe. These objects must not be presented as identical datasets. No genome-wide summary statistics or restricted participant data are downloaded. PMC and publisher supplement are access routes to one publication, not two independent validations. Later GWAS are outside this fixed annotation layer.

The main meta-analysis includes 60,674 cases/701,926 controls and is reported as 98% European and 2% East Asian ancestry. Functional SMR uses European GWAS subsets; endometrial QTL analyses remove Oxford and Melbourne from GWAS to limit internal overlap. The article's endometrial RNA-seq sample total is inconsistent (368 in Results, 369=206+163 in Methods); retain that discrepancy rather than resolve by guessing. QTL data measure bulk endometrium, uterus or blood, not purified proliferative stromal fibroblasts.

## Evidence hierarchy and boundaries

- Nearest gene and positional overlap indicate location only.
- A reported eQTL links a variant to expression in the specified tissue. It does not itself establish that the same variant drives disease or that changing expression treats disease.
- Published SMR/HEIDI is recorded as such; it is not renamed a Bayesian colocalization probability, perturbation experiment, or causal proof. No new MR, colocalization, fine mapping or gene association is performed.
- Chromatin contact is a separate mapping category. Endometrial cancer cell lines and immortalized cells are not primary endometrial stroma.
- Multi-trait MTAG findings are explicitly separated from the single-trait endometriosis GWAS and are not counted as independent replication.
- A tested, nonsignificant QTL/SMR result is inconclusive under that analysis, not evidence against all effects. No public mapping is `no_available_mapping_in_audited_source`, never a negative biological result. Disease-direction conflicts cannot be manufactured without aligned alleles, effect direction and a defensible mechanism.

No reuse of GSE266265 or E-MTAB-14039 is documented in the retrieved 2023 methods. Individual donor non-overlap across older UK/Australian tissue collections and the expression cohorts is **NOT VERIFIED**; no genotype check was performed. This is a separate published genetic source with residual participant-overlap uncertainty, not certified independent individuals. Multiple analyses reuse the same GWAS/QTL cohorts and count as correlated annotations from one publication.

## Retrieval and validation

Downloads/run paths remain on SRV01 under `/home/analyst/te/data/genetics`. Initial official PMC supplement requests returned HTML and are retained as failed attempts; HTML is never accepted as Excel. Publisher supplementary XLSX HEAD verified 120,578,101 bytes. Available disk before download was 63,536,336,896 bytes; a streaming 130 MiB limit retains >20 GiB host reserve. ZIP members are read selectively with stdlib XML parsing; no environment package mutation is needed. Source hashes, byte counts, URLs, dates, worksheet names and exact cells are retained. Failed direct GWAS Catalog REST metadata request (HTTP 410) does not negate accession verification in the primary paper.

Engineering checks establish complete 40-gene coverage and reproducible transcription of source cells. They do not prove gene causality, therapeutic validity, stromal specificity, or independent biological replication. Variant association, gene mapping, tissue applicability and causal intervention direction are separate claims.
