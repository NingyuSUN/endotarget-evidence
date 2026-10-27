# Public data and readiness audit

**Historical pre-analysis audit.** A later metadata environment, compact count route and author-confirmed donor crosswalk were used for the completed analysis; see the [method card](METHOD_CARD.md) and [input provenance](INPUT_PROVENANCE.md). Container execution was verified later; see the [local execution record](LOCAL_EXECUTION.md).


Checked: 2026-10-13. Labels: **VERIFIED** = source record or directly inspected metadata; **PROPOSED** = project choice; **NOT VERIFIED** = not yet inspected or executed. Public catalog presence alone does not make data analysis-ready.

## 1. Discovery: GSE266265

**VERIFIED, reported by source:** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265) lists an endometrial scRNA-seq study with 60 donors (35 endometriosis, 25 non-endometriosis) and 466,371 cells. The [2026 primary paper](https://www.nature.com/articles/s41467-026-73020-4) reports excluding exogenous hormone use within three months. Its proliferative-phase disease comparison used 12 cases and 11 controls after its filtering. These are original-study counts, not our analyzed sample size.

GEO lists `GSE266265_EndoAtlas.rds` at approximately 38.9 GB. Its smaller epithelial Monocle object does not substitute for stromal counts. Per-sample GEO pages refer to the pooled Seurat object rather than providing a verified small stromal export.

**NOT VERIFIED:** the selected RNA count layer, stromal labels, donor join, age/batch completeness, eligible donor count after our QC, actual import RAM, source-data redistribution terms. Begin with metadata and memory planning, not automatic download of the entire object. The main analysis excludes secretory and ambiguous-phase samples even though they belong to the source atlas.

[Author code](https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas) is a methodological reference. Its README states no formal code license; do not copy implementation. It also contains an accession text/link inconsistency. Use the GEO record and primary publication to resolve source identity and implement independently.

## 2. Replication: HECA newly generated nuclei

Use the newly generated nuclei associated with [E-MTAB-14039](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039), not the entire combined HECA cells atlas. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11387200/) and [official portal](https://www.reproductivecellatlas.org/endometrium_reference.html).

**VERIFIED, direct metadata audit:** the [public SDRF](https://www.ebi.ac.uk/biostudies/files/E-MTAB-14039/E-MTAB-14039.sdrf.txt) has 103 library records. Deduplication by the declared individual field within the nuclei assay identifies 63 donors. The proliferative subset has **16 endometriosis and 6 normal-labeled donors**, before stromal-cell counts, age/comorbidity and count-layer checks. This is a feasible metadata starting point, not a completed independent validation.

The root process independently reproduced this count from 61,604 bytes of public metadata. The [metadata-only receipt](../receipts/metadata_audit_2026-10-13.json) records SHA256 `cdcd25f5519b3ddac7b68a41140ae577b929a856c35732e598fc95fcb0bbf36b`, field names and donor-stratum totals. The downloaded table was held in temporary storage; no patient-level expression matrix was fetched.

Some file-name suffixes do not equal the SDRF individual field. This may reflect multiplexing/naming; it is not yet classified as a source error. Joining by guessed filename donor IDs is prohibited. Compare author-demultiplexed cell metadata with the SDRF; retain conflicts for adjudication.

Recruitment sources differ from the discovery study (Oxford ENDOX/FENOX versus Bern). Confirm the source/donor crosswalk before claiming independence. HECA new cells and new nuclei partially share donors; do not count these as two independent replication studies. Review other benign pathology/control indications from the public supplement.

**VERIFIED, anonymous HEAD access via author portal on audit date:**

| Object | Bytes | Intended use |
|---|---:|---|
| `endometriumAtlasV2_nuclei.h5ad` | 8,739,186,814 | Candidate replication input; includes more donors than primary subset |
| `endometriumAtlasV2_cells_with_counts.h5ad` | 16,012,377,883 | Not selected for v1 |
| `endometriumAtlasV2_nuclei_immune.h5ad` | 656,626,446 | Not a substitute for stromal nuclei |

Base URL: `https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/`.

**NOT VERIFIED:** H5AD count-layer contents, backed/range access performance, stromal cell-level donor mapping and eligible sample counts. The [official BioStudies JSON](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-14039) also describes unnormalized processed count files; assess sample-selective download after resolving library-to-donor mapping. Controlled-access raw reads are not a dependency.

## 3. Public genetics as a separate evidence layer

[GCST90205183](https://www.ebi.ac.uk/gwas/studies/GCST90205183) is the accessible summary-statistics accession cited by the [2023 endometriosis genetics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/). The public version excludes 23andMe. V1 prioritizes the original study's public locus-to-gene/QTL evidence tables and metadata, avoiding an unnecessary full-summary-statistics download when only annotations are used.

**PROPOSED:** add the public mapped evidence to a frozen expression candidate pool. Record ancestry, mapping method and tissue. No nearest-gene-as-causal shortcut, no new MR claim, and no drug efficacy inference. A newer [2026 multi-ancestry study](https://www.nature.com/articles/s41588-026-02582-2) exists; it is not silently combined as independent evidence because constituent cohorts may overlap. Version 1 deliberately pins one genetics source for tractability and provenance. A later evidence update requires a separate overlap audit.

## 4. Other resources and alternatives

- [STRING documentation](https://string-db.org/cgi/help?subpage=api): network resource candidate. Freeze human download version, experimental/database channel semantics, mapping and use terms before ranking. No network file has been fetched.
- [GSE179640](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640) and [GSE213216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216): possible future lesion-context resources, not v1 replication cohorts. Their mixed tissue/treatment contexts complicate the primary same-compartment contrast; both contribute to the pooled HECA cells resource. **Metadata-only audit, 2026-10-27** ([report](../reports/cohort_extension_feasibility.md), [receipt](../receipts/cohort_extension_metadata_audit.json)): natural-cycle eutopic single-cell donors are 3 ENDO / 0 non-ENDO proliferative and 2 / 0 secretory in GSE213216 and 0 / 1 proliferative in GSE179640, so neither study can extend the eutopic donor-level contrasts; GSE179640 offers 8 matched eutopic–lesion pairs under hormonal treatment and GSE213216 4 natural-cycle pairs. No expression data downloaded.
- No matched public drug/genetic perturbation cohort has been selected. V1 uses genetics; perturbation prediction is an extension, not a hidden requirement.

## 5. Readiness register

| Requirement | Current status | Required completion evidence |
|---|---|---|
| Two public source studies and same-tissue contrast | VERIFIED at source level | Freeze exact subset after cell/QC audit |
| Replication phase/disease donor counts | VERIFIED metadata only | Reconcile with cell-to-donor mapping |
| Raw counts in selected objects | NOT VERIFIED | Layer audit and small count readback |
| Independent selected donors | PROPOSED, supported by distinct sources | Written provenance/overlap ledger |
| Primary donor minima after exclusions | NOT VERIFIED | Donor-flow table with reasons |
| Large-object conversion feasible | NOT VERIFIED | Resource measurement and cost estimate |
| Source licenses/terms documented | PARTIAL | Per-resource use/redistribution register |
| Scientific replication and target dossiers | NOT EXECUTED | Fixed-list results, including failures |

No expression data was downloaded for this design. Do not label Stage 0 passed while any required gate remains open.
