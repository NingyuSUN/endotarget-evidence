# Clinical M0 gate audit — source-only pass

Status: **VERIFIED source extraction; final clinical eligibility NOT VERIFIED; M0 remains OPEN.**
Scope: metadata already stored on SRV01. No expression matrix, disease effect, published disease-gene list or target ranking is an input. Clinical XLSX integration is a separate pending pass; this snapshot must not supersede a later joined clinical audit.

## What the actual metadata establishes

| Item | Observed result | Interpretation |
|---|---|---|
| E-MTAB-14039 SDRF | 103 source records; 68 snRNAseq records; 63 declared nuclei individual IDs | Records and libraries are not independent people. |
| Fixed proliferative nuclei candidates | 24 records, 22 declared IDs: 16 endometriosis, 6 normal | Candidate count only, before age, hormone, identity, lineage/cell and design gates. |
| Filename suffix versus individual | 28/68 nuclei records differ; 12/24 candidate records differ, affecting 12/22 candidate IDs | Naming discrepancy is unresolved evidence, not proof that either label is correct. |
| BioStudies clinical record concordance | All 68 nuclei records have a matching individual/library/assay record and concordant disease/phase/fine-phase values | Same-deposition concordance does not independently validate donor identity. |
| Exact age | 0/22 candidates have exact age in SDRF; all 103 rows say “20 to 45” | Recruitment range is not individual age. Never use midpoint 32.5 or endpoints as a covariate. |
| Hormone exposure | No explicit donor-level exposure/washout field in SDRF | Proliferative is not direct evidence of no exogenous hormones. |
| Repeated candidate ID | FX1101 and FX1236 each appear in two libraries | Count each at most once; relation between records is not established. |
| Discovery GSE266265 | 60 GSM entries: 35 Endometriosis, 25 No Endometriosis | Study-level count agrees with paper; a clinical-table join is still needed. |
| Discovery GEO clinical keys | tissue, endometriosis status, endometriosis grade, Sex | No donor-level phase or age field in these GEO characteristics; cannot count discovery proliferative eligibility from GEO alone. |

Candidate table: [clinical_candidates.tsv](../manifests/clinical_candidates.tsv). One row is one **declared SDRF individual**, not a certified biological donor. Every row retains source names, suffixes, libraries and original SDRF line numbers. Nuclei row ledger: [clinical_gate_sdrf_rows.tsv](clinical_gate_sdrf_rows.tsv). GEO ledger: [clinical_gate_discovery_gsm.tsv](clinical_gate_discovery_gsm.tsv). Counts and source hashes are in [clinical_gate_summary.json](clinical_gate_summary.json) and [clinical_gate_audit.json](../receipts/clinical_gate_audit.json).

## Concrete identity risks to resolve before joining cells

Examples from the unmodified SDRF:

- UA_Endo11282193_EX0540 has declared individual FX1146 and disease normal; UA_Endo11282193_FX1146 has declared individual EX0540 and disease endometriosis.
- UA_Endo11596522_FX1104 has declared individual FX9005 and disease endometriosis; UA_Endo11596522_FX9005 has declared individual FX1104 and disease normal.
- UA_Endo11460962_FX1170 has declared individual FX1141, normal, proliferative; the paired source ending FX1141 has declared individual FX1170, endometriosis, secretory.
- FX1236 appears in UA_Endo11596521 and UA_Endo12865962; FX1101 appears in UA_Endo10281925 and UA_Endo10281926.

Thus a source-name suffix join can change both disease and phase. Preserve all alternatives until the author donor/demultiplex crosswalk determines the mapping. Do not “fix” by swapping labels, majority vote, file ordering, source prefix, or inferred expression.

HECA Methods Sec21 describes two-donor multiplexing, genotype matching with cellSNP/vireoSNP, and removal of doublet/unassigned barcodes. Required join is: **unique library + barcode -> author demultiplex identity -> clinical individual -> biopsy -> eligibility**. A barcode without its library is not necessarily unique.

## Study-level clinical evidence and its boundaries

Primary HECA source: [Mareckova et al., DOI](https://doi.org/10.1038/s41588-024-01873-w); local data/reference/heca_article.xml, Methods Sec10–12, Sec21. BioStudies data/metadata/E-MTAB-14039.json protocol P-MTAB-143093 repeats the recruitment description.

- Natural cycles require no hormonal treatment for at least three months before collection. The resource includes both natural cycles and hormone-treated donors. A particular candidate's natural-cycle membership must be supported by a linked clinical field; a proliferative label alone is insufficient.
- Endometrial biopsies include Oxford ENDOX/FENOX recruitment. The paper also describes Sanger Human Cell Atlasing and Immunology/Subfertility controls. Do not infer a candidate's recruitment study solely from FX/EX/SE prefixes.
- ENDOX/FENOX presence/absence was evaluated at laparoscopy. Additional controls have different ascertainment (clinical/medical history); “normal” is not synonymous with “healthy”. The article notes later histology revealed benign uterine/endometrial pathologies in 15 study participants. Candidate-specific control indication/comorbidity needs the supplement.
- Menstrual phase was assigned by two pathologists, with transcriptomic/cellular-profile fallback when histology was unavailable. Donor-specific ascertainment must be retained; author phase labels are not uniformly independent histology measurements.
- Nuclei in this study are described as snap-frozen superficial endometrial biopsies. SDRF says endometrium; individual biopsy depth/site is not recorded there.
- Five donors overlap between the newly generated scRNA and snRNA datasets. Excluding scRNA from replication avoids counting those assays as extra replication donors, but donor identity and cross-study independence still need explicit provenance.

Primary discovery source: [Duempelmann et al., DOI](https://doi.org/10.1038/s41467-026-73020-4); local data/reference/discovery_article.xml, Sec3 (cohort description), Sec14 (donor criteria), Sec22 (phase evaluation).

- The atlas description reports 60 women (35 cases, 25 controls) excluding exogenous hormones within three months. This directly supports an atlas-level exposure criterion; GEO does not contain individual exposure records. Link each GSM/title to that atlas clinical table rather than pretending GEO contains hormone measurements.
- Minimum recruitment age was 18. Surgical indications were heterogeneous; non-endometriosis donors with adenomyosis were excluded. Individual age, indication and phase are absent from the 60 GSM characteristics.
- Proliferative assignment used serum progesterone below 2 nmol/L and cycle day <=15. Broader stage classification also used histology and transcriptional profiles. Borderline/missing measurements need author labels and ascertainment preserved; do not recompute author phase from expression.
- No expression results were used to infer any eligibility label.

## PROPOSED operational rules — not a frozen analysis cohort

These proposals implement the existing protocol without asserting that eligibility is settled. Resolve against the clinical supplements and author crosswalk, then freeze exact IDs before effects.

1. **Conflicting identity or clinical fields:** keep both originals with source references and mark unresolved. Quarantine the affected donor/library association from primary aggregation until a documented author mapping resolves it. If clinical sources disagree, document the adjudication evidence; do not prefer whichever source yields enough controls. Count a declared ID as a candidate only.
2. **Hormone/phase:** retain the fixed proliferative-only primary contrast. Confirm donor-level no-exogenous-hormone membership and its >=3-month window where the study definition applies. Report whether evidence is a linked study inclusion rule, explicit clinical record, or unknown. Unknown phase/exposure cannot pass primary eligibility. Retain phase method and cycle day/progesterone values; predeclare a histology/clinical-only sensitivity if ascertainment permits, without deriving new labels from expression.
3. **Repeated libraries/biopsies:** sum technical libraries only after the same biopsy and resolved donor are documented. For genuinely repeated biopsies, proposed primary rule is earliest documented qualifying biopsy by collection date; break exact-date ties using a stable biopsy identifier, independent of disease-effect results. If chronology/biopsy identity is unavailable, leave primary inclusion unresolved rather than selecting the library with more cells. Publish excluded repeat records. Each person contributes at most one primary biological unit.
4. **Age and model:** do not convert “20 to 45” to numeric age. If exact donor ages become available, use complete-case age-adjusted analysis as specified by protocol and separately report unadjusted sensitivity; recheck final donor minima/design after missing-age exclusion. If age remains unavailable across a whole study, explicitly record residual age confounding and freeze the resulting design. Partial age availability must not silently trigger dropping age to improve eligibility.
5. **Independence and control provenance:** require recruiting study/site, recruitment period, stable donor/biopsy mapping, and known cross-resource overlap. Different accessions or different countries are useful context but do not prove donor independence. Keep surgical indication/comorbidities and ascertainment; final batch/disease overlap must be checked after all exclusions.

The protocol requires discovery >=8 cases and >=8 controls, replication >=5 cases and >=5 controls, each with >=50 eligible stromal cells/nuclei, plus a full-rank model and residual degrees of freedom >=5. The initial six replication control candidates allow loss of only one before the control minimum fails. No candidate here has passed the cell, design or complete clinical gates.

## Exact annotation fields still needed

Request/find an existing public author annotation file, without initiating contact in this task:

- Clinical crosswalk: accession, original sample title/source, original clinical donor ID, stable person ID, biopsy ID/date or relative order, repeated-biopsy/technical-replicate relation, assay, recruitment study/site/period.
- Clinical eligibility: disease status/stage and ascertainment, tissue site/eutopic status, age at biopsy (or documented unavailable), original phase and ascertainment, cycle day, progesterone with units, natural-cycle/hormone exposure/drug/last-use interval, menopausal/pregnancy/lactation status where recorded, control indication and benign pathology.
- Cell join: original unique cell/nucleus ID, library ID, raw barcode, author demultiplex donor label, genotype-to-clinical donor map, doublet/unassigned status, author QC inclusion, author lineage/cell-type labels and label provenance. Public QC annotation is sufficient; no raw genotypes or private identifiers are needed.
- Batch and independence: processing/capture/sequencing batch, biopsy source study, documented overlaps with other published cohorts. No donor minima may be recomputed as final before this join and the downstream count/lineage checks.

## Verification and limitations

Executed on SRV01: python3 scripts/clinical_gate_audit.py. Structural validation checks required SDRF headers, unique GSM identifiers, per-row BioStudies field concordance and reproducible output hashes. The receipt identifies exact source bytes. These are extraction checks, not biological validation.

The clinical workbook may supply ages, hormone categories, source cohorts or repeat relations; this source-only pass does not claim those facts are absent from all public data. Even a complete workbook cannot alone resolve barcode identity unless its identifiers connect to the author cell metadata. No M0 pass and no inferential scientific result is claimed.

YELLOW understanding status: not assessed. The user should be able to defend donor-versus-library counting, clinical-versus-expression-derived phase evidence, why age ranges cannot support age adjustment, and why a label conflict may reverse disease assignment.
