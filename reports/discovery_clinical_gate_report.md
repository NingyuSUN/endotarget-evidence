# Discovery clinical gate — archived primary supplement joined to GEO

**Result: clinical metadata supports a numerically feasible proliferative comparison, but the final cohort and M0 gate are NOT passed.** A conservative, explicitly proposed clinical-conflict scenario has 11 cases and 9 controls before cell/QC and model-design exclusions. This is newly extracted metadata, not an expression result.

## Sources and what was actually joined

The input is the 17,332-byte archived copy named 41467_2026_73020_MOESM3_ESM.xlsx, publicly distributed inside [OSF DOI 10.17605/OSF.IO/YJM7H](https://doi.org/10.17605/OSF.IO/YJM7H). Its member SHA256 is 5749f05690ffc1c7b6ff7492eba047b5870610c3f6387a47d7342e67538bce69, verified against the archive author's file manifest; the ZIP member CRC was verified by the retrieval agent. Direct byte equivalence with the publisher-hosted original was NOT VERIFIED. The full archive SHA was declared by OSF and was not independently recomputed.

This is an attributed archived primary supplement copy, not a new primary clinical study. Its usable evidence need not be discarded merely because the publisher CDN could not be accessed; its provenance limitation remains visible. The OSF archive's CC-BY notice applies to its authors' own material; original third-party rights remain attached to the reproduced supplement.

Clinical method context comes from the official XML of [Duempelmann et al.](https://doi.org/10.1038/s41467-026-73020-4), Sec3, Sec14 and Sec22, already on SRV01. Original GEO SOFT comes from GSE266265.

- The sheet has 67 rows with nonempty patient identifiers: 60 atlas samples plus 7 PID IHC records. IHC records were excluded from this scRNA cohort; empty formatting rows were ignored.
- Each of the 60 GEO titles matches exactly CTL_sample<number> or ENDO_sample<number>; the sample<number> portion matches one unique Patient sample ID in Sheet1. No fuzzy or positional join was used.
- Disease status and grade agree for all 60/60 GEO-to-workbook joins.
- The selected clinical columns include original lesion/pathology flags, age, hormonal treatment/IUD, day, serum progesterone, pathology phase, author final phase, capture/sequencing batches and minor clinical exclusions. Published analysis-inclusion flags, pseudotime and effect columns were not used as selection criteria.
- All 60 atlas rows have exact age and Hormonal treatment/IUD = FALSE. The paper states atlas donors with exogenous hormones within three months were excluded. Together these provide explicit donor-linked exposure evidence plus a study-level washout rule; we do not claim individual stop dates were measured.
- The OSF secondary PRIMARY_CLINICAL_MANIFEST.csv contains 23 samples. Their GSM, disease, phase, age, batch and hormone values agree with the archived clinical source. Its 23 identifiers exactly equal the original proliferative rows after removing the five flagged minor clinical exclusions. This agreement is a derivation check within one archive, not independent clinical validation. No synthetic QA manifests were used.

Deliverables: [60-row clinical ledger](../manifests/discovery_clinical_candidates.tsv), [counts and IDs](discovery_clinical_gate_summary.json), [source/column/hash receipt](../receipts/discovery_clinical_gate_join.json), [reproducible parser](../scripts/discovery_clinical_gate.py). Each ledger row retains original labels, GSM/title, sheet row, source hash, clinical flags and NOT_VERIFIED final eligibility.

## Pre-effect clinical scenarios, not a frozen analysis cohort

All exclusions below are proposed clinical judgments before reading disease effects. The existing protocol does not automatically incorporate the secondary author's 23-person selection.

| Scenario | Cases | Controls | >=8 per group numerically? |
|---|---:|---:|---|
| Original author proliferative labels | 15 | 13 | Yes |
| Remove five minor clinical flags | 12 | 11 | Yes |
| Also quarantine the contradictory control disease annotation | 12 | 10 | Yes |
| Also quarantine remaining histology/final-phase conflicts | 11 | 9 | Yes |
| Additionally remove adenomyosis cases, sensitivity only | 7 | 9 | **No** |

The first exclusion removes sample81 (control, endometrial polyp), sample108 (case, endometritis), sample115 (control, cervicitis), sample119 (case, endometrial polyp), and sample121 (case, endometrial intraepithelial neoplasm). Retain their original condition text and reasons. This is a consistent proposed removal of flagged co-pathology, not a claim the protocol already prohibited every listed condition. The secretory sample87 is also flagged minor cervicitis but never enters the proliferative candidate denominator.

Concrete unresolved contradictions:

- **sample73 / control:** non-ENDO in both GEO and workbook, yet Endometrioma = 1. Do not relabel as ENDO. Proposed quarantine pending clarification; the 23-person secondary subset includes it.
- **sample76 / case:** final phase proliferative and low progesterone/day9, but pathology Mid Secretory.
- **sample123 / control:** final phase proliferative and low progesterone/day11, but pathology Secretory.
- sample108 and sample121 also have secretory pathology versus proliferative final phase, but are already removed by co-pathology flags.
- sample81 is final proliferative with cycle day16, although paper Sec22 describes day<=15 plus progesterone<2nmol/L. It is already removed by the polyp flag. Do not silently change the phase.
- sample115 has a mixed LateSecEarlyProlif histology string and is already removed by cervicitis. The phase-conflict flag in the machine-readable ledger detects explicit “secretory” wording; it is not an exhaustive pathology-language ontology.

The most conservative 20-person scenario retains cases sample83, sample84, sample91, sample94, sample106, sample111, sample125, sample126, sample130, sample134, sample135; controls sample70, sample71, sample82, sample95, sample96, sample110, sample112, sample114, sample116. These are candidates, not final included participants.

Four of these 11 cases (sample91, sample111, sample126, sample135) carry Adenomyosis = 1. Their presence can confound interpretation of an endometriosis association. Removing all four leaves only seven cases, below the prespecified discovery minimum. This sensitivity must be reported as failing eligibility if used; do not relax the minimum or hide the limitation. The current protocol preserves control conditions and does not automatically exclude all adenomyosis cases.

## Age, hormones, tissue, controls and design

Candidate exact ages are 22–40 years for cases and 25–40 years for controls in the conservative scenario, with no missing age. Age adjustment is therefore feasible in principle, unlike the earlier GEO-only snapshot. No age imputation was performed.

The paper and GEO identify endometrial biopsies obtained before surgery. Individual surgical indication is not in the selected supplement columns; symptoms and disease/lesion fields do not make these controls healthy volunteers. No evidence of repeated patient identifiers appears among the 60 atlas rows, and the primary paper describes 60 women. A cell-level join remains necessary; the ledger does not assume every cell in a secondary object belongs to an eligible donor or lineage.

Conservative candidate capture batches are uneven: only batch14 and batch16 contain both cases and controls; the other seven retained capture batches contain one group. Both sequencing batches contain both groups (first: 5/5; second: 6/4). This is not perfect disease/batch confounding at the metadata level, but adjustment would rely heavily on the two mixed capture batches. Do not equate batch labels with a completed rank/design validation. Freeze the intended batch representation and check full rank and residual degrees of freedom AFTER all cell and missingness exclusions. Sequencing batch is nested within capture batches and must not be added redundantly.

## Required next checks and practical implication

The discovery clinical unknowns have narrowed substantially. It is now possible to defend explicit donor-linked age/exposure and calculate clinical scenarios without reading expression. However, the final author/secondary cell-to-donor and lineage mapping, >=50 eligible stromal cells per donor, raw-count provenance and complete final model design remain unchecked here.

The proposed conflict quarantine still leaves only one control above the minimum, so losing two controls at later cell/QC gates would fail discovery. If the conservative scenario is adopted, it should be documented before effects together with the co-pathology tradeoff and the failed adenomyosis-free sensitivity. Root integration will adjudicate and freeze; this report does not mark any scenario FINAL.

The HECA author repository search also succeeded. At pinned commit/tree 244c282c20c3efd53707e3eb9f58bbd1cbd56ad7 it exposes utils/metadata-donors.csv, utils/metadata-libraries.csv and utils/metadata-libraries_deconvoluted.csv. These three metadata files were downloaded under data/metadata/heca_author and handed to the HECA integration owner. They contain direct age/hormone/phase and donor/library fields; this discovery report does not pre-empt that separate adjudication. Repository tree and download hashes are recorded in receipts/discovery_clinical_gate_heca_*.json.

Validation consists of exact 60-ID joins, disease/grade agreement, selected clinical-column schema/formula checks, secondary23 comparison, preserved original hash and deterministic output replay. These establish extraction integrity only. No expression or clinical-effect validation has occurred.

YELLOW understanding status: not assessed. The important judgment is why “enough candidates” does not guarantee final power/eligibility, why ambiguous disease/phase labels should be resolved before effects, and why removing a co-pathology can make a sensitivity analysis impossible under a fixed donor minimum.
