# Final clinical/sample and design gate — executed before disease effects

Clinical selection and metadata model checks have been executed on SRV01. **Discovery: 22 donors (12 ENDO, 10 non-ENDO); HECA replication: 22 donors (15 ENDO, 7 controls), 52,744 selected primary-lineage nuclei.** Both chosen primary model matrices pass the fixed donor minima, full-rank and residual-df>=5 checks. This is a bounded clinical/sample/design gate, **not completion of M0, scientific validation, or a disease-effect result**.

Inputs began at verified clean Git HEAD 4d05b02. Rules implement docs/PRE_EFFECT_METHOD_DECISIONS.md plus the authorized deterministic one-library-per-person rule. No expression dataset was opened by this script. Existing original H5 obs row order and author metadata are retained.

## Exact selection and lineage

The committed ledgers are manifests/final_heca_clinical_samples.tsv (all 68 H5 samples), final_heca_clinical_donors.tsv (all 63 H5 genotype IDs), final_heca_clinical_selected.tsv (22 selected samples), and final_discovery_clinical_selected.tsv (22 selected discovery donors). They retain original donor/sample/library identifiers, H5 disease/age/phase, clinical pathology and hormone strings, source clinical conflicts, technical fields, reasons and sensitivity membership.

HECA primary labels: eStromal_MMPs, eStromal, eStromal_cycling, dStromal_early, dStromal_mid, dStromal_late, Fibroblast_basalis. HOXA13, mural/perivascular, epithelial, endothelial and immune labels remain excluded. Only author discard_lowQC=no nuclei contribute. These seven labels implement the adopted anatomical boundary; no marker-expression results were inspected.

Author Binary Stage=Proliferative is primary, with original clinical/histological/expression-assisted phase retained. FX1141 has fine Stage=Proliferative Menstrual and only45 primary nuclei: exclude the mixed menstrual category as well as noting its count failure. It is not rescued in the30-nucleus sensitivity. FX1215 has0 primary nuclei and is excluded. Both exclusions precede effects. FX0052's author Proliferative Disordered category is retained with its histology conflict and fibroid/polyp flags visible; it is not silently renamed normal histology.

Actual ages18 (EX0795) and46 (FX0052) are retained. No20–45 age cutoff is imposed from the historical SDRF placeholder.

## Repeated samples: one observed library per donor

| Donor | Selected author sample | Selected primary nuclei | Excluded sample | Excluded primary nuclei |
|---|---|---:|---|---:|
| FX1101 | UA_Endo10281926_FX1101 | 4,657 | UA_Endo10281925_FX1101 | 2,619 |
| FX1236 | UA_Endo11596521_FX1236 | 3,518 | UA_Endo12865962_FX1236 | 122 |

Rule: among clinically eligible samples meeting the current nucleus threshold, take the library with most primary-lineage nuclei; ties use lexicographic library ID, then sample ID. The same rule is rerun for30/100 sensitivities. Biopsy relationship remains NOT_VERIFIED; nothing is called a proven technical replicate or summed across these records.

This discards2,741 primary nuclei and10,578 nuclei of all labels from the two alternative libraries. The number of people does not change. The selected primary total is52,744; minimum per person203. A larger library can be more representative technically but selection by recovery can still favor a particular biological/technical composition; this is a pre-effect deterministic choice, not proof against selection bias.

Exact source indices: work/final_clinical/heca_primary_h5_obs_indices.npy, int64, zero-based ascending H5 obs row positions,52,744 entries, SHA256 d02ed9b8cee2c789c4f7aec7a3792a8f9444307163531933b632845cecf864ef. The file is ignored generated data and is not committed. The selected-sample ledger SHA256 is1a0ba5813046f9cb6e18862e17ca7fe64c31a39605d048309f1e3522447046cc. The count-ingestion owner received both paths, hashes and label definitions. Recreate indices with envs/m0-h5/bin/python scripts/final_clinical_gate.py.

## Clinical coding and what it does not mean

All26 author-proliferative H5 sample records have no explicit drug conflict; H5 hormonal treatment remains literal nan and clinical CSV remains NA. The eligibility basis is **study_defined_natural_cycle_group**, as adopted before effects from the paper's>=3-month hormone-free definition and separate hormone-treated group. Missing drug strings are never recoded FALSE, and per-person stop dates remain unavailable.

Pathology retains the source header Endometrial pathology_A-F-H-P and its original token strings. The operational crosswalk interprets A as adenomyosis, F as fibroids, H as hyperplasia, and P as polyps, supported by the author header and the primary paper Sec11's list of those benign pathologies. No separate prose code-key or new pathology revalidation was obtained. Only exact token A is used for the adenomyosis covariate; it is present in EX0540, FX0060 and FX1236, all cases. C/E tokens are retained as given, not substituted for independent diagnostic evidence.

Among seven selected controls, FX0052 has F+P and FX1104/FX9027 have P; the other four have C. They are source-labeled non-endometriosis controls with recorded conditions, not healthy volunteers. Selected cases also include F/H/P co-pathology. These conditions are preserved; this branch does not add a post hoc blanket exclusion. Discovery minor-clinical flags follow its explicitly adopted exclusion policy, so pathology handling is source-specific and leaves residual comparability limits.

One nonprimary donor, FX1158, has two author CSV rows with differing Stage and pathology (F versus E+F). Both original records are saved in the sample ledger's clinical_duplicate_rows_json; the donor is menstrual and ineligible regardless. The script initially stopped on this duplicate instead of silently overwriting it, then explicitly quarantined the conflict. An unrelated historical Fonseca duplicate in the all-study deconvolution CSV was outside the68 target H5 samples; the join now checks uniqueness only for the actual target sample IDs. These changes do not alter the selected22 HECA donors.

## Technical covariates and executed model checks

The official author library table supplies **Lane**,10x kit,Treatment,Processing and library identifiers. It does not separately document a capture-batch/sequencing-batch field or a physical definition of the numeric Lane values. Preserve Lane literally and use it as an author-recorded categorical library grouping, not a numeric dose, date or invented sequencing lane. All selected Processing values are Frozen.

Seven Lane groups remain. Five contain both disease groups:42196 (1/1),43103 (1/1),43281 (4/2),43937 (3/1),43996 (5/1), stated as cases/controls. Lane20 has one case;42058 has one control.10x kit and preparation Treatment are both fully nested in Lane here, so adding their dummy columns would be redundant. The primary adjustment uses Lane to retain that source-defined technical grouping.

An individual library is not automatically an appropriate batch factor:16 distinct library IDs for22 donors produce a model of rank19 with only3 residual degrees of freedom, below the required5. That alternative was executed and rejected for insufficient residual df; disease identifiability alone would not make it acceptable.

| Executed subset | Cases / controls | Model terms beyond intercept | Rank / columns | Residual df | Structural gate |
|---|---:|---|---:|---:|---|
| Discovery primary |12 /10|disease, age, capture batch, adenomyosis|12 /12|10|Pass|
| HECA primary |15 /7|disease, age, author Lane group, adenomyosis|10 /10|12|Pass|
| Discovery histology-conflict exclusion |11 /9|same terms|12 /12|8|Pass|
| HECA histology-conflict exclusion |10 /3|same terms, observed factor levels|8 /8|5|**Fail: controls<5**|
| Discovery adenomyosis-free |8 /10|disease, age, capture batch|11 /11|7|Pass|
| HECA adenomyosis-free |12 /7|disease, age, author Lane group|9 /9|10|Pass|
| Discovery strict phase + adenomyosis-free |7 /9|disease, age, capture batch|11 /11|5|**Fail: cases<8**|

“Histology-conflict exclusion” removes explicit secretory disagreements, including the mixed pathology wording in FX0052; it does not require a positive histology result for every remaining donor. Missing histology remains missing. It cannot be described as a uniformly histology-confirmed subset.

The30/100 primary-nucleus sensitivity selections were recomputed from all clinically eligible HECA samples using the same one-library rule; both retain15/7 and the same model rank/df because the primary minimum is203 and mixed-phase FX1141 remains clinically excluded. Discovery cells all exceed100 under the supplied P077 labels; the separate ingestion owner is freezing its source-label crosswalk.

All models use centered age and an ENDO-minus-non-ENDO coefficient. Adenomyosis is included when it varies, omitted when constant after the corresponding cohort exclusion. The primary cases include3 HECA and4 discovery adenomyosis-positive donors. The disease coefficient is an adjusted association under an additive model, not isolation of a causal endometriosis effect. Full rank and df demonstrate mathematical estimability, not power or removal of all clinical/technical confounding.

Actual ordered design matrices are reports/final_clinical_design_*.tsv. reports/final_clinical_gate_summary.json contains donor order, columns, rank, degrees of freedom, overlap, exclusions and sensitivity failures. Failed sensitivities must not be fitted as confirmatory analyses or rescued by changing minima.

## Verification, outputs and remaining work

Executed: envs/m0-h5/bin/python scripts/final_clinical_gate.py. Independent metadata validation rechecked52,744 source indices against sample IDs/labels, one selected sample per person, deterministic replay, both primary ranks, the fixed-minimum failure, repeated-library choices, preserved missing-hormone labels and unchanged source bytes. A deliberately disease-confounded technical grouping is rejected by the identifiability guard; individual-library df failure is also retained.

receipts/final_clinical_gate.json records source/output hashes, the source-decision document hash, script hash and runtime. receipts/final_clinical_validation.json records validation outcomes. These are extraction and model-algebra checks.

The source-specific clinical/sample/design packets are ready for root integration and final freeze. External raw counts, common gene identifiers, complete annotation crosswalk integration, independence assessment and network/ranking freezes remain separate M0 requirements. No disease effects, rankings or full M0 pass are claimed.

YELLOW understanding status: not assessed. Defend one person versus two libraries, source-defined hormone eligibility versus measured drug absence, author technical grouping versus invented batch, and a full-rank model versus biological validity. No new user decision or approval is required for this authorized bounded implementation.
