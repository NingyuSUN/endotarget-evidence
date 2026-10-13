# Fixed-panel lineage plausibility audit

**Executed, descriptive only:** the unchanged fixed 23-marker panel was assessed for the 22 discovery donor fibroblast pseudobulks and 22 HECA donor stromal-fibroblast pseudobulks. The same 22 donors' source-defined mural remainder was available after exact checks. The aggregate pattern is broadly consistent with the inherited fibroblast boundary, while residual signals and shared marker expression prevent any claim of anatomical purity or independent cell reannotation. No disease labels, group effects, p-values, target rankings or new exclusion rules were used.

## What was measured

The panel and interpretation were fixed at 2026-10-13T08:40:50Z before marker-value inspection (receipts/lineage_marker_panel_freeze.json). Each donor's CPM is one million times the marker's count divided by the sum of **all retained stable-gene counts in that donor and population (26,109 discovery; 26,125 HECA)**. Every donor then has equal weight. These are relative RNA abundances, not percent cells or purity. Discovery's original fractional estimated UMIs remain unrounded. HECA uses genuine downloaded integer UMIs, without inverse reconstruction; no count input was changed by this audit.

All 23 markers map uniquely to the pinned HGNC stable IDs and are present in all three population matrices. Each summary below covers 22 donors. Quartiles use the inclusive sample-quantile convention; full donor values, mean, median, quartiles, nonzero-count denominators and ranges are in lineage_marker_audit_donors.tsv and lineage_marker_audit_summary.tsv. An absent feature would be reported unavailable, never replaced by zero.

## Discovery: all predeclared markers

| Marker | Primary fibroblast median CPM [Q1,Q3] | Source mural remainder median CPM [Q1,Q3] |
|---|---:|---:|
| DCN | 880.30 [726.27, 1008.69] | 540.82 [453.41, 633.49] |
| LUM | 2254.15 [1853.21, 2796.21] | 1385.88 [1180.01, 1531.60] |
| COL1A1 | 6789.99 [6124.32, 8108.70] | 6699.88 [5827.39, 8105.15] |
| COL1A2 | 2952.31 [2542.01, 3467.97] | 4332.08 [3916.19, 5269.00] |
| PDGFRA | 114.71 [106.05, 140.24] | 15.62 [11.32, 23.95] |
| PRL | 1.70 [1.01, 1.95] | 5.87 [2.68, 9.64] |
| IGFBP1 | 0.22 [0.11, 0.32] | 0.00 [0.00, 0.00] |
| C7 | 2.18 [1.97, 2.80] | 0.19 [0.00, 1.42] |
| CXCL14 | 10.70 [7.05, 16.89] | 95.38 [69.16, 141.51] |
| RGS5 | 21.54 [18.15, 27.83] | 3061.29 [2514.98, 3550.39] |
| PDGFRB | 120.44 [104.25, 139.40] | 411.02 [352.59, 502.63] |
| MYH11 | 2.51 [2.13, 4.85] | 383.69 [299.33, 427.12] |
| ACTA2 | 354.06 [262.95, 443.25] | 7486.37 [6813.46, 8018.33] |
| TAGLN | 524.97 [327.95, 670.93] | 3985.47 [3542.15, 4537.00] |
| EPCAM | 20.31 [16.57, 33.66] | 15.21 [12.06, 23.32] |
| KRT8 | 15.25 [11.61, 20.76] | 35.99 [28.02, 49.18] |
| KRT18 | 23.06 [17.15, 33.29] | 158.96 [124.17, 211.37] |
| PECAM1 | 3.58 [2.11, 5.14] | 5.69 [1.74, 8.48] |
| VWF | 3.45 [2.41, 5.20] | 2.38 [1.56, 10.85] |
| PTPRC | 4.85 [3.56, 7.49] | 4.44 [3.48, 6.74] |
| LST1 | 1.01 [0.87, 1.34] | 0.55 [0.00, 1.58] |
| HOXA13 | 6.43 [4.00, 10.98] | 3.92 [1.65, 6.56] |
| KRT5 | 0.19 [0.06, 0.36] | 0.00 [0.00, 0.00] |

## Interpretation within the fixed limits

DCN, LUM, COL1A1, COL1A2 and PDGFRA have recorded counts in all 22 primary donors. PDGFRA is descriptively more prominent in the primary population, whereas RGS5, MYH11 and contractile markers are more prominent in the source mural remainder. This supports broad mesenchymal/fibroblast versus mural population plausibility. It is not a formal between-population statistical test.

The result also demonstrates why markers cannot be exclusive gates. Collagen genes are abundant in both populations, and CXCL14 is more prominent in the mural remainder despite its presence in an inherited fibroblast subtype name. Primary ACTA2/TAGLN remain substantial; the fixed primary includes dS1-myofibroblast. None of these observations changes eligibility.

Primary HOXA13 has median 6.43 CPM and range 1.43–29.26 across 22/22 donors with nonzero counts; KRT5 is nonzero in 19/22, median 0.19 CPM and maximum 0.94. EPCAM, KRT8/KRT18, endothelial and immune context markers are also nonzero. Their presence does **not** prove or quantify cervical, epithelial or immune cell contamination. It leaves unresolved intrinsic/shared expression, ambient RNA, doublets, minority cells and imperfect inherited labels. Discovery has no independently verified HOXA13-negative per-cell exclusion; do not claim that all cervical-origin cells were excluded merely because P077 matches fibroblast labels. HECA has its separate fixed source label exclusion, not a zero-HOXA13 requirement.

PRL/IGFBP1 are low or variable (one primary donor has zero IGFBP1; its maximum is 40.09 CPM). No donor is rephased or excluded: these are decidual context markers, and proliferative donors need not express them.

## Comparator integrity and scope

Every one of the 34,600 canonical source features across all 23 original compact donors was checked by exact Decimal P005-P077 subtraction (795,800 feature/donor pairs): no negative remainder. The two matrices have identical unique feature/donor identifiers. The original annotation stroma partitions exactly into the frozen fibroblast and mural labels, and its per-donor P005/P077 cell membership matches the source population manifest. The 26,109 accepted stable-gene P077 values exactly equal the prepared primary matrix for every selected donor. Retained P005=P077+remainder mass is conserved exactly per donor; no clipping, rounding or rescaling was used.

This is a remainder of the same secondary OSF aggregation of inherited author labels, not an independent purified cell reference. Exact arithmetic and membership checks establish compatibility, not biological purity. Full counts/identifiers, input/output hashes and validation results are in receipts/lineage_marker_audit.json.

## HECA: same frozen panel, no disease grouping

The 22 HECA donors contribute 52,744 previously selected nuclei. All 23 panel genes are available in the 26,125-gene matrix. The denominator remains each donor's full retained-gene count total, and every donor has equal summary weight. These snRNA abundances are not calibrated to discovery scRNA values: different assays, tissue compartments, populations and retained-gene universes prevent a cross-cohort purity comparison.

| Marker | HECA median CPM [Q1,Q3] | Nonzero donors / 22 | Range CPM |
|---|---:|---:|---:|
| DCN | 212.69 [184.55, 291.12] | 22/22 | 159.53–804.73 |
| LUM | 100.93 [68.92, 135.39] | 22/22 | 49.67–773.45 |
| COL1A1 | 1825.58 [1156.35, 2817.60] | 22/22 | 581.08–4154.51 |
| COL1A2 | 1261.10 [934.70, 1729.40] | 22/22 | 674.80–2482.15 |
| PDGFRA | 254.98 [207.96, 313.26] | 22/22 | 168.03–781.04 |
| PRL | 0.26 [0.16, 0.33] | 19/22 | 0.00–0.65 |
| IGFBP1 | 0.00 [0.00, 0.62] | 9/22 | 0.00–26.51 |
| C7 | 1.99 [1.66, 3.58] | 22/22 | 0.44–23.32 |
| CXCL14 | 1.17 [0.56, 2.09] | 22/22 | 0.11–30.98 |
| RGS5 | 7.21 [5.68, 8.09] | 22/22 | 4.45–12.98 |
| PDGFRB | 112.86 [87.90, 151.76] | 22/22 | 72.74–518.95 |
| MYH11 | 21.17 [18.31, 31.13] | 22/22 | 11.19–51.84 |
| ACTA2 | 54.87 [51.53, 74.46] | 22/22 | 38.94–263.49 |
| TAGLN | 23.75 [18.32, 32.05] | 22/22 | 9.53–757.09 |
| EPCAM | 43.87 [35.74, 52.93] | 22/22 | 26.19–74.52 |
| KRT8 | 14.38 [12.36, 23.84] | 22/22 | 3.56–33.59 |
| KRT18 | 27.73 [19.38, 40.49] | 22/22 | 5.06–60.88 |
| PECAM1 | 3.35 [2.60, 5.14] | 22/22 | 2.16–8.64 |
| VWF | 5.20 [3.64, 7.51] | 22/22 | 2.51–10.47 |
| PTPRC | 12.76 [8.61, 20.29] | 22/22 | 4.00–41.37 |
| LST1 | 1.73 [1.29, 2.23] | 22/22 | 0.58–4.43 |
| HOXA13 | 1.80 [0.86, 4.06] | 20/22 | 0.00–12.52 |
| KRT5 | 0.32 [0.00, 1.02] | 14/22 | 0.00–3.02 |

DCN, LUM, COL1A1, COL1A2 and PDGFRA have recorded counts in all 22 donors, supporting broad stromal/fibroblast population plausibility. Mural/contractile signals also remain: MYH11 median 21.17 CPM, ACTA2 54.87, TAGLN 23.75 with a 9.53–757.09 range. The visible donor variability is preserved, without treating it as disease association or removing a donor.

Nonzero epithelial, endothelial and immune context signals remain (EPCAM median 43.87 CPM; PTPRC 12.76). HOXA13 is nonzero in 20/22 donors, median 1.80 CPM and maximum 12.52; KRT5 is nonzero in 14/22, median 0.32 and maximum 3.02. Excluding the source HOXA13-labelled cluster therefore does not imply zero HOXA13 expression in every retained nucleus. Aggregation cannot separate shared expression, ambient RNA, minority cells, doublets or misannotation. The assay supplies no purity fraction or independent anatomical proof.

PRL is nonzero in 19/22 and IGFBP1 in 9/22; low/absent decidual markers are compatible with the predeclared proliferative context and trigger no rephasing. No new threshold or eligibility change follows from any of these observations. The bounded marker follow-up is complete; source identity, independent per-cell validation and purity limitations remain explicit rather than being converted into a scientific validation pass.

## Gate meaning and remaining limits

The executed audit completes the agreed **limited population-level interpretation** of the protocol's fixed broad-lineage marker check for both cohorts. Independent cell-level lineage validation, anatomical purity and disease-specific cell-mixture effects remain unestablished. Those limitations must be carried into interpretation; no new cell exclusions or outcome-guided thresholds are justified by this audit. The replication check has now been executed with the same unchanged script and frozen panel using --heca-counts. The final HECA count SHA256 is 77121542de07e346ea2112f80bfee6248993c42f90c393ec3154b0a9de529af0. All original discovery/mural table rows are byte-for-byte preserved (receipts/lineage_marker_heca_extension.json).

The source crosswalk is now explicitly qualified in config/discovery_lineage.json as secondary OSF grouping of original author labels. The unchanged discovery preparation script was rerun once to refresh its receipt, and both source-count and stable-count SHA256 hashes were unchanged (receipts/lineage_config_receipt_refresh.json). Source-level participant independence for the final 22+22 ledgers is supported in cohort_independence_audit.md, with no genotype proof claimed.

YELLOW understanding not yet assessed: the user should explain the donor-level CPM denominator, equal weighting, shared marker expression and why neither nonzero contaminant markers nor correct source arithmetic provides a cell-purity percentage. No additional user decision is required in this authorized round; root adjudicates M0 using these explicit limits.
