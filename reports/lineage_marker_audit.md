# Fixed-panel lineage plausibility audit

**Executed, descriptive only:** the fixed 23-marker panel was assessed for the 22 discovery donor fibroblast pseudobulks. The same 22 donors' source-defined mural remainder was available after exact checks. The aggregate pattern is broadly consistent with the inherited fibroblast boundary, while residual signals and shared marker expression prevent any claim of anatomical purity or independent cell reannotation. No disease labels, group effects, p-values, target rankings or new exclusion rules were used.

## What was measured

The panel and interpretation were fixed at 2026-10-13T08:40:50Z before marker-value inspection (receipts/lineage_marker_panel_freeze.json). Each donor's CPM is one million times the marker's count divided by the sum of **all 26,109 retained stable-gene counts in that donor and population**. Every donor then has equal weight. These are relative RNA abundances, not percent cells or purity. Discovery's original fractional estimated UMIs remain unrounded.

All 23 markers map uniquely to the pinned HGNC stable IDs and are present in both population matrices. Each summary below covers 22 donors. Quartiles use the inclusive sample-quantile convention; full donor values, mean, median, quartiles, nonzero-count denominators and ranges are in lineage_marker_audit_donors.tsv and lineage_marker_audit_summary.tsv. An absent feature would be reported unavailable, never replaced by zero.

## All predeclared markers

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

## Gate meaning and remaining limits

The executed audit completes the agreed **limited population-level interpretation** of the protocol's fixed broad-lineage marker check for discovery. Independent cell-level lineage validation, anatomical purity and disease-specific cell-mixture effects remain unestablished. Those limitations must be carried into interpretation; no new cell exclusions or outcome-guided thresholds are justified by this audit. The replication marker check remains pending until the HECA integer-count stable matrix is ready; the same frozen panel and disease-blind script can process it via --heca-counts.

The source crosswalk is now explicitly qualified in config/discovery_lineage.json as secondary OSF grouping of original author labels. The unchanged discovery preparation script was rerun once to refresh its receipt, and both source-count and stable-count SHA256 hashes were unchanged (receipts/lineage_config_receipt_refresh.json). Source-level participant independence for the final 22+22 ledgers is supported in cohort_independence_audit.md, with no genotype proof claimed.

YELLOW understanding not yet assessed: the user should explain the donor-level CPM denominator, equal weighting, shared marker expression and why neither nonzero contaminant markers nor correct source arithmetic provides a cell-purity percentage. No additional user decision is required in this authorized round; root adjudicates M0 using these explicit limits.
