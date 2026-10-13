# Fixed lineage marker audit — pre-effect decision

This decision adds a bounded **population-level plausibility audit** to protocol section3's broad-lineage marker check. The marker panel is fixed before inspecting marker count values. It does not claim independent per-cell validation, train a classifier, or modify any cell/donor eligibility.

## Inputs, methods and rationale

Read the existing22-donor discovery stable-gene count matrix. Use the already pinned official HGNC snapshot to link exact panel symbols to stable IDs. The primary population is the fixed P077-derived fibroblast set. Original Alevin fractional counts remain unchanged. The same panel can be applied to the separately prepared HECA donor count matrix when available; this is not a reason to download more data.

The23-marker panel is in config/lineage_marker_panel.json:

- Fibroblast/ECM context: DCN,LUM,COL1A1,COL1A2,PDGFRA. Decidual context: PRL,IGFBP1. Retained-label context: C7,CXCL14.
- Mural/contractile context: RGS5,PDGFRB,MYH11,ACTA2,TAGLN.
- Epithelial: EPCAM,KRT8,KRT18. Endothelial: PECAM1,VWF. Immune: PTPRC,LST1.
- Cervical context: HOXA13,KRT5.

These are fixed from standard broad-lineage biological rationale, the original labels and already cached author methods, not from disease-associated gene lists or observed audit values. Evidence levels are explicit per marker: the cached discovery Methods Sec20 describes Symphony/Tan transfer, confidence filtering and original author DotPlot confirmation, but does not provide this complete marker list. The HECA Sec3 anatomical paragraph specifically supports HOXA13/KRT5 cervical concern; HECA epithelial spatial descriptions support EPCAM. Do not claim the remaining general markers were explicitly validated by those paragraphs.

For each donor and population compute gene CPM =1,000,000*gene count/sum of all retained stable-gene counts in that donor/population. The denominator is the complete mapped matrix, not the marker panel or cell number. Report every donor, plus equal-donor-weight mean, median, quartiles, min/max and number with nonzero recorded counts. More cells or a larger library do not give a donor more summary weight. Zero is an observed zero only when the feature is present; absent IDs are unavailable.

## Optional cached comparator

The original stromal annotation rows partition into the fixed fibroblast labels and source mural labels; no new biological labels are inferred. P005-P077 can serve as a **source-defined mural remainder** only after exact unique feature/donor correspondence, nonnegative Decimal subtraction at every feature/donor, and per-donor total conservation. Use the same accepted stable-gene mapping and retained-gene denominator as the primary population. Do not compare panels with inconsistent feature universes. If any condition fails, report why and omit the comparator; no negative-value clipping, pseudocount correction or count rounding is allowed.

The comparator is derived from the same atlas and inherited annotations, not an independently purified reference. Per-donor matrix comparison is descriptive only; no disease labels, p-values, donor-class effect, target ranking or cell-level threshold enters the audit.

## Interpretation fixed in advance

A pattern broadly consistent with fibroblast labels is supporting population-level context, not purity proof. Mural/epithelial/immune/cervical signals may arise from intrinsic/shared expression, ambient RNA, doublets, rare cells or annotation error; pseudobulk cannot distinguish them. MYH11/ACTA2/TAGLN are not automatic exclusion triggers, particularly with retained dS1-myofibroblast. Low PRL/IGFBP1 in proliferative donors is not an identity failure. HOXA13 or KRT5 alone cannot establish anatomical contamination. Absence of signal does not prove absence of a population.

No numerical purity pass/fail cutoff, new eligibility filter or marker-driven reclassification will be derived from these values. Unexpected signals will be reported for interpretation and, if severe, explicit protocol adjudication; they will not silently alter the frozen candidate population. Cross-cohort scRNA/snRNA differences further prevent direct purity calibration.

This bounded check satisfies the **limited broad-population plausibility** interpretation of the marker-audit requirement. It does not satisfy a claim of independent single-cell lineage validation. That scope clarification must remain visible in M0 and the final report. Full-RDS retrieval or cell reannotation is not required for this audit, and no outcome-guided expansion is permitted.

YELLOW understanding status: not assessed. The user should be able to explain equal-donor weighting, the CPM denominator, why a marker is not an exclusive identity label, and why inherited annotation plus an aggregate check cannot establish anatomical purity.
