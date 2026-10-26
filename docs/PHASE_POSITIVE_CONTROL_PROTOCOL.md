# Phase positive-control protocol — v1.0

Date: 2026-10-26. Status: **PRESPECIFIED AND FROZEN BEFORE ANY SECRETORY-PHASE
EXPRESSION DATA WAS OBTAINED.** Machine parameters: [`config/phase_positive_control.json`](../config/phase_positive_control.json).
Freeze receipt: [`receipts/phase_positive_control_freeze.json`](../receipts/phase_positive_control_freeze.json).

## 1. Purpose and boundaries

The completed v1 disease analysis ([method card](METHOD_CARD.md)) found zero
BH-significant genes in each 22-donor cohort and no transfer advantage for network
propagation. A null of that kind has two readings: the disease effect in proliferative
stromal fibroblasts is too small for 22 donors, or the pipeline is insensitive. This
protocol separates them by running the **identical frozen donor-level pipeline** on a
contrast whose biology is known to be large: the proliferative-to-secretory transition
of endometrial stromal fibroblasts (decidualization), within the same two cohorts.

What this is not: it does not revise, rerun or reinterpret v1 outputs; it does not
reopen disease thresholds; it does not nominate targets; a secretory-phase gene is
not an endometriosis gene. Its results are reported whatever they are, including
failure.

## 2. Timing disclosure

Donor eligibility was tabulated from the already-cached clinical and cell metadata of
both studies before this freeze — the same metadata files v1 used. **No secretory-phase
count matrix exists locally at freeze time.** The discovery matrix has 22 proliferative
donor columns and the HECA matrix 22 proliferative donors; their hashes are recorded in
the freeze receipt. Obtaining secretory counts requires new downloads (section 10), and
none has started. This is a retrospective protocol frozen before effects, not a
literature-blind preregistration: the expected direction of decidualization markers is
public knowledge and is the point of a positive control.

## 3. Estimand and population

**Coefficient:** `SEC_minus_PRO` — adjusted secretory minus proliferative normalized log2
expression in donor pseudobulk of the v1 primary stromal-fibroblast lineage. Disease,
age, batch and adenomyosis are nuisance covariates. The sign convention is fixed:
positive means higher in secretory.

**Lineage:** unchanged from v1 — discovery: the 12 refined labels in
[`config/discovery_lineage.json`](../config/discovery_lineage.json); HECA: the seven author
labels `eStromal`, `eStromal_MMPs`, `eStromal_cycling`, `dStromal_early`, `dStromal_mid`,
`dStromal_late`, `Fibroblast_basalis`. Several labels are phase-specific stromal states
(`dS1`–`dS4`, `fib C7 secretory*`, `dStromal_*`). The lineage-level contrast therefore
includes shifts in stromal-state composition as well as within-state expression change.
For a positive control that is the intended biology; it is declared, not corrected, and
donor-level state proportions are reported descriptively.

**Proliferative arm:** exactly the v1 primary donors, no additions or removals.
Discovery: sample70, 71, 76, 82, 83, 84, 91, 94, 95, 96, 106, 110, 111, 112, 114, 116,
123, 125, 126, 130, 134, 135 (12 ENDO / 10 non-ENDO). HECA: EX0045, EX0288, EX0540,
EX0795, FX0011, FX0052, FX0060, FX1097, FX1101, FX1104, FX1119, FX1125, FX1129, FX1136,
FX1146, FX1150, FX1160, FX1236, FX9005, FX9020, FX9024, FX9027 (15 ENDO / 7 non-ENDO).

**Secretory arm, discovery (GSE266265):** author `MenstrualCyclePhase` in
{`secretory_early`, `secretory_mid`, `secretory_late`}; hormonal treatment/IUD FALSE;
minor clinical-exclusion flag 0 (this removes sample87, cervicitis, by the same rule v1
applied to the proliferative arm); no disease-field inconsistency; age present;
≥50 primary-lineage cells after extraction. Author `periovulatory` donors are excluded
as phase-ambiguous and are not used anywhere. Metadata-level candidates, 13 donors:
ENDO sample75, 85, 100, 127, 128, 132, 136; non-ENDO sample88, 104, 113, 117, 131, 133.
Three ENDO candidates (sample75, 128, 132) carry adenomyosis = 1.

**Secretory arm, HECA (E-MTAB-14039):** author `Binary Stage` = `Secretory`; `Stage`
not `Hormones`, `Menstrual` or `Unknown`; hormonal-treatment field missing, read as the
study-defined natural-cycle group exactly as v1 coded it; superficial eutopic
endometrium; nuclei with `discard_lowQC = no`; ≥50 primary nuclei; one sample per donor
by the v1 rule (largest primary-nucleus count, then lexical library then sample ID; no
summing across libraries); adenomyosis coded from the author pathology `A` token as in
v1; donors with conflicting duplicate clinical rows excluded as in v1. Metadata-level
candidates, 20 donors: ENDO EX0277, EX0366, EX0661, EX0789, FX0033, FX1087, FX1100,
FX1130, FX1156, FX1163, FX1170, FX1205, FX1223, FX1233, FX9006, FX9022; non-ENDO
FX0028, FX1106, FX9008, FX9011. FX1170, FX9006 and FX1223 have two secretory samples
each; the largest is taken. The second FX1223 sample has 49 primary nuclei and is not
summed.

**Minimums and design checks:** ≥8 donors per phase group in discovery, ≥5 per phase
group in HECA (the v1 group minima applied to the phase factor); full-rank design;
residual degrees of freedom ≥5. All checks apply to the FINAL matrix after count
extraction and covariate completeness, not only to this metadata list. A failed check
yields `NOT_ESTIMABLE`; no minimum is relaxed and no donor is reclassified to pass it.

## 4. Fixed design hierarchy

Columns are numeric treatment contrasts in the v1 design-file format.

| Cohort | Primary design | Fallback if the primary design is rank-deficient |
|---|---|---|
| Discovery | intercept, `SEC_minus_PRO`, `ENDO_minus_non_ENDO`, age centered on the arithmetic mean of the analyzed donors, `capture_batch` dummies versus batch11, `adenomyosis` | replace capture batch by `sequencing_batch` (second versus first) |
| HECA | intercept, `SEC_minus_PRO`, `ENDO_minus_non_ENDO`, age centered as above, `library_group_Lane` dummies versus the lowest level present, `adenomyosis` | omit `library_group_Lane` |

Which level of the hierarchy applies is decided once, by the metadata-level rank and
residual-df check in the feasibility step, before any count download. The decision is
recorded in the feasibility receipt and may not be changed after effects. If the
fallback is also rank-deficient the cohort is `NOT_ESTIMABLE`. A constant adenomyosis
column is omitted with a receipt, as in v1. No covariate is removed to improve a result.

## 5. Estimator — identical to v1

All settings in [`config/statistical_analysis.json`](../config/statistical_analysis.json)
apply unchanged except the coefficient name and the grouping used by `filterByExpr`,
which is the phase factor instead of disease: edgeR `filterByExpr` (min.count 10,
min.total.count 15, large.n 10, min.prop 0.7), recomputed library sizes, TMM; limma
`voom` with fixed span 0.5; least-squares `lmFit`; robust `eBayes`; two-sided tests; BH
over the complete expression-pass family of each cohort; 95% conditional intervals.
Gene identifiers pass the v1 stable-ID mapping
([`config/gene_mapping_policy.json`](../config/gene_mapping_policy.json)); discovery
fractional Alevin estimates are retained unrounded; HECA integer counts are required.
Cohorts are modeled separately on original counts. No integration or embedding supplies
model input.

## 6. Order of operations and freezes

| Step | Action | Freeze |
|---|---|---|
| F0 | This protocol, config and metadata-level candidate IDs | hashes + git head in `receipts/phase_positive_control_freeze.json` |
| 1 | Feasibility: candidate ledgers, design matrices, rank/df check, download manifest with sizes, resource plan. **Metadata only; no counts read** | feasibility receipt |
| 2 | Discovery secretory counts from the public GEO object `GSE266265_EndoAtlas.rds` (section 10), same per-donor conservation, barcode-order, label-membership and fractional-retention checks as the v1 preparation | counts hash |
| 3 | Discovery `SEC_minus_PRO` fit; B0 = descending \|moderated t\|; B1 and G1 on the frozen v1 network with v1 parameters, seeded by normalized \|t\| | **F1:** full DE table, top-20/top-50 lists for B0/B1/G1, saved coefficient signs, marker-panel results |
| 4 | HECA secretory archives download and preparation with the v1 HECA count pipeline | counts hash |
| 5 | HECA `SEC_minus_PRO` fit; evaluation with the v1 definitions in [REPLICATION_EVALUATION.md](REPLICATION_EVALUATION.md) | results |

No HECA phase effect may be computed before F1 exists. If step 2 proves infeasible
within the resource plan, the protocol degrades to a declared reduced form: HECA-only
P1 and S1; P2 and S2 are reported as `NOT_EVALUABLE`, not as failures and not as
passes. The roles of the cohorts are not swapped.

## 7. Endpoints and pass criteria, fixed in advance

**P1 — within-cohort detection.** Number of genes with BH q < 0.05 for `SEC_minus_PRO`,
reported for each cohort with the family size. Pass: **≥100 in discovery and ≥100 in
HECA**. Counts ≥1000 are additionally labeled "strong" for description only. The exact
numbers are reported whatever they are.

**P2 — cross-cohort transfer, v1 primary endpoint.** The frozen discovery B0 top-20 is
evaluated in HECA exactly as v1 evaluated disease lists: same-direction / evaluable /
selected, denominator 20, no refill, missing-outcome bounds, and coverage-limited if
fewer than 16 of 20 are evaluable. Pass: **≥16 of 20 selected genes have the same sign
in HECA, and ≥10 of 20 have the same sign with HECA q < 0.05.** Top-50 is reported as
secondary.

**S1 — known-marker panel, direction only.** Seven genes with an expected sign, each
attributed to a public source, are read from the full DE tables; genes that fail the
expression filter count as failures and are not replaced.

| Gene | Expected `SEC_minus_PRO` sign | Source |
|---|---|---|
| IGFBP1 | + | canonical decidual-cell marker; Gellersen & Brosens 2014 review |
| PRL | + | canonical decidual-cell marker; Gellersen & Brosens 2014 review |
| DKK1 | + | HECA: `dStromal_mid` upregulates DKK1 |
| LEFTY2 | + | HECA: `dStromal_late` premenstrual marker |
| PLCL1 | + | HECA: `dStromal_early` progesterone-induced gene |
| ESR1 | − | HECA: `dStromal_mid`/`late` downregulate ESR1 |
| PGR | − | HECA: `dStromal_mid`/`late` downregulate PGR |

Pass: **≥6 of 7 with the expected sign in discovery (q < 0.05) and ≥6 of 7 with the
expected sign in HECA.** Five of the seven come from the HECA authors' own
characterization of the labels used to define the HECA lineage, so for HECA this
endpoint is partly dependent on the source annotation; for discovery it is independent.
IGFBP1 and PRL are independent of both cohorts.

**S2 — graph baselines (optional, resources permitting).** B1 and G1 top-20 transfer
and the paired G1 − B0 500-draw external-donor bootstrap, exactly as v1 defined them.
No superiority statement unless the interval excludes zero and coverage rules hold.
S2 is not required for the positive control to be complete.

## 8. Interpretation rules, fixed in advance

- P1 and P2 pass: the v1 pipeline detects and transfers a known stromal program at
  these sample sizes. The v1 disease null is then not attributable to pipeline
  insensitivity; it indicates disease effects below what 22 donors per cohort resolve
  in this lineage. No gene becomes a target because of this result.
- P1 passes, P2 fails: within-cohort detection works; cross-platform transfer
  (scRNA to snRNA, different label vocabularies) is limited, so the v1 transfer endpoint
  had bounded sensitivity. Both facts are reported.
- P1 fails in either cohort: a pipeline, eligibility or count-preparation problem is
  presumed until explained. No further disease analysis is started on this pipeline
  before the cause is documented.
- S1 is explanatory, not decisive. S2 contributes nothing to pass/fail.
- No positive-control number enters any v1 table, figure, dossier or brief; v1 frozen
  artifacts remain byte-identical.

## 9. Prespecified sensitivities, all descriptive

Secretory sub-phase composition (early/mid/late) and the discovery/HECA effect scatter
stratified by sub-phase; donor-level stromal-state proportions by phase; the
scRNA-versus-snRNA measurement difference. Leave-one-donor-out and matched reference
sets are not run unless S2 is executed, in which case the v1 definitions apply.

## 10. Resource plan to be verified in the feasibility step

- Discovery: the OSF compact export covers only the 23 author-selected proliferative
  donors; the only public source of secretory-donor counts is the full GEO object
  (`GSE266265_EndoAtlas.rds`, 41,729,517,596 bytes, SHA-256 recorded by the author
  route). Extraction uses the archived author stream parser after the documented Linux
  fixes (RSS units, memory ceiling, disk reserve), with per-donor conservation against
  `nCount_RNA`, barcode/metadata order and author-label membership checks. Required
  disk and peak memory are measured before download.
- HECA: one public BioStudies archive per selected secretory sample, the same route v1
  used; about 20 archives. Sizes are read from the server before download.
- Nothing is downloaded, and no model is fitted, in the feasibility step.

## 11. Sources

According to PubMed: Marečková et al. 2024, *Nat Genet*
([DOI 10.1038/s41588-024-01873-w](https://doi.org/10.1038/s41588-024-01873-w); HECA labels and
dStromal marker statements); Wang et al. 2020, *Nat Med*
([DOI 10.1038/s41591-020-1040-z](https://doi.org/10.1038/s41591-020-1040-z); widespread
decidualization feature of stromal fibroblasts across the cycle); Gellersen & Brosens
2014, *Endocr Rev* ([DOI 10.1210/er.2014-1045](https://doi.org/10.1210/er.2014-1045);
decidualization markers); Talbi et al. 2006, *Endocrinology*
([DOI 10.1210/en.2005-1076](https://doi.org/10.1210/en.2005-1076); cycle-phase
transcriptome structure). Discovery cohort: Duempelmann et al. 2026
([DOI 10.1038/s41467-026-73020-4](https://doi.org/10.1038/s41467-026-73020-4)).
