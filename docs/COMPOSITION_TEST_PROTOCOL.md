# Stromal-state composition and disease — protocol

Status: rules frozen before any composition statistic is computed (freeze receipt
[`composition_test_freeze.json`](../receipts/composition_test_freeze.json)). Written
2026-11-02, after the v1 results and the five follow-on analyses were known. The published
composition tables had been read descriptively (donor proportions, group means) in the
[composition audit](../reports/stromal_composition.md); no disease model had been fitted on
them.

## 1. Purpose

The v1 disease contrast is a lineage-level pseudobulk: one observation per donor over all
primary stromal fibroblasts. Its null leaves open the objection that disease might change
the **mixture of stromal states** within the lineage while leaving the lineage average
unchanged, or that a within-state change is diluted by mixing. This analysis tests the first
part directly: does the proportion of author-labelled stromal states differ between
endometriosis and non-endometriosis donors in the proliferative phase, with the same
donors, covariates and design as v1?

## 2. Inputs

| Input | Content | Pin |
|---|---|---|
| `reports/stromal_composition_donors.tsv` | per donor and author state: cells, primary cells, proportion (complete grids: 22 × 12 discovery, 22 × 7 HECA) | SHA-256 in [`config/composition_test.json`](../config/composition_test.json); listed in [`summary_figures_v1.json`](../receipts/summary_figures_v1.json) |
| `reports/stromal_composition_denominators.tsv` | per donor primary-lineage cell count | same |
| `reports/final_clinical_design_{discovery,heca}_primary.tsv` | the v1 primary design matrices (intercept, disease, centered age, batch contrasts, adenomyosis) | SHA-256 in the configuration |

No count matrix, label or design is recomputed. Author labels are inherited as published;
nothing is re-annotated.

## 3. Coarse partition, fixed by label name

The two atlases use different vocabularies; the only harmonization is the authors' own
prefix convention, applied before any value is read:

| Partition | Discovery (`AnnotationRefined`) | HECA (`celltype`) |
|---|---|---|
| decidualized-like (d) | `dS1-myofibroblast`, `dS2`, `dS3`, `dS4` | `dStromal_early`, `dStromal_mid`, `dStromal_late` |
| non-decidualized (e) | `eF1`, `eF2`, `eF3`, `eF4-CXCL14` | `eStromal`, `eStromal_MMPs`, `eStromal_cycling` |
| other | `fib C7`, `fib C7 secretory`, `fib C7 secretory_late`, `fib C7-SFRP2` | `Fibroblast_basalis` |

The mapping follows the published names (d = decidualized stroma, e = endometrial
fibroblast or non-decidualized stroma); it is an assumption about the authors' naming, stated
here and not revised after the data are seen.

## 4. Primary endpoint and model

Per cohort: the **decidualized-like share** of each donor's primary cells,
d / n, on the logit scale with a continuity correction,
logit* = log((d + 0.5) / (n − d + 0.5)), regressed by ordinary least squares on the v1
primary design (intercept, `ENDO_minus_non_ENDO`, centered age, batch contrasts,
adenomyosis). Each donor is one observation with equal weight. The estimand is the disease
coefficient on the logit scale, exp(coefficient) being the odds ratio of decidualized-like
versus other primary cells. Two-sided t test on the residual degrees of freedom (10 in
discovery, 12 in HECA), 95% t interval, and the minimum detectable coefficient (80% power,
α 0.05, exact non-central t) from the realized SE and df, reported also as an odds ratio.

## 5. Secondary endpoints

- Every author state: centered log-ratio (CLR) over all states in the cohort with 0.5 cell
  added to each count, regressed on the same design; BH over the states within the cohort
  (12 and 7 tests).
- The non-decidualized share and the "other" share on the logit scale, as for the primary.

## 6. Sensitivities (reported, never change a label)

- Arcsine square-root transform of the share instead of the logit (propeller-style).
- Reduced design: intercept, disease, centered age.
- Donors with fewer than 500 primary cells excluded.
- Adenomyosis donors excluded.
- Leave-one-donor-out: range and sign stability of the primary coefficient.

## 7. Fixed labels

| Label | Rule |
|---|---|
| per cohort `composition_shift_detected` | primary two-sided p < 0.025 (sign reported) |
| per cohort `no_detectable_shift` | otherwise, with the interval and the detectable odds ratio as the bound |
| pair `replicated_shift` | both cohorts p < 0.05 with the same sign |
| pair `no_replicated_shift` | otherwise |
| `NOT_ESTIMABLE` | design rank-deficient or residual df < 5 after exclusions |

## 8. What this does and does not establish

A detected shift would show that disease changes the state mixture within the stromal
lineage at these donor numbers; its absence bounds how large such a shift could be. The
analysis cannot test within-state expression changes, inherits author labels whose
vocabularies and resolutions differ (12 states in cells, 7 in nuclei), has no positive
control of its own, and rests on 22 donors per cohort. Nothing here nominates a gene or a
cell state as a target.

## 9. Order of operations

Freeze (hash protocol, configuration, code, tests; pin inputs; refuse if outputs exist) →
one run of `scripts/composition_test_analysis.py` (Python 3.10, NumPy 2.2.6, SciPy 1.15.2,
local) → published tables, summary, receipt and report. Defects found before any value is
read are handled by a superseding freeze that names the old one; afterwards the results
stand and a correction is published separately.
