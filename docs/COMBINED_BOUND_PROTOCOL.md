# Two-cohort combined bound on the proliferative-phase stromal disease contrast — protocol

Status: rules frozen before any combined statistic is computed (freeze receipt
[`combined_bound_freeze.json`](../receipts/combined_bound_freeze.json)). Written 2026-11-02,
after the v1 results, the phase positive control, the minimum detectable effect (MDE)
analysis and the secretory program contrast were known. No combined estimate, p value or
count had been computed when this protocol was frozen.

## 1. Purpose

The v1 analysis fitted the disease contrast (ENDO minus non-ENDO, author-adjudicated
proliferative-phase eutopic stromal fibroblasts, donor pseudobulk) separately in the
discovery cohort (GSE266265, 12 / 10 donors) and the external cohort (HECA E-MTAB-14039,
15 / 7 donors) and found no gene at BH q < 0.05 in either. This analysis asks one
prespecified question of the same two fits: **what does the evidence of both cohorts
together exclude, gene by gene, and does any gene reach the family-wide threshold when
their precision is pooled?** It is a bound-tightening analysis with a fixed discovery rule,
not a search for targets, and it adds no donor.

## 2. Inputs

Only the two published per-gene tables of the MDE analysis, which carry the saved v1 fit
quantities, pinned by SHA-256 in [`config/combined_bound.json`](../config/combined_bound.json)
and in the [MDE results receipt](../receipts/mde_analysis_results.json):

| Input | Columns used | Origin |
|---|---|---|
| `reports/mde_analysis/v1_discovery_mde_per_gene.tsv.gz` | `gene_id`, `logFC`, `moderated_SE`, `df_total`, `fit_estimable` | v1 discovery disease fit (limma voom, robust eBayes), 17,733 genes |
| `reports/mde_analysis/v1_heca_mde_per_gene.tsv.gz` | same | v1 HECA disease fit, 17,643 genes |

`logFC` is the coefficient `ENDO_minus_non_ENDO` on the log2 scale, `moderated_SE` the
moderated standard error and `df_total` the residual plus prior degrees of freedom of the
moderated t test. No count matrix, design or model is touched; nothing is refitted.

## 3. Family

The **combined family** is the set of genes present and `fit_estimable` in both tables.
Genes present in one table only are listed with their single-cohort values and are not
combined. Multiplicity is controlled by Benjamini–Hochberg over the combined family.

## 4. Primary estimator

For each gene, with cohort estimates β₁, β₂, standard errors s₁, s₂ and degrees of freedom
ν₁, ν₂:

- weights wᵢ = 1 / sᵢ²; combined estimate β̂ = Σ wᵢ βᵢ / Σ wᵢ; combined standard error
  ŝ = (Σ wᵢ)^(−1/2) (fixed-effect inverse-variance combination);
- degrees of freedom by Welch–Satterthwaite for the weighted combination,
  ν̂ = (Σ aᵢ² sᵢ²)² / Σ (aᵢ² sᵢ²)² / νᵢ with aᵢ = wᵢ / Σ wᵢ (equal precision gives ν₁ + ν₂);
- two-sided test of β = 0 on t = β̂ / ŝ with ν̂ degrees of freedom; 95% interval
  β̂ ± t₀.₉₇₅,ν̂ ŝ; family-wide (Bonferroni) 95% interval with α / family size;
- heterogeneity: Cochran Q = Σ wᵢ (βᵢ − β̂)² with one degree of freedom and its p value;
  I² = max(0, (Q − 1) / Q), descriptive. With two studies no between-study variance is
  estimable; **no random-effects model is fitted**, and the fixed-effect estimate is read
  as a precision-weighted average of two platform-specific effects, not as a common effect.

## 5. Bound quantities

- **Minimum detectable effect** of the combined test per gene: the smallest |β| detected with
  80% power at the realized ŝ and ν̂, solved exactly on the non-central t distribution, for
  α = 0.05 (one prespecified gene) and for the BH boundaries at exactly k = 1, 10, 100 and
  1,000 rejections in the combined family (α = 0.05 k / m). The closed-form approximation
  (t_crit + t₀.₈) ŝ is reported alongside.
- **Interval exclusion**: the number and share of genes whose per-gene and family-wide
  intervals lie within ±0.585 (1.5-fold), ±1 (2-fold) and ±2 (4-fold) log2.
- Comparison with the single-cohort values from the same tables (median MDE and exclusion
  shares), to state how much the pooling tightens the bound.

## 6. Descriptive agreement

Share of family genes with the same sign of logFC in both cohorts, and the Pearson and
Spearman correlations of logFC across cohorts over the family. Under the null each is
about 50% and 0; they are reported as description, not tested.

## 7. Sensitivities (reported, never change a label)

- Stouffer weighted z: zᵢ from each cohort's two-sided p value with the sign of βᵢ, weights
  1 / sᵢ, combined z = Σ wᵢ zᵢ / (Σ wᵢ²)^(1/2), normal p value, BH over the family.
- Normal-approximation inverse-variance test (same β̂ and ŝ, standard normal reference),
  to show the size of the degrees-of-freedom correction.
- Primary estimator restricted to genes with Q p ≥ 0.05 (descriptive: how many genes are
  removed; the label is not recomputed).

## 8. Fixed labels

| Label | Rule |
|---|---|
| `no_combined_discovery` | 0 genes at BH q < 0.05 in the primary test |
| `combined_discovery` | ≥ 1 gene at BH q < 0.05 in the primary test; all such genes are listed with both cohort estimates, Q and the interval. A listed gene is a hypothesis generated by pooling the same two cohorts, not a replicated finding and not a target |

The 20 genes with the largest |t| are published descriptively in every case.

## 9. Order of operations and re-issue rule

1. Freeze: hash the protocol, configuration, scripts and tests, verify the input pins
   against the configuration and the MDE results receipt, record the git head; refuse if any
   combined output exists.
2. Run `scripts/combined_bound_analysis.py` once on the pinned tables (Python 3.10, NumPy
   2.2.6, SciPy 1.15.2, the versions of the v1 runtime; executed locally, versions recorded).
3. Publish per-gene table, cohort-only list, summary, top-20 table, figure, results receipt
   with every output hash, and the results report.

If a defect in the code is found before any value is read, the freeze is superseded by a new
receipt that names the old one; if found after, the results stand as reported and a
correction is published separately.

## 10. What this establishes and does not

It establishes, under rules fixed in advance, the tightest bound the two v1 cohorts give
together on proliferative-phase stromal disease differences, and whether any gene crosses
the family-wide threshold once both cohorts' precision is used. It does not add data, does
not estimate a common effect across a single-cell and a single-nucleus platform beyond a
precision-weighted average, does not separate composition from within-state change, and
does not nominate a target. The MDE values remain conditional on the realized moderated
variances of the original fits.
