# Published bulk disease signatures in the single-cell stromal results — results

**Labels under the frozen rules: T2 `concordant` for the Poli-Neto 2020 stage I–II and the
Bunis 2022 proliferative signatures, `no_detectable_concordance` for the stage III–IV and the
consensus sets; T1 `no_enrichment` for all four sets; T3 `no_detectable_difference` in both
cohorts for all four sets, pair `not_replicated`.** The published whole-tissue directions of
4,661 and 2,383 family genes agree with the pooled stromal disease estimates in 58% and 55% of
genes against 48% and 47% for expression-matched genes carrying the same labels (matched
percentile 100, binomial p < 10⁻⁴), so a small part of the tissue-level signal points the same
way within stroma. Its size is small: the donor-level stromal signature score differs between
endometriosis and control donors by −0.01 [−0.08, +0.07] (discovery) and +0.08 [−0.06, +0.23]
(HECA) log2 units per gene for the largest signature, where 1 would be the full published
twofold effect; the intervals exclude a stromal share above about a fifth of the published
per-gene effect in either direction. The signature genes are not more strongly disease-associated
than other genes (T1), and the two cohorts' score differences point in opposite directions.

Completed 2026-11-02 under the [protocol](../docs/PUBLISHED_SIGNATURE_PROTOCOL.md) frozen the
same day after the sets were built and before any single-cell statistic was read for them
(freeze [`published_signature_test_freeze.json`](../receipts/published_signature_test_freeze.json),
git head a185b41; results receipt
[`published_signature_test_results.json`](../receipts/published_signature_test_results.json)).
Run once on the analysis host in the v1 statistical runtime (Python 3.10.22, NumPy 2.2.6,
SciPy 1.15.2); the pseudobulk count matrices were hashed against their receipted pins at run
time.

## Signature sets

| Set | Source and contrast | Published genes | In the 15,435-gene family (up / down in disease) |
|---|---|---:|---|
| `polineto2020_proliferative_stage_I_II` | Poli-Neto et al. 2020 meta-analysis (GSE4888, GSE6364, GSE7305, GSE51981), proliferative phase, stage I–II endometriosis versus healthy; adjusted p < 0.05 and twofold | 4,757 | 4,661 (2,022 / 2,639) |
| `polineto2020_proliferative_stage_III_IV` | same, stage III–IV | 237 | 224 (55 / 169) |
| `bunis2022_proliferative_disease_vs_control` | Bunis et al. 2022, GSE51981 reanalysis, proliferative phase, disease versus control; direction flipped from the published control-over-disease logFC | 2,547 | 2,383 (524 / 1,859) |
| `consensus_proliferative` | genes in both of the above with the same direction (all 1,630 shared genes agree) | 1,630 | 1,599 (368 / 1,231) |

Build receipt [`published_signature_build.json`](../receipts/published_signature_build.json);
set table [`signature_sets.tsv`](published_signature_test/signature_sets.tsv); 49 unmapped
symbols in [`unmapped_or_conflicting.tsv`](published_signature_test/unmapped_or_conflicting.tsv).
The three published sets overlap in source data and all contrast disease with healthy controls.

## T1 — unsigned enrichment (pooled |t|; per-cohort secondary)

| Set | Mean \|t\| set / rest | Matched reference | Matched percentile | Wilcoxon p | Label | Per cohort (matched percentile, discovery / HECA) |
|---|---|---:|---:|---:|---|---|
| Poli-Neto I–II | 0.96 / 0.95 | 0.98 | 0.5 | 0.36 | `no_enrichment` | 3 / 67 |
| Poli-Neto III–IV | 0.80 / 0.96 | 0.98 | 0.0 | 0.9995 | `no_enrichment` | 99.8 / 99.3 |
| Bunis | 0.91 / 0.96 | 0.98 | 0.0 | 0.997 | `no_enrichment` | 13 / 30 |
| Consensus | 0.91 / 0.96 | 0.99 | 0.0 | 0.98 | `no_enrichment` | 43 / 11 |

Pooled |t| of the signature genes is at or below that of matched genes. The stage III–IV set
is enriched within each cohort separately (percentiles 99.8 and 99.3) but not in the pooled
statistic, because the two cohorts' estimates for these genes disagree in sign (see T3).

## T2 — directional concordance (primary transfer statistic)

| Set | Signed mean t (set / matched) | Matched percentile, two-sided p | Sign agreement (n) / matched | Binomial p (one-sided, > 0.5) | Spearman published logFC vs pooled t | Label |
|---|---|---|---|---:|---:|---|
| Poli-Neto I–II | **+0.246** / −0.053 | 100, 0.0002 | **58.4%** (4,661) / 48.2% | < 10⁻⁴ | 0.19 | `concordant` |
| Poli-Neto III–IV | +0.142 / +0.006 | 95.3, 0.094 | 61.2% (224) / 50.6% | 0.0005 | 0.18 | `no_detectable_concordance` |
| Bunis | **+0.110** / −0.089 | 100, 0.0002 | **54.8%** (2,383) / 47.2% | < 10⁻⁴ | 0.085 | `concordant` |
| Consensus | +0.064 / −0.108 | 100, 0.0002 | 52.8% (1,599) / 46.5% | 0.014 | 0.08 | `no_detectable_concordance` |

The matched reference share lies below 50% because the pooled stromal estimates lean positive
while most published directions are "down in disease"; the signature genes exceed that
reference by 6–11 percentage points. The stage III–IV set misses the percentile rule and the
consensus set the binomial rule; both lean the same way. The rank correlation between the
published tissue-level fold change and the pooled stromal statistic is 0.08–0.19. Per-gene
values: [`set_members.tsv`](published_signature_test/set_members.tsv).

## T3 — donor-level signature score on the v1 primary design

Score per donor = mean over signature genes of direction × centered pseudobulk log-CPM (log2
units; 1 = every gene shifted by one log2 unit in the published disease direction). Disease
coefficient from ordinary least squares on the v1 primary design (intercept, disease, age,
batch, adenomyosis; 22 donors per cohort: discovery 12 ENDO / 10 non-ENDO, HECA 15 / 7).

| Set | Cohort (genes scored) | Coefficient [95% CI] | p | MDE (80% power) | AUC | Leave-one-donor-out range | Label |
|---|---|---|---:|---:|---:|---|---|
| Poli-Neto I–II | discovery (4,709) | **−0.007** [−0.081, +0.066] | 0.83 | 0.10 | 0.56 | −0.03 to +0.02 | `no_detectable_difference` |
| | HECA (4,676) | **+0.085** [−0.055, +0.225] | 0.21 | 0.20 | 0.56 | +0.05 to +0.09 | `no_detectable_difference` |
| Poli-Neto III–IV | discovery (230) | −0.35 [−1.09, +0.40] | 0.32 | 1.04 | 0.50 | −0.51 to −0.16 | `no_detectable_difference` |
| | HECA (228) | +0.40 [−0.25, +1.05] | 0.21 | 0.91 | 0.70 | +0.12 to +0.40 | `no_detectable_difference` |
| Bunis | discovery (2,431) | −0.052 [−0.167, +0.062] | 0.33 | 0.16 | 0.59 | −0.09 to −0.03 | `no_detectable_difference` |
| | HECA (2,402) | +0.072 [−0.046, +0.191] | 0.21 | 0.17 | 0.62 | +0.04 to +0.07 | `no_detectable_difference` |
| Consensus | discovery (1,616) | −0.060 [−0.176, +0.057] | 0.28 | 0.16 | 0.58 | −0.09 to −0.04 | `no_detectable_difference` |
| | HECA (1,605) | +0.079 [−0.046, +0.203] | 0.19 | 0.17 | 0.61 | +0.04 to +0.08 | `no_detectable_difference` |

Pair label `not_replicated` for every set. The logFC-weighted projection score (sensitivity)
gives the same picture (largest set: −0.008, p 0.75, and +0.066, p 0.20). Every signature gene
was published at a twofold or larger tissue-level change, so a stromal lineage carrying the
full published effect would score about 1 or more; the observed intervals for the two large
sets lie within ±0.23. Per-donor scores: [`donor_scores.tsv`](published_signature_test/donor_scores.tsv).

## Readings

- **Direction transfers weakly, magnitude does not.** The published whole-tissue directions
  are shared by the stromal pooled estimates at 55–58% of genes, about 10 points above the
  matched expectation, with rank correlations below 0.2. At donor level the stromal score moves
  by at most a few hundredths to tenths of a log2 unit, with intervals that exclude a fifth of
  the published per-gene effect in both cohorts for the two large sets. The stroma of the v1
  donors carries a trace of the tissue signature's direction, not its size.
- **The two cohorts lean opposite ways.** Every discovery coefficient is negative and every
  HECA coefficient positive; none is significant, and the stage III–IV set shows the largest
  spread (−0.35 versus +0.40, MDE about 1). This is the same between-cohort disagreement that
  gave the per-gene combined bound no discovery.
- **Healthy versus surgical controls.** The signatures were derived against healthy controls,
  the contrast that drives the GSE51981 proliferative signal (bulk replication, E3), whereas the
  v1 cohorts use surgical controls. Whole-tissue composition (immune and epithelial
  compartments) and the control definition are the likely carriers of the published signal;
  within stroma, what remains is a weak directional echo.

## What this establishes and does not

The published bulk proliferative-phase disease signatures do not transfer in magnitude to the
stromal lineage of the two single-cell cohorts; their directions are weakly concordant with the
pooled stromal estimates. This bounds the stromal share of the whole-tissue signal at these
donor numbers; it does not locate the signal in another compartment, does not test secretory
phases, and does not evaluate the signatures' own validity. The sets overlap in source data
(GSE51981 enters all three) and contrast disease with healthy controls. The Tamaresis 2014
other-pathology contrast and the Burney 2007 lists were not available and are not tested.
Nothing here nominates a target.

## Files

- Protocol [`PUBLISHED_SIGNATURE_PROTOCOL.md`](../docs/PUBLISHED_SIGNATURE_PROTOCOL.md);
  configuration [`published_signature_test.json`](../config/published_signature_test.json).
- Code: [`published_signature_build.py`](../scripts/published_signature_build.py) (sets from the
  downloaded supplements), [`published_signature_test_analysis.py`](../scripts/published_signature_test_analysis.py)
  (T1–T3); tests [`test_published_signature_test.py`](../tests/test_published_signature_test.py)
  (planted signature recovered by all three tests, random set unlabelled, flipped directions
  labelled discordant, freeze and once-only refusals).
- Receipts: supplement downloads in [`literature_receipts.tsv`](../manifests/literature_receipts.tsv)
  (route `europepmc_supplementary`), build, freeze and results receipts as linked above.
- Outputs in [`published_signature_test/`](published_signature_test/): `summary.json`,
  `set_results.tsv`, `set_members.tsv`, `donor_scores.tsv`, the set table and the unmapped
  list. Publisher supplements and pseudobulk counts are not committed.
