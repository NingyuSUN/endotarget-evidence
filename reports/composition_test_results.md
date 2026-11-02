# Stromal-state composition and disease — results

**Labels under the frozen rules: `no_detectable_shift` in both cohorts; pair
`no_replicated_shift`.** The share of decidualized-like states among a donor's primary
stromal fibroblasts does not differ detectably between endometriosis and non-endometriosis
donors in the proliferative phase: odds ratio 0.75 (95% interval 0.48–1.17, p 0.18) in
discovery and 1.08 (0.64–1.80, p 0.75) in HECA, with the v1 covariates. The designs would
have detected an odds ratio of about 1.9 (or 1/1.9) in discovery and 2.1 in HECA with 80%
power. No author state shifts at BH q < 0.05 within either cohort (smallest q 0.16 and
0.21).

Completed 2026-11-02 under the [protocol](../docs/COMPOSITION_TEST_PROTOCOL.md) frozen the
same day before any model was fitted on the composition tables (freeze
[`composition_test_freeze.json`](../receipts/composition_test_freeze.json), git head aae7faf;
results receipt [`composition_test_results.json`](../receipts/composition_test_results.json)).
Inputs: the published composition grids of the [composition audit](stromal_composition.md)
(22 donors × 12 author states in discovery, 22 × 7 in HECA, complete with observed zeros)
and the v1 primary design matrices, pinned by hash. Run once locally (Python 3.10.21, NumPy
2.2.6, SciPy 1.15.2). The coarse partition (decidualized-like `dS*` / `dStromal_*`,
non-decidualized `eF*` / `eStromal*`, other) follows the authors' label prefixes and was fixed
before any value was read.

## Primary endpoint: decidualized-like share

| Cohort | Donors (ENDO / non-ENDO) | Mean share, ENDO | Mean share, non-ENDO | Disease coefficient (logit) | Odds ratio [95%] | p (df) | Detectable odds ratio (80% power) | Label |
|---|---|---:|---:|---:|---|---|---|---|
| Discovery (scRNA) | 12 / 10 | 12.8% | 13.2% | −0.29 (SE 0.20) | 0.75 [0.48, 1.17] | 0.18 (10) | 1.86 | `no_detectable_shift` |
| HECA (snRNA) | 15 / 7 | 7.6% | 7.4% | +0.08 (SE 0.24) | 1.08 [0.64, 1.80] | 0.75 (12) | 2.06 | `no_detectable_shift` |

Model: logit of (decidualized-like cells + 0.5) / (other primary cells + 0.5), ordinary least
squares on the v1 primary design (disease, centered age, batch contrasts, adenomyosis), one
observation per donor.

## Secondary endpoints

| Cohort | Non-decidualized share, logit coefficient [95%], p | "Other" share, coefficient [95%], p |
|---|---|---|
| Discovery | +0.41 [−0.13, +0.94], 0.12 | −0.90 [−1.90, +0.11], 0.076 |
| HECA | −0.04 [−0.52, +0.44], 0.85 | −0.38 [−1.97, +1.21], 0.61 |

Per author state (centered log-ratio on the same design; BH within cohort over 12 and 7
states; full table in [`per_state_clr.tsv`](composition_test/per_state_clr.tsv)): the smallest
adjusted values are discovery `fib C7` (coefficient −0.69, p 0.013, q 0.16; mean share 3.6%
in ENDO versus 4.9%) and HECA `dStromal_mid` (+0.51, p 0.031, q 0.21; 2.3% versus 1.7%).
Nothing reaches q < 0.05. The largest raw difference in group means, HECA `eStromal` 57%
versus 77% with `eStromal_MMPs` 24% versus 2%, has p 0.40 and 0.48: it is carried by a few
donors and the donor-level variance is wide, which is exactly why a donor-level test is
needed.

## Sensitivities (none changes a label)

| Sensitivity | Discovery coefficient [95%], p | HECA coefficient [95%], p |
|---|---|---|
| Arcsine square-root share | −0.045 [−0.119, +0.029], 0.20 | +0.009 [−0.065, +0.082], 0.81 |
| Reduced design (disease, age) | −0.11 [−0.35, +0.14], 0.38 | −0.11 [−0.54, +0.33], 0.61 |
| Donors with < 500 primary cells excluded | identical (none excluded) | +0.19 [−0.37, +0.74], 0.47 (2 excluded) |
| Adenomyosis donors excluded | −0.43 [−0.86, +0.01], 0.053 (4 excluded, df 7) | −0.003 [−0.55, +0.54], 0.99 (3 excluded) |
| Leave-one-donor-out (22 refits) | −0.46 to −0.11, sign constant, smallest p 0.02 | −0.03 to +0.17, sign not constant, smallest p 0.49 |

The discovery estimate is consistently negative (fewer decidualized-like cells in
endometriosis donors) across refits and sharpens when adenomyosis donors are excluded, as the
secretory-program analysis also saw for its own estimand; it does not reach the per-cohort
threshold, and HECA shows nothing in the same direction. It is a hypothesis for a larger
cohort, not a finding.

## What this establishes and does not

With the same donors and covariates as v1, disease does not detectably change the mixture of
author-labelled stromal states within the primary lineage; the data are compatible with
odds ratios for the decidualized-like share between about 0.5 and 1.2 (discovery) and 0.6 and
1.8 (HECA). This closes the composition half of the objection that pseudobulk hides a
within-lineage change at these donor numbers; it cannot address within-state expression
changes, inherits author labels with different vocabularies and resolutions, has no positive
control of its own, and rests on 22 donors per cohort. No state and no gene is nominated.

## Files

[`reports/composition_test/`](composition_test/): [donor shares](composition_test/donor_shares.tsv),
[per-state CLR table](composition_test/per_state_clr.tsv), [summary](composition_test/summary.json).
Code: [`scripts/composition_test_analysis.py`](../scripts/composition_test_analysis.py),
[`scripts/followon_freeze.py`](../scripts/followon_freeze.py), [tests](../tests/test_composition_test.py).
