# External evaluation results — English summary

This page summarizes the generated bilingual report
[replication_v2_results.md](replication_v2_results.md) in English. It adds no new
analysis: every number below is copied from that report, its receipts or the
[research brief](executive_brief_en.md). The generated report remains the
authoritative record and lists all traceable inputs.

**Main result.** Of the frozen top-20 lists, expression ranking (B0) had 14
same-direction genes among 20 evaluable and 20 selected; network propagation (G1)
had 6 / 19 / 20. In the external cohort, 0 of 17,643 tested genes reached
BH q < 0.05. These results do not confirm a therapeutic target or an intervention
direction.

## Question

Within author-adjudicated proliferative-phase eutopic endometrial stromal cells,
do rankings frozen in the discovery cohort — expression (B0), network degree (B1)
and network propagation (G1) — preserve the direction of expression change in an
independent cohort? Counts are summed per donor and each person is one
statistical unit. Rankings, directions, top-20/50 lists and matched reference sets
were frozen before any external disease effect was computed.

External cohort: 15 cases and 7 non-endometriosis controls (22 donors); 17,643
genes entered the full tested family and none reached BH q < 0.05. The smallest
unadjusted p was 6.993 × 10⁻⁵ and the smallest BH q was 0.9699. Each cohort was
modeled separately with library-size/composition normalization and donor-level
weighted fits, adjusting for the frozen age, batch and adenomyosis covariates.
A log2 difference is an association that can include cell-state composition and
unmeasured confounding; controls can have other gynecological conditions.

## Fixed-list results

| Method | k | Same direction | Evaluable | Selected | Same direction among evaluable | Missing-outcome bounds | Same direction and full-family q < 0.05 |
|---|---:|---:|---:|---:|---:|---:|---:|
| B0 | 20 | 14 | 20 | 20 | 70.0% | 70.0%–70.0% | 0 |
| B0 | 50 | 36 | 49 | 50 | 73.5% | 72.0%–74.0% | 0 |
| B1 | 20 | 3 | 13 | 20 | 23.1% | 15.0%–50.0% | 0 |
| B1 | 50 | 16 | 35 | 50 | 45.7% | 32.0%–62.0% | 0 |
| G1 | 20 | 6 | 19 | 20 | 31.6% | 30.0%–35.0% | 0 |
| G1 | 50 | 26 | 47 | 50 | 55.3% | 52.0%–58.0% | 0 |

Direction agreement only asks whether evaluable genes keep the same sign; small,
imprecise effects can agree by chance. The missing-outcome bounds treat every
unevaluable gene as all-disagreeing or all-agreeing; lists are never refilled and
missing genes are never counted as opposite. Top 20 is the primary descriptive
endpoint and top 50 is secondary. Seven B1 genes and one G1 gene failed external
expression filtering. B1 and G1 share 9 of their top-20 genes.

G1 minus B0 (top 20) was −0.384, with missing-outcome bounds [−0.40, −0.35]. Both
methods needed at least 16/20 evaluable genes, and no advantage could be claimed
if the bounds included zero or could reverse. No evidence met all predefined
conditions for G1 outperforming B0.

![Same-direction, opposite-direction and unevaluable genes in the fixed top-20 and top-50 lists](figures/replication_v2_counts.png)

## Matched reference sets

| Method | Status | Adequately covered sets | Observed percentile | Observed list coverage |
|---|---|---:|---:|---:|
| B0 | ESTIMABLE | 200/200 | 96.5 | 20/20 |
| B1 | ESTIMABLE | 200/200 | 0.5 | 13/20 |
| G1 | ESTIMABLE | 200/200 | 0.5 | 19/20 |

Each method has 200 frozen sets matched on expression and degree, with directions
taken from each gene's own discovery effect. The percentile is a descriptive
position; the sets share genes and donors, so they are neither 200 independent
biological replicates nor a calibrated p value.

## Donor sensitivity

Paired external-donor resampling completed 500 of 500 draws: 487 valid and 13
failed (design not full rank). All methods share the same donor draws, fixed
discovery lists/directions and original evaluable denominators; failures are not
replaced. The conditional 95% interval for the G1 − B0 top-20 difference was
[−0.390, +0.087] (type-7 percentile of valid draws). It describes external
donor-sampling variation for the already selected lists only.

All 22 discovery leave-one-donor-out refits completed (0 failures). Each refit
repeats filtering, model fitting and ranking without one donor.

| Method | Valid / total folds | Median top-20 retention [range] | Median Jaccard [range] |
|---|---:|---:|---:|
| B0 | 22/22 | 0.725 [0.3, 1] | 0.569 [0.176, 1] |
| B1 | 22/22 | 1 [1, 1] | 1 [1, 1] |
| G1 | 22/22 | 0.975 [0.8, 1] | 0.952 [0.667, 1] |

Retention is the intersection divided by the original 20; Jaccard is the
intersection divided by the union. Stable network lists did not yield better
external sign agreement.

![Paired bootstrap distribution of the G1 − B0 top-20 agreement difference](figures/replication_v2_bootstrap.png)

![Leave-one-donor-out list stability](figures/replication_v2_lodo.png)

![Discovery versus external expression effects](figures/replication_v2_scatter.png)

## Notes for interpretation

These observations follow directly from the reported numbers; they were not
prespecified and have not been tested further.

- **The bootstrap interval is not centered on the observed difference.** The
  observed G1 − B0 value (−0.384) lies at the lower edge of the interval
  [−0.390, +0.087], and most valid draws fall near −0.2. Resampling external
  donors adds noise to each gene's external sign, which pulls both methods'
  agreement toward 50% and shrinks their difference. The interval therefore
  overstates how often the difference could reach zero. This does not change the
  conclusion: G1 showed no advantage, and its observed agreement was lower than
  B0's.
- **B1 and G1 agreed less often than their matched reference sets.** Both
  network lists sit at the 0.5th percentile of their 200 matched sets, and 13 of
  19 evaluable G1 genes changed direction. This is a descriptive
  lower-than-reference pattern, not only an absence of improvement; the
  percentile is not a calibrated p value. Possible contributors, such as systematic
  single-cell versus single-nucleus differences for highly connected genes or
  cell-state composition, were not investigated.

## Candidate pool and dossiers

The union of the frozen B0 and G1 top-20 lists contains 40 genes, all retained in
the supplementary tables. The deterministic rule allowed up to three
same-direction genes with full-family q < 0.05, one precise opposite case and one
missing/inconclusive case, each chosen by full B0 rank and stable ID. It selected
one dossier and zero supported same-direction dossiers.

**PARP9 (ENSG00000138496).** Expression inconclusive: discovery log2 difference
+0.4659; external −0.1415, unadjusted 95% interval [−0.4126, 0.1296], full-family
q = 0.9699. No mapping in the audited genetics source. Cell-state specificity and
intervention direction are not established. Recommendation: **hold** pending
better expression precision or missing evidence; no drug or clinical
recommendation.

Genetic evidence comes from one prespecified 2023 publication and its public
mapping table; several annotations from that paper are not independent
confirmation, and a missing mapping is not a negative genetic test. See the
[genetic evidence report](genetic_evidence_report.md).

## Limits

- The 30/50/100-cell thresholds select the same 22 + 22 donors, libraries and
  designs, so the primary fits are reused; this is not three independent
  validations. All 44 donors have complete ages.
- Discovery reuses author-provided secondary Alevin estimated (fractional) counts;
  the original full RDS was not recomputed. In HECA, raw and stored H5 totals for
  52,744 nuclei differ by 9,607 / 360,765,431 counts (0.00266%): 1,359 are
  explained by globally missing features and 8,248 remain unassigned.
- Inherited phase and cell labels with marker checks support group-level
  plausibility, not purity or independent re-annotation. Separate recruitment
  sources do not establish genotype-confirmed non-overlap of individuals.
- Clinical subset analyses are reported separately in the
  [clinical sensitivity report](clinical_sensitivity_report.md) (Chinese); the
  [research brief](executive_brief_en.md) summarizes them in English.
