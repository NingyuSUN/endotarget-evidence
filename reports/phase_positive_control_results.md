# Phase positive control — results and decision

**Decision: the v1 pipeline detects and transfers a known stromal expression program at these
sample sizes. Both primary endpoints pass. The v1 disease null therefore cannot be attributed
to an insensitive pipeline; it reflects disease effects below what 22 donors per cohort
resolve in proliferative stromal fibroblasts. No gene becomes a target because of this
result, and nothing here is disease evidence.**

Completed 2026-10-27 under the protocol frozen on 2026-10-26
([protocol](../docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md); freezes
[F0](../receipts/phase_positive_control_freeze.json),
[F1](../receipts/phase_positive_control_f1_discovery_freeze.json); evaluation receipt
[`phase_positive_control_evaluation.json`](../receipts/phase_positive_control_evaluation.json)).
The estimator, network, ranking code and external-evaluation definitions are the v1
implementations unchanged; the only prespecified differences are the contrast
(`SEC_minus_PRO`) and cycle phase as the grouping variable. Executed on the native analysis
host (SRV01) in the v1 statistical runtime.

## Cohorts and fits

| Cohort | Donors (secretory / proliferative) | Primary cells or nuclei | Genes tested | BH q < 0.05 |
|---|---|---:|---:|---:|
| Discovery, GSE266265 (scRNA) | 13 / 22 | 225,207 | 17,341 | **7,367** (3,611 up, 3,756 down) |
| External, HECA E-MTAB-14039 (snRNA) | 20 / 22 | 94,150 | 17,243 | **3,986** (1,825 up, 2,161 down) |

Designs: phase, disease, centered age, batch (capture batch; library lane), adenomyosis;
rank 13 in both, residual df 22 and 29. The proliferative arms are exactly the v1 donors
([discovery counts](phase_positive_control_discovery_counts.md),
[HECA counts](phase_positive_control_heca_counts.md)). Compare v1's disease fits on the same
pipeline: 0 genes at q < 0.05 in each cohort.

## Endpoints against the prespecified thresholds

| Endpoint | Rule | Observed | Result |
|---|---|---|---|
| **P1** within-cohort detection | ≥100 genes at q < 0.05 in each cohort | 7,367 and 3,986 | **pass** (both above the 1,000 "strong" label) |
| **P2** cross-cohort transfer of the frozen B0 top-20 | ≥16/20 same direction in HECA and ≥10/20 also HECA q < 0.05; ≥16 evaluable | **20/20 same direction, 20/20 at q < 0.05**; missing-outcome bounds 1.00–1.00 | **pass** |
| S1 marker panel (explanatory) | ≥6/7 expected sign with q < 0.05 in discovery; ≥6/7 expected sign in HECA | discovery 5/7; HECA 6/7 | discovery component fails, HECA passes: **not met**, as reported at F1 |
| S2 network baselines (optional) | descriptive; no superiority claim unless the paired interval excludes zero with adequate coverage | G1 − B0 = −0.37 in the fixed sample; 499/500 valid draws; conditional 95% interval **[−0.53, −0.21]** | G1 below B0 across the interval |

Secondary top-50: B0 49/49 evaluable genes same direction, all at HECA q < 0.05 (one gene
filtered); G1 36/48 (13 significant); B1 19/34 (16 filtered, coverage-limited).

**The frozen B0 top-20 in HECA** (discovery sign; HECA log2 FC; HECA q): BCAT1 (+; +2.19;
1.4 × 10⁻⁷), CFD (+; +3.87; 1.2 × 10⁻³), ISOC1 (−; −1.80; 7.1 × 10⁻⁸), EVA1C (+; +2.33;
2.3 × 10⁻¹⁰), MAOA (+; +3.59; 8.0 × 10⁻⁶), PDLIM1 (+; +1.96; 1.3 × 10⁻⁷), APOD (+; +2.54;
8.6 × 10⁻⁴), ZBTB16 (+; +6.55; 5.0 × 10⁻⁶), TRIB2 (−; −2.08; 5.5 × 10⁻⁸), SCARA5 (+; +4.82;
2.5 × 10⁻⁴), ENOX2 (+; +1.29; 1.2 × 10⁻⁸), MTHFD2 (+; +1.84; 5.3 × 10⁻⁶), MYOCD (+; +3.76;
1.0 × 10⁻³), CRYAB (+; +2.41; 1.8 × 10⁻³), ZEB1 (+; +1.56; 4.0 × 10⁻⁹), MASP1 (+; +2.90;
2.8 × 10⁻⁸), IL1R1 (+; +1.28; 3.2 × 10⁻⁴), FMO1 (−; −3.34; 5.7 × 10⁻⁵), DKK1 (+; +4.58;
3.0 × 10⁻⁵), ALDH1A3 (+; +4.88; 7.5 × 10⁻⁸). Every gene keeps its direction across platform
(scRNA to snRNA), study, country and label vocabulary.

## Network methods under a strong signal

| Top-20 list | Evaluable | Same direction | Same direction and HECA q < 0.05 |
|---|---:|---:|---:|
| B0 expression | 20 | 20 | 20 |
| B1 degree | 13 (7 hubs fail HECA expression filter) | 8 | 2 |
| G1 propagation | 19 | 12 | 4 |

v1 observed G1 − B0 = −0.38 with a conditional interval of [−0.39, +0.09] that included zero;
with a strong program to propagate, the same comparison gives −0.37 with an interval of
[−0.53, −0.21] that excludes zero. Network propagation on this frozen functional-association
graph ranks hub genes (TP53, UBC, ACTB, histone genes) whose phase effects are small, and
those lists transfer worse than the expression ranking. This is a statement about this
graph, this seeding rule and this transfer endpoint, not about graph methods in general;
the bootstrap conditions on the fixed lists and excludes discovery-selection uncertainty.
One of 500 draws was invalid (rank-deficient resampled design) and was retained, not
redrawn.

## Marker panel detail

HECA: six of six fitted genes carry the expected sign (DKK1, PLCL1, ESR1, PGR at q < 0.05;
IGFBP1 q = 0.063; LEFTY2 q = 0.27); PRL fell out of the HECA expression filter. Discovery:
six of six fitted genes carry the expected sign, five at q < 0.05; IGFBP1 fell out of the
discovery filter. The two filtered genes and the two imprecise ones are the panel's
late-secretory and secreted-protein members, in cohorts whose secretory arms are mostly early
and mid secretory and, for HECA, measured in nuclei. The panel was prespecified as
explanatory and its thresholds were not revisited.

## What this does and does not establish

- The pipeline, as frozen for v1, is sensitive to a large donor-level program and transfers
  it across the two cohorts. The v1 disease null is informative: it bounds the effect size
  of proliferative-phase stromal disease differences below what this design detects.
- Secretory-phase genes are not endometriosis genes. The B0 list is a positive-control
  input and well-known cycle biology, not a shortlist of any kind.
- The lineage contrast includes shifts among stromal states (decidualized versus
  non-decidualized) by design; it is not a within-state effect.
- HECA's secretory arm is 80% endometriosis cases; disease is adjusted additively, and the
  phase coefficient remains an association.
- Everything was run once, under frozen rules, with every intermediate hashed; none of
  the v1 artifacts was modified (their hashes are recorded in the receipts).

## Files

Published copies of the frozen and evaluated products:
[`reports/phase_positive_control_f1/`](phase_positive_control_f1/) (discovery fit, rankings,
markers) and [`reports/phase_positive_control_f2/`](phase_positive_control_f2/) (HECA fit,
evaluation tables, HECA markers, bootstrap draws and summary). Preceding reports:
[feasibility](phase_positive_control_feasibility.md),
[HECA counts](phase_positive_control_heca_counts.md),
[discovery counts](phase_positive_control_discovery_counts.md),
[discovery fit and F1](phase_positive_control_f1_discovery.md).
