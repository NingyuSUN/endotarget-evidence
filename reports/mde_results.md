# Minimum detectable effect analysis — results

**What the v1 null means in numbers.** For the median tested gene, the v1 disease design had
80% power for a true difference of about 1.5-fold as a single prespecified test, but needed
about 2.2-fold for that gene to appear among 100 BH discoveries and about 3.3-fold to be the
only discovery. The observed intervals exclude twofold disease differences for more than
half of all tested genes even under a family-wide (Bonferroni) interval, and exclude 1.5-fold
differences for about half of the genes per gene but for almost none family-wide. A study of
22 donors recovers a program as large as the secretory-phase program inconsistently: at
11 + 11 donors with the primary design, the median draw detected 45% of the discovery phase
program and 13% of the HECA phase program, and the HECA 2.5th percentile was 3 genes.

Completed 2026-10-27 under the protocol frozen the same day
([protocol](../docs/MDE_ANALYSIS_PROTOCOL.md); freeze
[`mde_analysis_freeze.json`](../receipts/mde_analysis_freeze.json); results receipt
[`mde_analysis_results.json`](../receipts/mde_analysis_results.json); published products
[`reports/mde_analysis/`](mde_analysis/)). No disease model was fitted; the v1 and
positive-control artifacts were hashed before and after the run and are unchanged. The
first freeze (commit 5005d76) was superseded before any value was computed, because the
host run stopped in its verification step on a defect in the verifier; the re-issued
freeze records this.

## Part A — minimum detectable effect from the saved v1 fits

MDE is the smallest true log2 fold change that the v1 moderated test would have rejected
with 80% probability, at the realized precision of each gene (moderated SE and total
degrees of freedom), solved exactly on the non-central t distribution. Scenarios:
`alpha_0.05` is a single prespecified gene; `bh_k` is the Benjamini–Hochberg boundary when
exactly *k* genes are rejected in a family of about 17,700.

| Cohort (disease contrast) | Scenario | Per-gene α | MDE median (fold) | MDE q25 – q75 | Genes with MDE ≤ 1 | ≤ 2 |
|---|---|---|---|---|---|---|
| Discovery, 17,733 genes, residual df 10 | `alpha_0.05` | 0.05 | **0.57 (1.5×)** | 0.36 – 1.06 | 73% | 95% |
| | `bh_k1000` | 2.8 × 10⁻³ | 0.87 (1.8×) | 0.56 – 1.62 | 56% | 83% |
| | `bh_k100` | 2.8 × 10⁻⁴ | **1.12 (2.2×)** | 0.72 – 2.09 | 45% | 73% |
| | `bh_k10` | 2.8 × 10⁻⁵ | 1.40 (2.6×) | 0.90 – 2.63 | 32% | 64% |
| | `bh_k1` | 2.8 × 10⁻⁶ | **1.73 (3.3×)** | 1.11 – 3.24 | 19% | 56% |
| HECA, 17,643 genes, residual df 12 | `alpha_0.05` | 0.05 | **0.63 (1.5×)** | 0.41 – 0.97 | 76% | 97% |
| | `bh_k1000` | 2.8 × 10⁻³ | 0.94 (1.9×) | 0.61 – 1.46 | 53% | 88% |
| | `bh_k100` | 2.8 × 10⁻⁴ | **1.20 (2.3×)** | 0.77 – 1.86 | 40% | 79% |
| | `bh_k10` | 2.8 × 10⁻⁵ | 1.47 (2.8×) | 0.95 – 2.28 | 28% | 68% |
| | `bh_k1` | 2.8 × 10⁻⁶ | **1.78 (3.4×)** | 1.15 – 2.76 | 17% | 57% |

Precision varies about eightfold across genes (moderated SE 5th to 95th percentile 0.08 to
0.67 in discovery, 0.09 to 0.58 in HECA), and it is driven by expression: in the lowest
expression quintile the discovery median MDE is 1.30 (`alpha_0.05`) and 2.58 (`bh_k100`),
in the highest 0.36 and 0.71. The robust prior contributes only 2.6 to 3.1 degrees of
freedom to the residual 10 and 12, so the per-gene tests rest on 12.6 to 15.1 total degrees
of freedom. The closed-form approximation (t<sub>crit</sub> + t<sub>0.8</sub>) · SE agrees
with the exact solution within 0.005 log2 at `alpha_0.05` but departs by up to 2.5 log2 for
the least precise genes at the smallest α; the exact values are reported throughout.

**What the observed intervals exclude.** Genes whose confidence interval for the disease
coefficient lies entirely within ±Δ log2:

| Cohort | Interval | ±0.585 (1.5×) | ±1 (2×) | ±2 (4×) |
|---|---|---:|---:|---:|
| Discovery (17,733) | per-gene 95% | 8,632 (49%) | 12,717 (72%) | 16,687 (94%) |
| | family-wide Bonferroni 95% | 123 (0.7%) | 3,734 (21%) | 10,145 (57%) |
| HECA (17,643) | per-gene 95% | 8,356 (47%) | 13,788 (78%) | 17,150 (97%) |
| | family-wide Bonferroni 95% | 77 (0.4%) | 3,584 (20%) | 10,610 (60%) |

The per-gene counts are not simultaneous statements; the family-wide counts are conservative.
Between them: the v1 data rule out large (fourfold) proliferative-phase stromal disease
differences for most genes under any reading, rule out twofold differences for most genes
per gene and for a fifth of genes family-wide, and say little about 1.5-fold differences
family-wide. Nothing here says that effects below these bounds are absent.

**The phase fits under the same lens.** The positive-control fits (35 and 42 donors, residual
df 22 and 29) have median MDE 0.50 and 0.52 at `alpha_0.05` and 0.89 and 0.90 at `bh_k100`.
Their detected programs have median |log2FC| 0.79 (discovery, 7,367 genes) and 0.84 (HECA,
3,986 genes), with quartiles 0.45 to 1.50 and 0.50 to 1.39.

**Bridge (descriptive).** Share of the phase program whose |log2FC| exceeds the disease fit's
median MDE, and, gene by gene, that gene's own disease MDE:

| Cohort | Scenario | above the disease median MDE | above the gene's own disease MDE (gene-matched n) |
|---|---|---:|---:|
| Discovery | `alpha_0.05` | 65% | 66% (7,291) |
| | `bh_k100` | 35% | 18% |
| | `bh_k1` | 20% | 4.6% |
| HECA | `alpha_0.05` | 65% | 63% (3,936) |
| | `bh_k100` | 32% | 14% |
| | `bh_k1` | 14% | 3.0% |

A program of the secretory phase's size and shape would, gene by gene, clear the v1 single-
discovery boundary for only 3 to 5% of its genes and the 100-discovery boundary for 14 to
18%. This compares two different contrasts on overlapping donors; it does not convert phase
effects into expectations about disease effects.

## Part B — donor-subsampling sensitivity curve on the phase signal

Each draw takes *n* secretory and *n* proliferative donors without replacement from the
positive-control cohorts and is analyzed as a new study with the v1 estimator (filter, TMM,
voom and robust empirical Bayes re-estimated; BH over the draw's family), under the primary
design (with batch contrasts) and a reduced design (without). 100 draws per level; identity
check of the subsampling estimator against the frozen full-cohort fits: maximum difference
5.6 × 10⁻¹⁴ in log2FC (discovery) and exactly zero (HECA). Median over valid draws with the
2.5th to 97.5th percentile in brackets; recall is the share of the full-cohort q < 0.05 set
recovered with the same sign; sign disagreement with the full-cohort estimate was ≤ 0.13%
at every level.

**Discovery** (full cohort 13 + 22 donors: 7,367 genes at q < 0.05):

| n per arm | Design | Valid draws | Genes at q < 0.05 | Recall | Top-20 overlap | Typical MDE `alpha_0.05` / `bh_k100` |
|---:|---|---:|---|---|---:|---|
| 5 | primary | 0 / 100 | — | — | — | — |
| 5 | reduced | 100 | 1,510 [41, 4,013] | 0.19 [0.01, 0.47] | 5 | 1.03 / 2.29 |
| 6 | primary | 0 / 100 | — | — | — | — |
| 6 | reduced | 100 | 2,519 [107, 5,108] | 0.32 [0.01, 0.59] | 6.5 | 0.88 / 1.84 |
| 8 | primary | 30 / 100 | 1,310 [81, 4,013] | 0.17 [0.01, 0.50] | 7 | 0.98 / 2.21 |
| 8 | reduced | 100 | 3,893 [1,054, 5,989] | 0.50 [0.14, 0.69] | 9 | 0.74 / 1.43 |
| 10 | primary | 100 | 2,539 [415, 4,832] | 0.34 [0.06, 0.62] | 10 | 0.81 / 1.69 |
| 10 | reduced | 100 | 4,987 [3,357, 6,533] | 0.64 [0.44, 0.78] | 11 | 0.64 / 1.17 |
| **11** | **primary** | 100 | **3,386 [1,346, 5,544]** | **0.45 [0.18, 0.70]** | 11.5 | 0.74 / 1.48 |
| 11 | reduced | 100 | 5,649 [3,934, 6,948] | 0.71 [0.52, 0.83] | 11 | 0.60 / 1.08 |
| 13 | primary | 100 | 4,705 [3,484, 5,770] | 0.63 [0.47, 0.77] | 14 | 0.64 / 1.19 |
| 13 | reduced | 100 | 6,444 [5,386, 7,492] | 0.80 [0.70, 0.86] | 12 | 0.54 / 0.96 |

**HECA** (full cohort 20 + 22 donors: 3,986 genes at q < 0.05):

| n per arm | Design | Valid draws | Genes at q < 0.05 | Recall | Top-20 overlap | Typical MDE `alpha_0.05` / `bh_k100` |
|---:|---|---:|---|---|---:|---|
| 5 | primary | 0 / 100 | — | — | — | — |
| 5 | reduced | 100 | 26 [0, 1,310] | 0.006 [0, 0.23] | 3 | 1.08 / 2.46 |
| 6 | primary | 0 / 100 | — | — | — | — |
| 6 | reduced | 100 | 98 [0, 1,961] | 0.02 [0, 0.35] | 3 | 0.95 / 2.00 |
| 8 | primary | 66 / 100 | 70 [0, 1,263] | 0.02 [0, 0.20] | 4.5 | 1.09 / 2.50 |
| 8 | reduced | 100 | 792 [40, 2,301] | 0.17 [0.01, 0.44] | 7 | 0.76 / 1.45 |
| 10 | primary | 96 / 100 | 235 [0, 2,545] | 0.05 [0, 0.48] | 6 | 0.95 / 1.97 |
| 10 | reduced | 100 | 1,339 [359, 3,488] | 0.30 [0.08, 0.66] | 9 | 0.66 / 1.22 |
| **11** | **primary** | 98 / 100 | **531 [3, 2,069]** | **0.13 [0.001, 0.46]** | 7 | 0.84 / 1.64 |
| 11 | reduced | 100 | 1,565 [513, 3,350] | 0.34 [0.11, 0.62] | 10 | 0.64 / 1.16 |
| 13 | primary | 100 | 1,215 [224, 3,259] | 0.29 [0.06, 0.67] | 10 | 0.73 / 1.36 |
| 13 | reduced | 100 | 2,362 [1,234, 4,242] | 0.49 [0.29, 0.73] | 12 | 0.58 / 1.02 |
| 16 | primary | 100 | 2,374 [951, 4,297] | 0.55 [0.22, 0.82] | 14 | 0.64 / 1.14 |
| 16 | reduced | 100 | 3,301 [2,161, 4,342] | 0.67 [0.43, 0.80] | 12.5 | 0.51 / 0.88 |
| 20 | primary | 100 | 3,682 [2,938, 4,134] | 0.90 [0.72, 0.96] | 18 | 0.54 / 0.93 |
| 20 | reduced | 100 | 4,180 [3,901, 4,528] | 0.79 [0.75, 0.83] | 14 | 0.46 / 0.78 |

Invalid draws under the primary design were all structural: with 10 or 12 donors spread over
9 capture batches or 9 library lanes, the batch contrasts were rank-deficient or left fewer
than 5 residual degrees of freedom in every draw (discovery n = 5: 58 rank-deficient, 39
below the df minimum, 3 without a batch containing both phases); at n = 8 the primary design
was estimable in 30 of 100 discovery draws and 66 of 100 HECA draws. These rates are a
property of the design, not of the data drawn. Figures:
[discovery](mde_analysis/curve_discovery.png), [HECA](mde_analysis/curve_heca.png); MDE
distributions per fit: [v1 discovery](mde_analysis/v1_discovery_mde_ecdf.png),
[v1 HECA](mde_analysis/v1_heca_mde_ecdf.png).

Three readings of the curves, all descriptive:

- **Detection of a strong program is far from saturated at 22 donors.** At 11 + 11 with the
  primary design, the median draw recovers 45% (discovery) and 13% (HECA) of the full-cohort
  program, with 2.5th percentiles of 18% and 0.1%; HECA draws at that size detected between
  3 and 2,069 genes. The 2.5th-percentile HECA draw at 11 + 11 is nearly a null result on a
  program of 3,986 genes.
- **The covariate structure costs as much as several donors.** At every level the reduced
  design detects more and recovers more than the primary design; in HECA at 20 + 20 the
  reduced design detects more genes than the full 42-donor primary fit, which is the price
  of eight lane contrasts. The protocol does not change the primary design on this basis;
  it records the cost.
- **Variability between draws of the same size is large.** The 95% ranges of recall span 0.5
  or more at every level up to 13 per arm. Any single small study's gene count is a weak
  guide to the program's size.

## What this establishes and what it does not

The v1 disease analysis was sensitive to 1.5-fold effects per gene and to 2.2- to 3.4-fold
effects at the multiplicity it actually faced, with wide gene-to-gene variation set by
expression level; its intervals exclude fourfold proliferative-phase stromal disease
differences for most genes and twofold differences for a fifth of genes family-wide. A
22-donor study of this design recovers even a several-thousand-gene program inconsistently
and sometimes barely at all. These are statements about the design's resolution and about
one known program's behaviour under subsampling. They do not estimate disease effects, do
not say any specific gene's effect is small, do not rank anything, and do not choose a
future sample size; they are the quantities such a choice would use, computed once under
rules fixed in advance. The MDE values are conditional on the realized variance estimates
and precision weights of the original fits, and the subsampling levels are limited to the
donors available (13 secretory in discovery, 20 in HECA).

## Execution record

Analysis host, v1 statistical runtime (R 4.4.3, edgeR 4.4.0, limma 3.62.1). Part A: four
fits, about 40 s each. Part B: 600 and 800 draws, two designs each, 12 workers per cohort,
11.5 min in total. Protected artifacts (both v1 fit directories, both phase fit directories,
both positive-control count directories, 47 files): hashes identical before and after.
Untracked working copies were removed from the host checkout afterwards; results stay under
the host's `work/mde_analysis/`. Published copies and receipts are listed in
[`mde_analysis_results.json`](../receipts/mde_analysis_results.json).
