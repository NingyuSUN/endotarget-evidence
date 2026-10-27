# Minimum detectable effect analysis — prespecified protocol

Version 1.0, 2026-10-27. Status: frozen before any minimum-detectable-effect value was
computed (receipt `receipts/mde_analysis_freeze.json`). Machine-readable parameters:
[`config/mde_analysis.json`](../config/mde_analysis.json). This protocol fixes every
definition, threshold scenario, resampling rule and interpretation sentence below; nothing
is chosen after seeing a result.

## 1. Purpose and scope

The v1 disease analysis ([method card](METHOD_CARD.md)) found no gene at BH q < 0.05 in
either cohort, and the [phase positive control](PHASE_POSITIVE_CONTROL_PROTOCOL.md) showed
that the same estimator detects and transfers a large stromal program at the same donor
numbers. This analysis turns the null into quantitative statements, in two parts:

- **Part A — minimum detectable effect (MDE) from the saved v1 fits.** For every tested
  gene, the smallest true log2 fold change that the v1 model would have detected with 80%
  probability, under the realized precision of that gene, at several significance scenarios.
  The same computation is applied to the two phase positive-control fits as a reference.
- **Part B — donor-subsampling sensitivity curve on the phase signal.** How the number of
  detected genes, the recall of the full-cohort phase program and the typical MDE change
  when a smaller study is drawn from the positive-control cohorts and analyzed from scratch
  with the v1 estimator.

The analysis fits **no disease model**, modifies **no v1 or positive-control artifact**, and
produces **no evidence for or against any gene**. It bounds what the v1 design could have
seen and measures how a known program degrades with fewer donors.

## 2. Inputs and their provenance chains

All inputs already exist and are hashed in committed receipts; the runs verify every link
before computing anything and record the checks in their provenance files.

| Input | Path on the analysis host | Verified against |
|---|---|---|
| v1 discovery fit | `work/statistical/discovery_primary/de/model.rds` | object hash pinned by `receipts/discovery_rank_freeze.json` (`artifacts`), and its `provenance.json` (hash pinned by `receipts/discovery_primary_execution.json`, field `DE_provenance_sha256`), which records the same object hash under its original task path |
| v1 HECA fit | `work/statistical/heca_primary/de/model.rds` | object hash pinned by `receipts/heca_primary_execution.json` (`output_sha256`) |
| Phase discovery fit | `work/phase_positive_control/discovery_fit/de/model.rds` | object hash pinned by the published `reports/phase_positive_control_f1/de_provenance.json` |
| Phase HECA fit | `work/phase_positive_control/heca_fit/de/model.rds` | object hash pinned by the published `reports/phase_positive_control_f2/heca_de_provenance.json` |
| Positive-control count matrices and donor ledgers | `data/processed/phase_positive_control/{discovery,heca}/` | hashes in the F1 freeze and the evaluation receipt |
| Estimator configuration | `config/statistical_analysis.json`, `config/phase_positive_control.json` | frozen bytes (hashed in this protocol's freeze) |

Each saved fit is the complete limma object: per-gene unscaled standard deviations of the
coefficients, posterior residual variances after robust empirical-Bayes moderation, per-gene
prior and total degrees of freedom, the design matrix and the voom precision weights.

## 3. Part A — minimum detectable effect from the saved fits

### 3.1 Definition

For gene *g* and coefficient *c* (the disease contrast `ENDO_minus_non_ENDO` in the v1 fits;
`SEC_minus_PRO` in the phase fits), the v1 moderated statistic is
t = β̂ / (u · s), with u the unscaled standard deviation of the coefficient, s² the
posterior variance, and ν the per-gene total degrees of freedom (residual plus robust
prior). SE = u · s is the `moderated_SE` column already in every `de_all.tsv`.

Under a true effect δ (log2 units), t is treated as non-central t with ν degrees of freedom
and non-centrality δ / SE, holding SE and ν at their realized values. A two-sided test at
level α rejects when |t| > t<sub>1−α/2, ν</sub>. **MDE(α)** is the δ > 0 at which the
rejection probability equals 0.80, solved numerically (R `uniroot`, tolerance 10⁻⁸, with
`pt(…, ncp)`); the closed-form approximation (t<sub>1−α/2,ν</sub> + t<sub>0.8,ν</sub>) · SE
is recorded alongside as a check column. Genes that were not estimable in the original fit
get no MDE.

This is a **design sensitivity under the realized variance estimates**, not an observed
power: it uses the moderated variance of each gene, not its observed effect, and it is
conditional on the voom weights and prior that the original data produced.

### 3.2 Significance scenarios (fixed)

With *m* the tested family size of the fit (all expression-pass genes, as the v1 BH
adjustment used):

| Scenario | Per-gene α | Meaning |
|---|---|---|
| `alpha_0.05` | 0.05 | a single prespecified gene, no multiplicity |
| `bh_k1` | 0.05 × 1 / m | the p-value a gene needed to be the only BH q < 0.05 call |
| `bh_k10` | 0.05 × 10 / m | the boundary if exactly 10 genes were called |
| `bh_k100` | 0.05 × 100 / m | the boundary if exactly 100 genes were called |
| `bh_k1000` | 0.05 × 1000 / m | the boundary if exactly 1,000 genes were called |

The BH scenarios are the exact rejection boundaries of the Benjamini–Hochberg step-up rule
when exactly *k* genes are rejected; they describe what one gene needed given a number of
companions, not the power of the procedure as a whole.

### 3.3 Reported quantities (per fit)

- Per-gene table: gene, average expression, observed log2FC, SE, ν, |t|, p, q, 95% CI,
  MDE at each scenario, approximate MDE at `alpha_0.05`.
- MDE quantiles across estimable genes (5, 10, 25, 50, 75, 90, 95%) per scenario, in log2
  units and as fold changes (2^MDE).
- Share of genes with MDE ≤ 0.25, 0.5, 0.585 (1.5-fold), 1 and 2 log2 units per scenario.
- Median MDE and median SE by quintile of average expression (quantile type 7 bins of
  `AveExpr` among estimable genes).
- **What the observed intervals exclude:** the number of genes whose per-gene 95% CI lies
  entirely within (−Δ, Δ) for Δ ∈ {0.25, 0.5, 0.585, 1, 2}, and the same count for a
  family-wide Bonferroni interval (critical value t<sub>1−0.05/(2m), ν</sub>). The per-gene
  count is not simultaneous; the Bonferroni count is conservative. Both are reported.
- Observed descriptors: genes at q < 0.05, maximum |t|, and for the phase fits the
  quantiles of |log2FC| among q < 0.05 genes (the program's effect-size distribution).

### 3.4 Bridge between the phase program and the disease MDE (descriptive)

For each cohort, with *P* the phase-fit genes at q < 0.05 ("phase program"): the share of
*P* whose |log2FC| exceeds the disease fit's median MDE under `alpha_0.05` and under
`bh_k100`; and, for genes present in both families, the share of *P* whose phase |log2FC|
exceeds that same gene's disease-fit MDE under each scenario (gene-matched). These are
descriptions of two different contrasts on overlapping donors; they do not convert phase
effects into expectations about disease effects.

## 4. Part B — donor-subsampling sensitivity curve

### 4.1 Source and reference

Source matrices and ledgers are the frozen positive-control cohorts: discovery 13 secretory
vs 22 proliferative donors; HECA 20 vs 22. The reference for each cohort is its frozen
full-cohort phase fit: the set *S* of genes at BH q < 0.05 with their signs, and the top-20
genes by |moderated t| among estimable genes.

Before any subsample is fitted, the subsampling estimator is run on **all** donors with the
primary design and must reproduce the frozen full-cohort table (log2FC within 10⁻⁸, t within
10⁻⁶, log10 p within 10⁻⁶). Failure stops the run; nothing is reported.

### 4.2 Draws (fixed)

- Equal arms, *n* donors per phase arm, drawn **without replacement** within arm
  (subsampling, not bootstrap), stratified by phase.
- Levels: discovery *n* ∈ {5, 6, 8, 10, 11, 13}; HECA *n* ∈ {5, 6, 8, 10, 11, 13, 16, 20}.
  *n* = 11 matches the v1 total of 22 donors; the largest level uses the whole secretory arm.
- 100 draws per level; Mersenne-Twister, seed 42, draws generated in level order and
  written to a ledger before any fit. Draws are never redrawn.

### 4.3 Analysis of each draw (fixed)

Each draw is analyzed as a new small study with the v1 estimator unchanged: expression
filter recomputed on the subset (`filterByExpr`, phase as group, same parameters), TMM
recomputed, voom and the robust empirical-Bayes fit re-estimated, BH over the subset's own
family. Two designs are fitted for every draw:

- **primary**: intercept, `SEC_minus_PRO`, `ENDO_minus_non_ENDO`, centered age, batch
  treatment contrasts over the sorted batch levels present in the draw (capture batch for
  discovery, library lane for HECA), adenomyosis;
- **reduced**: the same without the batch contrasts.

Constant nuisance columns are dropped with a receipt. A design that is rank-deficient or
leaves fewer than 5 residual degrees of freedom makes that draw **invalid** for that design;
invalid draws are retained and counted, not replaced. The primary design is expected to be
infeasible at small *n* given the batch structure; the validity rate is itself a result.

### 4.4 Metrics per valid draw

| Metric | Definition |
|---|---|
| `n_sig` | genes at BH q < 0.05 in the draw (and up / down counts) |
| `recall_full` | share of *S* that is q < 0.05 in the draw with the same sign |
| `sign_disagreement` | among the draw's q < 0.05 genes, share whose sign differs from the full-cohort estimate |
| `top20_overlap` | size of the intersection between the draw's top-20 by \|t\| and the full-cohort top-20 |
| `typical_mde_alpha_0.05`, `typical_mde_bh_k100` | exact non-central-t MDE at the draw's median SE and median ν (k = 100 uses the draw's family size) |
| design receipt | residual df, rank, dropped nuisance columns, number of batch levels, family size |

### 4.5 Summaries

Per level and design: valid draws, invalid reasons, and the median with 2.5% and 97.5%
percentiles (type 7) of every metric. One figure per cohort with `n_sig`, `recall_full` and
the typical MDE against *n*. Draw-level rows are published in full.

## 5. Fixed interpretation rules

1. Part A sentences take the form "for gene *g*, a true effect of at least MDE(α) would have
   been detected with 80% probability under the realized precision"; summaries take the
   form "for half of the tested genes, MDE under scenario *k* is at most *x* log2 units". No
   sentence claims that effects below the MDE are absent, and none claims power for the BH
   procedure as a whole.
2. The CI-containment counts are the only statements of the form "effects larger than Δ are
   excluded for *N* genes"; the per-gene and family-wide versions are always reported
   together.
3. Part B describes how a large, known program degrades with fewer donors under this
   estimator. It is not a power calculation for disease effects, whose sizes are unknown,
   and it does not extrapolate beyond the levels drawn.
4. The bridge statistics are descriptive and are reported with their denominators; they do
   not justify any choice of donor number.
5. There is no pass/fail endpoint. Numbers are reported as computed, once, including
   unfavorable or uninformative ones (for example, if the primary design is never valid at
   a level).
6. No v1 or positive-control artifact is modified; their hashes are recorded before and
   after the run.

## 6. Order of operations

1. Freeze this protocol, the configuration and the code (`scripts/mde_freeze.py`), with the
   git head of the published checkout, before any MDE value is computed on real fits. A
   synthetic fixture (`scripts/mde_synthetic_fixture.R`, invented counts) may be used to
   test the code before the freeze; synthetic outputs are labelled and never published as
   results. If a run stops in its verification step before computing anything, the code
   may be corrected and the freeze re-issued; the new receipt names the superseded freeze
   commit and the reason, and the git history keeps both.
2. On the analysis host, in the v1 statistical runtime: verify the freeze and every input
   chain; run Part A on the four fits; run the identity check and Part B on both cohorts.
3. Fetch outputs, verify their hashes, publish compact copies under `reports/mde_analysis/`
   and write `receipts/mde_analysis_results.json` (`scripts/mde_finalize.py`); write the
   report `reports/mde_results.md` from the receipt. Remove untracked working copies from
   the host checkout and record what was kept.

## 7. What this analysis does and does not establish

It establishes, per gene and in aggregate, the effect sizes the v1 design was sensitive to
and the effect sizes its observed intervals exclude, and it measures how detection of a
real program depends on donor number under the frozen estimator. It does not establish
that disease effects are small in any specific gene, does not estimate disease effects,
does not rank targets, and does not choose a sample size for future work; it provides
the quantities such a choice would use.
