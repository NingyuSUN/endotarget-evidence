# Secretory-phase program-level disease contrast — results

**Pair label under the frozen rules: `no_detectable_difference`.** In both cohorts the
secretory stromal program, defined in the other cohort, separates proliferative from secretory
donors almost perfectly, and its expression in secretory-phase endometriosis donors is not
measurably different from secretory-phase non-endometriosis donors. Discovery: θ = +0.04
program units, 95% interval −0.18 to +0.25. HECA: θ = −0.11, interval −0.40 to +0.19. The
designs could have detected an attenuation of 0.31 (discovery) and 0.42 (HECA) program units
with 80% power; the discovery interval excludes an attenuation larger than 14% of the observed
secretory shift, the HECA interval does not exclude anything short of a 71% loss. Gene-level
rotation tests on the same signatures agree in sign with the donor-level estimates and are not
significant. No gene reaches BH q < 0.05 for the interaction in either cohort.

Completed 2026-10-27 under the protocol frozen the same day before any donor was scored
([protocol](../docs/SECRETORY_PROGRAM_PROTOCOL.md); freeze
[`secretory_program_freeze.json`](../receipts/secretory_program_freeze.json); results receipt
[`secretory_program_results.json`](../receipts/secretory_program_results.json)). Executed once per
cohort on the analysis host (SRV01) in the v1 statistical runtime (R 4.4.3, edgeR 4.4.0, limma
3.62.1). Every pinned input hash matched; the recomputed expression filter, TMM factors and
design reproduced the frozen positive-control fit exactly (relative differences ≤ 5 × 10⁻¹⁵);
the 47 protected v1 and positive-control files were hash-identical before and after the run.

## Donors, signature and validity gate

| Cohort | Donors (secretory ENDO / non-ENDO) | Signature source | Signature genes in family (up / down) | Phase coefficient on the score (M0) | AUC secretory vs proliferative | Gate |
|---|---|---|---|---|---|---|
| Discovery, GSE266265 (scRNA) | 35 (7 / 6) | frozen HECA phase program | 3,580 of 3,986 (1,675 / 1,905) | +1.26 (SE 0.06; p 1 × 10⁻¹⁶) | 1.000 | pass |
| HECA, E-MTAB-14039 (snRNA) | 42 (16 / 4) | frozen discovery (F1) phase program | 6,570 of 7,367 (3,210 / 3,360) | +0.56 (SE 0.07; p 4 × 10⁻⁹) | 0.986 | pass |

Program projection score by group (mean; SD): discovery proliferative non-ENDO −0.01 (0.06),
proliferative ENDO +0.00 (0.14), secretory non-ENDO **1.28 (0.10)**, secretory ENDO **1.27
(0.16)**; HECA proliferative non-ENDO −0.05 (0.04), proliferative ENDO +0.02 (0.12), secretory
non-ENDO **0.56 (0.23)**, secretory ENDO **0.52 (0.24)**. A score of 1 is the signature
cohort's full phase shift; HECA nuclei project at about half of the discovery cell program,
with a wide spread across secretory donors (0.17 to 0.85) that follows sub-phase. The
equal-weight mean-signed-shift score correlates with the projection at 0.9998 and 0.9989.
Figures: [`discovery_scores.png`](secretory_program/discovery_scores.png),
[`heca_scores.png`](secretory_program/heca_scores.png).

## Primary estimand θ (M1: phase × disease interaction, all donors, positive-control design)

| Cohort | θ | 95% interval | SE | residual df | two-sided p | θ in secretory-score SD units | minimum detectable θ (80% power, α 0.05) | Label |
|---|---|---|---|---|---|---|---|---|
| Discovery | **+0.037** | −0.179 to +0.253 | 0.104 | 21 | 0.73 | +0.28 | 0.31 | `no_detectable_difference` |
| HECA | **−0.105** | −0.400 to +0.189 | 0.144 | 28 | 0.47 | −0.46 | 0.42 | `no_detectable_difference` |

Relative to each cohort's observed secretory shift (1.26 and 0.56 program units), the
intervals are compatible with an attenuation of at most 14% in discovery and 71% in HECA, and
with an enhancement of at most 20% and 34%. The proliferative-arm disease difference, the
descriptive baseline, is −0.04 (discovery) and +0.03 (HECA).

## Secondary estimands (secretory donors only)

| Cohort | M2: δ = ENDO − non-ENDO, age-adjusted | M2-sub: δ with sub-phase contrasts | Sub-phase effects in M2-sub |
|---|---|---|---|
| Discovery (13 donors) | −0.02 [−0.21, +0.16], p 0.81, df 10 | −0.06 [−0.24, +0.11], p 0.44, df 8 | mid vs early +0.20; late vs early +0.11 |
| HECA (20 donors) | +0.05 [−0.28, +0.37], p 0.76, df 17 | +0.10 [−0.14, +0.33], p 0.41, df 15 | Mid vs Early +0.39; unspecified vs Early +0.09 |

Sub-phase moves the score as declared in the protocol, mid-secretory donors scoring higher
than early ones in both cohorts; adjusting for it does not move δ away from zero.

## Sensitivities of θ (none changes a label)

| Sensitivity | Discovery θ [95%] | HECA θ [95%] |
|---|---|---|
| Reduced design, no batch contrasts | −0.00 [−0.19, +0.19] | −0.05 [−0.34, +0.24] |
| Mean-signed-shift score | +0.03 [−0.19, +0.26] | −0.12 [−0.43, +0.20] |
| Own-cohort program (circular) | +0.02 [−0.16, +0.20] | −0.15 [−0.59, +0.29] |
| Donors below 500 primary cells or nuclei excluded | identical (no donor excluded) | −0.18 [−0.59, +0.23]; 6 donors excluded, 3 non-ENDO secretory donors remain |
| Adenomyosis donors excluded | **+0.20 [+0.03, +0.37], p 0.026**; 7 donors excluded, 4 ENDO secretory donors remain | −0.10 [−0.44, +0.25]; 6 donors excluded |
| Leave-one-donor-out (35 and 42 refits) | θ from −0.03 to +0.10; largest change 0.07 (sample100); sign not constant | θ from −0.16 to −0.04; largest change 0.06 (FX0011); all negative, all p > 0.3 |

The discovery adenomyosis sensitivity is the only estimate whose interval excludes zero. It
points to a **larger** secretory shift in the remaining four endometriosis donors, not an
attenuation; it rests on a 4-versus-6 secretory comparison with 15 residual df, has a nominal
p of 0.026, is not reproduced in HECA under the same exclusion (−0.10), and by the frozen rules a
sensitivity never changes a label. The three excluded discovery ENDO secretory donors carry two
of the three lowest ENDO scores (1.08 and 1.12) and one above the mean (1.38). Whether
adenomyosis rather than endometriosis associates with a lower secretory program is a
hypothesis this analysis raises and cannot test.

## Gene-level counterpart (v1 estimator on the M1 design, coefficient `SEC_x_ENDO`)

| Cohort | Family | Genes at BH q < 0.05 | max \|t\| | fry on the signature (direction; p) | β-weighted projection of interaction log2 FC | Mean signed interaction log2 FC | Signature genes with interaction sign opposing β |
|---|---|---|---|---|---|---|---|
| Discovery | 17,341 | **0** | 5.08 | Up; p 0.60 | +0.036 | +0.040 | 46% |
| HECA | 17,243 | **0** | 5.30 | Down; p 0.38 | −0.081 | −0.107 | 63% |

The gene-level projections (+0.036 and −0.081) reproduce the donor-level θ (+0.037 and −0.105)
within the difference between voom-weighted and unweighted fitting, so the two routes are
measuring the same quantity. The published top-20 tables by |t| are descriptive; the smallest
interaction q value is 0.16 (HECA) and 0.60 (discovery).

## Readings

- **The hypothesis is not supported at program level in these cohorts, and not refuted
  either beyond a modest bound.** In discovery the secretory ENDO and non-ENDO donors sit on
  top of each other on a 3,580-gene program and the data exclude an attenuation larger than
  about one seventh of the secretory shift. In HECA the point estimate has the hypothesized
  sign in every leave-one-out refit, but four non-ENDO donors spanning the whole secretory
  range leave an interval that admits anything from a 71% loss to a 34% gain.
- **The signature transfers cleanly across platforms** (AUC 1.000 and 0.986 with the other
  cohort's gene set and effect sizes), so the null is not a failure of the program definition.
- **Sub-phase, not disease, is the visible source of score variation among secretory
  donors**, in both cohorts; the prespecified adjustment leaves δ at zero.
- **Design bound:** a program-level attenuation would have to exceed 0.31 (discovery) or 0.42
  (HECA) program units to be detected with 80% power at these donor numbers. Published
  whole-tissue findings of progesterone resistance describe gene-set-level dysregulation;
  these data can neither confirm nor exclude an effect of the size those studies imply within
  the stromal-fibroblast lineage, except that in discovery it would have to be smaller than
  14% of the secretory shift.

## What this does and does not establish

It establishes, under rules fixed in advance, that the secretory stromal program measured in
donor pseudobulk of the v1 lineage shows no detectable disease difference in two independent
cohorts, with the bounds above. It does not establish absence of progesterone resistance:
the HECA bound is wide, the discovery secretory arm has seven cases, sub-phase labels are the
authors' and may themselves absorb a delayed transition, and the lineage-level score mixes
stromal-state composition with within-state expression. No gene is nominated; the adenomyosis
observation is a sensitivity, not a finding. Nothing here revises a v1 or positive-control
number.

## Execution record and files

Freeze commit 7fd87ab (protocol, configuration and code), freeze receipt committed as
56c7e72; synthetic fixture test on the host before the freeze recovered a planted attenuation
(θ −0.43 [−0.50, −0.36] against a planted −0.5 after signature shrinkage) and gave no label to
the null fixture; host run 2026-10-27 06:27 UTC, 17 s per cohort; finalization verified every
provenance file and output hash on the published checkout. Published copies under
[`reports/secretory_program/`](secretory_program/): per-donor scores, signature gene lists,
leave-one-donor-out tables, interaction tables (full, compressed; top-20), identity checks,
summaries, provenance and figures for both cohorts. Host working copies were removed after the
run; the host keeps the outputs under its work tree.
