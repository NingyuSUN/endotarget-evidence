# Secretory-phase program-level disease contrast — prespecified protocol

Version 1.0, 2026-10-27. Status: frozen before any donor was scored and before any
disease-within-secretory quantity was computed (receipt
`receipts/secretory_program_freeze.json`). Machine-readable parameters:
[`config/secretory_program.json`](../config/secretory_program.json). Every definition,
model, sensitivity, label and interpretation sentence below is fixed here; nothing is chosen
after seeing a result.

## 1. Purpose and hypothesis

The v1 disease analysis ([method card](METHOD_CARD.md)) tested endometriosis against
non-endometriosis in **proliferative-phase** stromal fibroblasts and found no gene at
BH q < 0.05 in either cohort. The [phase positive control](PHASE_POSITIVE_CONTROL_PROTOCOL.md)
showed that the same estimator detects a large stromal program (proliferative to secretory,
7,367 and 3,986 genes) at the same donor numbers, and the
[minimum detectable effect analysis](MDE_ANALYSIS_PROTOCOL.md) bounded what the v1 design
could have seen (median gene: 1.5-fold at α 0.05; 2.2–3.4-fold at the multiplicity faced).

The most cited molecular hypothesis for eutopic endometrium in endometriosis is
**progesterone resistance**: an attenuated progesterone response in the secretory phase,
with impaired decidualization of the stroma. Transcriptome studies of whole eutopic
endometrium reported a dysregulated proliferative-to-secretory transition and an attenuated
progesterone-regulated program in endometriosis (Burney et al. 2007), and the hypothesis is
the subject of a large review literature (Patel et al. 2017); decidualization biology is
summarized by Gellersen and Brosens (2014). The v1 design could not address it: it had
no secretory donors. The positive control created, for other reasons, secretory-phase
donor matrices in both cohorts with both disease groups represented, and a frozen
definition of the secretory program in each cohort.

This protocol tests the hypothesis once, at the **program level**, with a single primary
estimand per cohort:

> **θ (theta)** — the phase-by-disease interaction on a donor-level program score: the
> difference between the secretory shift of endometriosis donors and the secretory shift of
> non-endometriosis donors, each shift measured relative to the cohort's proliferative
> donors. Under progesterone resistance θ < 0.

What this is not: it is not a gene-level screen, it nominates no target, it does not revise
v1 or positive-control outputs, and a positive result would be an association in eutopic
endometrium of surgically characterized donors, not a mechanism. The secretory arms are
small (13 and 20 donors; 4 non-endometriosis donors in HECA), so the main product may be a
bound rather than a detection; both are reported under the same rules.

## 2. Timing disclosure

- The secretory count matrices and donor ledgers exist locally since the positive control;
  their hashes are in the F1 freeze and the evaluation receipt. No new data are obtained.
- The positive-control fits adjusted the phase contrast for disease additively. Their
  disease coefficient was never extracted, tabulated or reported: the published tables carry
  the `SEC_minus_PRO` columns only, and the saved objects were read by the MDE analysis for
  that coefficient only. This protocol's estimand (an interaction, or a secretory-only
  contrast) does not exist in those fits.
- No program score, no secretory-only disease contrast and no interaction has been computed
  in either cohort. The donor composition in section 4 is metadata the positive-control
  documents already publish; it was examined to write the sensitivities, which is permitted
  design information, not outcome information.
- The expected direction is public knowledge; this is a literature-informed confirmatory
  protocol, frozen before effects, not a blind one.

## 3. Inputs and provenance chains

Every input already exists and is pinned by a committed receipt. The analysis verifies each
hash before computing anything and records the checks in its provenance file.

| Input | Path on the analysis host | Pinned by |
|---|---|---|
| Positive-control count matrices and donor ledgers | `data/processed/phase_positive_control/{discovery,heca}/` | F1 freeze (discovery), evaluation receipt (HECA) |
| Frozen discovery phase fit: full table, expression filter, normalization, design | `work/phase_positive_control/discovery_fit/de/` | `receipts/phase_positive_control_f1_discovery_freeze.json`, `artifacts` |
| Frozen HECA phase fit: full table, expression filter | `work/phase_positive_control/heca_fit/de/` | `receipts/phase_positive_control_evaluation.json`, `artifacts` |
| Frozen HECA normalization and design | `work/phase_positive_control/heca_fit/de/` | the same receipt through their published byte copies `reports/phase_positive_control_f2/heca_de_{normalization,design}.tsv` |
| Estimator and phase configuration | `config/statistical_analysis.json`, `config/phase_positive_control.json` | frozen bytes, hashed in this protocol's freeze |

The frozen phase tables are the complete `de_all` tables of the positive control (log2 fold
change, moderated t, BH q and estimability flag per gene). The published compressed copies
of both tables are hashed in the freeze as well.

## 4. Population and declared imbalances (metadata)

The donors are exactly the positive-control donors; nobody is added, removed or
reclassified. Proliferative arms are the v1 donors.

| Cohort | Proliferative ENDO / non-ENDO | Secretory ENDO / non-ENDO | Secretory sub-phase by disease | Secretory adenomyosis (ENDO) | Secretory age range ENDO / non-ENDO |
|---|---|---|---|---|---|
| Discovery (GSE266265, scRNA) | 12 / 10 | **7 / 6** | ENDO: early 2, mid 5, late 0; non-ENDO: early 2, mid 2, late 2 | 3 of 7 | 29–41 / 33–43 |
| HECA (E-MTAB-14039, snRNA) | 15 / 7 | **16 / 4** | ENDO: early 3, mid 7, unspecified 6; non-ENDO: mid 2, unspecified 2 | 3 of 16 | 22–40 / 18–32 |

Declared in advance, with the consequence each has for the analysis:

- **Sub-phase composition differs by disease** in discovery (late-secretory donors occur
  only among non-ENDO). A monotone program would score late donors higher, which could
  produce θ < 0 without any disease biology. The sub-phase-adjusted secretory model
  (section 7) is prespecified for this reason. Sub-phase labels are the authors' labels
  and are not re-derived.
- **Age differs by disease** among HECA secretory donors (non-ENDO younger). Age is a
  covariate in every model.
- **Four non-ENDO secretory donors in HECA.** Every HECA disease estimate rests on them.
  The leave-one-donor-out analysis is prespecified for this reason.
- **Depth varies widely in HECA** (113 to 7,873 primary nuclei per donor; four donors below
  500). The depth-minimum sensitivity is prespecified for this reason.
- **Batch** cannot be adjusted within the secretory arms (seven capture batches among 13
  donors; seven lanes among 20). The primary model uses all donors and the positive-control
  batch contrasts; the secretory-only models do not adjust for batch, and the secretory
  batch-by-disease table is published.
- **HECA stromal labels** include decidualized states, so the program partly measures
  stromal-state composition. This was declared for the positive control and holds here.

## 5. Signature: the other cohort's frozen program

For each test cohort the signature is defined in the **other** cohort, so that gene
selection is independent of the donors being tested:

- discovery is scored on the frozen HECA phase program; HECA is scored on the frozen
  discovery (F1) phase program;
- signature genes are those with `fit_estimable` TRUE and BH q < 0.05 for `SEC_minus_PRO`
  in the signature cohort's frozen table, restricted to genes in the test cohort's
  expression-pass family;
- each signature gene carries the signature cohort's log2 fold change β<sub>g</sub> as its
  weight and sign(β<sub>g</sub>) as its direction.

The test cohort's own program (its own q < 0.05 genes) is used only in a sensitivity
labeled circular, because it was selected on the same donors; it never contributes to a
label. The positive control's 20-of-20 transfer of the discovery top-20 into HECA is the
reason to expect the cross-cohort signature to carry the program; the validity gate in
section 7 checks it formally.

## 6. Expression scale and donor scores

Before any scoring, the analysis recomputes the positive-control expression filter (edgeR
`filterByExpr` with the phase group), the TMM factors and the primary design from the count
matrix and the ledger, and requires them to equal the frozen tables exactly (filter), within
10⁻¹⁰ relative (normalization) and within 10⁻⁸ (design). A failed identity check stops the
run.

Expression is the voom log-CPM, log2((count + 0.5) / (filtered library size × TMM factor + 1)
× 10⁶), on the expression-pass genes. For gene g, m<sub>g</sub> is the mean log-CPM over the
cohort's proliferative donors. Two scores per donor d:

- **Projection (primary):** Σ<sub>g</sub> β<sub>g</sub> (x<sub>gd</sub> − m<sub>g</sub>) /
  Σ<sub>g</sub> β<sub>g</sub>², over signature genes in the family. It is the least-squares
  regression of the donor's centered profile on the signature direction: 0 is the
  proliferative mean, 1 is the signature cohort's full phase shift. Platform differences
  (scRNA versus snRNA) scale the slope but do not affect a within-cohort disease contrast.
- **Mean signed shift (secondary):** mean<sub>g</sub> sign(β<sub>g</sub>)
  (x<sub>gd</sub> − m<sub>g</sub>), in log2 units, equal-weighted.

Scores are computed once per cohort and signature and are not re-estimated inside any
model or resample.

## 7. Models and estimands

All models are ordinary least squares on donor scores with numeric treatment-contrast
designs, two-sided t tests, 95% intervals, rank tolerance 10⁻⁷, and residual df ≥ 5; a
rank-deficient design, too few residual df, or a positive-control design without a batch
containing both phases yields `NOT_ESTIMABLE` for that model, with no relaxation. Age is
centered on the donors entering each model. A constant nuisance column in a subset is
dropped with a receipt.

| Model | Donors | Design | Coefficient | Role |
|---|---|---|---|---|
| **M0** | all | the positive-control primary design: intercept, `SEC_minus_PRO`, `ENDO_minus_non_ENDO`, centered age, batch contrasts, adenomyosis | `SEC_minus_PRO` | validity gate |
| **M1** | all | M0 plus `SEC_x_ENDO` = `SEC_minus_PRO` × `ENDO_minus_non_ENDO` | `SEC_x_ENDO` = **θ** | **primary** |
| **M2** | secretory | intercept, `ENDO_minus_non_ENDO`, centered age | `ENDO_minus_non_ENDO` = δ | secondary |
| **M2-sub** | secretory | M2 plus sub-phase treatment contrasts (discovery: early, mid, late; HECA: Early, Mid, unspecified; reference = first level present) | `ENDO_minus_non_ENDO` | secondary, sub-phase sensitivity |

**Validity gate (M0).** The projection score must separate phases: the `SEC_minus_PRO`
coefficient on the score is positive with p < 0.001, and the area under the ROC curve of the
score for secretory versus proliferative donors is ≥ 0.9. If the gate fails, the cohort is
`NOT_INTERPRETABLE` and no disease label is issued for it.

**Primary estimand θ** from M1 on the projection score. Reported with SE, t, residual df,
two-sided p, 95% interval, θ in units of the secretory-donor score SD, and the **minimum
detectable θ** at 80% power and α 0.05 from the realized SE and df (same noncentral-t
solver as the MDE analysis). The proliferative-arm disease difference (the
`ENDO_minus_non_ENDO` coefficient of M1) is reported as a descriptive baseline.

**Secondary estimand δ** from M2 and M2-sub, reported the same way. δ answers the direct
question (are secretory ENDO donors lower?) without the proliferative anchor and without
batch; M2-sub removes sub-phase composition at the cost of df (8 in discovery).

**Sensitivities of θ** (all M1, all reported, none changes a label): (i) the reduced
design without batch contrasts; (ii) the mean-signed-shift score; (iii) the own-cohort
program (circular); (iv) donors below 500 primary cells or nuclei excluded; (v) adenomyosis
donors excluded; (vi) leave-one-donor-out: M1 refitted without each donor in turn, with the
range of θ, the largest absolute change, the most influential donor and whether every
refit keeps the sign of the primary estimate.

## 8. Gene-level counterpart

The v1 estimator (filter, TMM, voom span 0.5, least-squares fit, robust empirical Bayes,
BH over the expression-pass family) is run once per cohort on the M1 design, and the
`SEC_x_ENDO` coefficient is reported for every gene: the number of genes at BH q < 0.05
(expected near zero at these sizes), the full table, and the top-20 by |t| labeled
descriptive. The program-level statement at gene level is a rotation gene-set test
(limma `fry`, directional, gene weights = sign(β<sub>g</sub>), contrast `SEC_x_ENDO`):
"Down" means the program is attenuated in endometriosis. Two descriptive aggregates
accompany it: the β-weighted projection of the interaction log2 fold changes onto the
signature (the gene-level analogue of θ) and the mean signed interaction log2 fold change.
This is a new disease fit by design; it is declared as such and modifies nothing frozen.

## 9. Multiplicity and evidence labels, fixed in advance

The confirmatory family is {θ<sub>discovery</sub>, θ<sub>HECA</sub>}. Labels are applied
mechanically by the finalization step from the frozen rules; secondary and sensitivity
results never change a label.

Per cohort:

| Label | Rule |
|---|---|
| `attenuation_detected` | θ < 0 and two-sided p < 0.025 (Bonferroni over the two cohorts) |
| `enhancement_detected` | θ > 0 and two-sided p < 0.025 |
| `no_detectable_difference` | otherwise; the 95% interval and the minimum detectable θ are the reported bound |
| `NOT_INTERPRETABLE` | validity gate failed |
| `NOT_ESTIMABLE` | M1 not estimable |

Pair:

| Label | Rule |
|---|---|
| `replicated_attenuation` | θ < 0 with p < 0.05 in both cohorts |
| `replicated_enhancement` | θ > 0 with p < 0.05 in both cohorts |
| `single_cohort_signal` | exactly one cohort at p < 0.025, either direction, not replicated |
| `no_detectable_difference` | neither cohort at p < 0.025 and not replicated |
| `NOT_EVALUABLE` | a cohort is `NOT_INTERPRETABLE` or `NOT_ESTIMABLE`; the other cohort's label stands alone |

## 10. Interpretation rules, fixed in advance

- `replicated_attenuation`: the secretory stromal program is expressed less fully in
  eutopic endometrium of endometriosis donors in two independent cohorts and platforms,
  consistent with progesterone resistance at the program level. This is an association;
  sub-phase, age and depth sensitivities say whether composition could produce it. No gene
  is nominated; the gene-level table is descriptive.
- `single_cohort_signal`: reported as such, with the other cohort's interval; not
  replicated, and not escalated by sensitivities.
- `no_detectable_difference`: an attenuation larger than the interval bound is unlikely in
  these cohorts; the minimum detectable θ states what this study could have seen. With 13
  and 20 secretory donors the bound is expected to be wide, and the result is a design
  statement, not evidence of absence.
- `replicated_enhancement` or an `enhancement_detected` cohort: reported as the opposite
  direction, with the same caveats; not reinterpreted as support.
- A gate failure is a transfer or measurement problem, documented before anything else is
  concluded about disease; the cohort's disease label is withheld, not inferred.
- The gene-level fry test and aggregates are secondary; they agree or disagree with θ and
  are reported either way.
- No number from this analysis enters a v1 table, figure, dossier or brief; v1 and
  positive-control artifacts remain byte-identical and their hashes are rechecked after the
  run.

## 11. Order of operations and freezes

1. **Freeze** (this document, the configuration, the code, the receipts that pin the
   inputs, the count matrices and ledgers, the published phase-table copies):
   `scripts/secretory_program_freeze.py` cross-checks every pinned hash against the
   committed receipts and refuses to run if results exist. The code is tested on a synthetic
   fixture (`scripts/secretory_program_synthetic_fixture.R`: invented cohorts with a planted
   program and a planted attenuation in one cohort) before the freeze; no real donor is
   scored before the freeze.
2. **Host run**, once per cohort, on the analysis host in the v1 statistical runtime:
   freeze verification, pinned-input verification, identity checks, scores, models,
   sensitivities, gene-level fit, provenance with output hashes. Protected artifacts are
   hashed before and after the run.
3. **Finalization** on the published checkout: `scripts/secretory_program_finalize.py`
   verifies the freeze and every provenance file, applies the pair label, publishes copies
   under `reports/secretory_program/` and writes `receipts/secretory_program_results.json`.
4. The results report (`reports/secretory_program_results.md`) presents the numbers under
   the labels above.

If the host run stops before any score or model is produced because of a defect in the
code, the defect is corrected, the freeze is re-issued with a `supersedes` field naming the
earlier freeze commit and the reason, and both freezes remain in the history. No rule,
threshold, label or sensitivity may change in a re-issue.

## 12. What this does and does not establish

It establishes whether, at the program level and in these two cohorts, the secretory
stromal program differs between endometriosis and non-endometriosis eutopic endometrium,
with a fixed decision rule and a bound when nothing is detected. It does not establish a
mechanism, a within-state effect (the lineage includes decidualized states), a gene, or a
target; it cannot separate sub-phase composition from disease beyond the prespecified
adjustment; and the HECA estimate rests on four non-endometriosis donors.

## 13. Sources

According to PubMed: Burney et al. 2007, *Endocrinology*
([DOI 10.1210/en.2006-1692](https://doi.org/10.1210/en.2006-1692); dysregulated
proliferative-to-secretory transition and attenuated progesterone response in eutopic
endometrium of endometriosis); Patel et al. 2017, *Acta Obstet Gynecol Scand*
([DOI 10.1111/aogs.13156](https://doi.org/10.1111/aogs.13156); progesterone resistance
review); Gellersen and Brosens 2014, *Endocr Rev*
([DOI 10.1210/er.2014-1045](https://doi.org/10.1210/er.2014-1045); decidualization);
Marečková et al. 2024, *Nat Genet*
([DOI 10.1038/s41588-024-01873-w](https://doi.org/10.1038/s41588-024-01873-w); HECA);
Duempelmann et al. 2026 ([DOI 10.1038/s41467-026-73020-4](https://doi.org/10.1038/s41467-026-73020-4);
discovery cohort).
