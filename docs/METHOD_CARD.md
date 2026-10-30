# EndoTarget Evidence — method card

**Use:** auditable observational benchmark of cross-study stromal expression prioritization in endometriosis. **Result:** no supported shortlist; no evidence that network propagation outperforms expression ranking. Research follow-up only. **Follow-on analyses (October 2026):** the same pipeline detects and transfers a known cycle-phase program at these sample sizes; the v1 null is bounded to effects below about 2.2- to 3.4-fold for the median gene at the multiplicity faced; the secretory-phase stromal program shows no detectable disease difference; no further public eutopic single-cell cohort can extend the design ([details](#follow-on-analyses-what-the-null-resolves)).

## Inputs and analysis population

| Component | Actual input and retained scope |
|---|---|
| Discovery | GSE266265 author secondary Alevin aggregation; 22 selected donors (12 ENDO / 10 non-ENDO), 150,788 cells; fractional estimated counts retained without rounding; full source RDS not recomputed |
| External evaluation | HECA newly generated nuclei, E-MTAB-14039; 22 donors (15 / 7), 52,744 nuclei; original integer counts, audited donor/library selection |
| Population | Author-adjudicated proliferative eutopic endometrium, primary stromal-fibroblast boundary; ≥50 primary cells/nuclei per donor; author QC and labels inherited |
| Network | STRING 12.5 human functional associations; direct experiments/database channels, maximum channel score ≥0.7; no transferred, text-mining or coexpression channels; 11,749 induced non-isolated genes, 123,313 undirected edges |
| Genetics | One prespecified 2023 GWAS publication and its public mapping/QTL supplement; separate annotation of all 40 frozen pool genes; no de novo GWAS, colocalization or Mendelian randomization |

[Clinical ledgers](../manifests/final_discovery_clinical_selected.tsv), [HECA ledger](../manifests/final_heca_clinical_selected.tsv), [input provenance](INPUT_PROVENANCE.md), [clinical eligibility](../reports/final_clinical_gate_report.md) and [network audit](../reports/string_gene_reference_audit.md) define the inputs. Controls are surgical non-ENDO donors, not necessarily healthy. Different recruiting sources support source separation; individual non-overlap is not genotype-confirmed.

## Estimand and donor model

One donor pseudobulk is one observation. The disease coefficient is adjusted ENDO minus non-ENDO normalized log2 expression. It does not separate changes within states from changes in subtype proportions or identify a causal disease effect.

The analysis used edgeR `filterByExpr` (min.count 10, min.total.count 15, large.n 10, min.prop 0.7; disease-group size, not full-design leverage), recalculated library totals, TMM normalization, then limma-voom (fixed span 0.5), least-squares fit and robust empirical Bayes moderation. Fixed numeric designs contain disease, centered age, capture batch in discovery / literal library-group category in HECA, and adenomyosis. Discovery design rank/residual df = 12/10; HECA = 10/12. No covariate is removed to rescue an unfavorable primary result.

BH adjustment uses each complete expression-pass family: 17,733 discovery and 17,643 external genes, both with zero q < 0.05. Gene-level 95% intervals are conditional model intervals, not simultaneous or causal confidence statements. [Full decisions](STATISTICAL_IMPLEMENTATION_DECISIONS.md) and [machine parameters](../config/statistical_analysis.json) specify all options and failure handling.

## Split, baselines and evaluation

This is a discovery/external-cohort split, not a random cell split. Donor IDs, discovery universe, network, lists, directions, parameters and 200 expression/degree-matched reference sets per method were frozen before external disease-effect computation. External count detectability does not alter discovery lists. The retrospective protocol is not preregistered or literature-blind.

| Method | Ranking on identical discovery universe |
|---|---|
| B0 | Descending absolute moderated discovery disease t statistic |
| B1 | Descending weighted network degree; topology control |
| G1 | Restart random walk seeded by normalized absolute discovery t; restart 0.5, L1 tolerance 1e−10, maximum 10,000 iterations; no training/tuning |

All directions come from the saved discovery disease coefficient; unsigned propagation is not activation/inhibition. [Ranking configuration](../config/ranking_algorithm.json) and [freeze](../receipts/discovery_rank_freeze.json) retain the definitions.

Primary endpoint: top-20 sign agreement, reported as agree/evaluable/selected without refill. Missing-outcome bounds retain denominator 20; fewer than 16 evaluable genes is coverage-limited. B0 = 14/20/20, B1 = 3/13/20, G1 = 6/19/20. No same-sign selected gene has external q < 0.05. Top-50 results, effect scatter and all missing reasons remain in the [evaluation report](../reports/replication_v2_results.md).

External bootstrap: same disease-stratified donor draw for every method; 500 draws, fixed discovery lists/directions and original external masks; renormalize/refit without refiltering. The 13 singular failures were retained; 487 valid draws yield a conditional G1−B0 95% percentile interval [−0.38947, +0.08684]. Discovery leave-one-donor-out reruns filtering/model/rankings: 22/22 valid. Matched sets are descriptive, share genes/donors and do not supply a calibrated p value. [Resampling definitions](RESAMPLING.md) distinguish these uncertainty questions.

## Decision, limits and reproducibility

The frozen 40-gene pool yields one limited SKIDA1 mapping and 39 unavailable mappings; all intervention directions are unknown. Only inconclusive PARP9 is selected by the deterministic dossier rule: hold. Unknown evidence is not a negative test. No target efficacy, antibody-design, GNN-learning, causal or clinical claim follows.

Limitations include source labels and hormone/phase ascertainment, scRNA/snRNA differences, residual confounding, heterogeneous author-state proportions, discovery estimated-count uncertainty and the unresolved 8,248-count component of the HECA H5 discrepancy. Three viable clinical sensitivities were completed; two were not estimable. Their implementation followed primary results using previously recorded subset definitions, not preregistration.

The original scientific analysis ran in the native SRV01 environment: R 4.4.3,
edgeR 4.4.0, limma 3.62.1, statmod 1.5.2; Python 3.10.22, NumPy 2.2.6 and
SciPy 1.15.2. The [explicit lock](../requirements/statistical-runtime.explicit.txt)
records the scientific packages. The completed Docker build and
[manifest check](../receipts/local_container_manifest.json) fixed the reproduction
environment. [Synthetic checks](../receipts/local_container_nextflow.json)
completed five fresh and five cached tasks. A
[real-primary replay](../receipts/local_primary_container_replay.json) completed
five fresh tasks and reproduced all 34 scientific TSVs and eight JSON result
metadata files. Discrete selections matched exactly; continuous scientific
values met absolute tolerance 1e-8 with zero relative tolerance (largest
difference 1.10e-13). The original 500 bootstrap and 22 donor-exclusion
checkpoints were checked, not rerun.

Freeze receipts record [pre-effect implementation](../receipts/statistical_pre_effect_freeze.json),
[discovery](../receipts/discovery_rank_freeze.json) and
[pre-external implementation](../receipts/replication_implementation_freeze.json).
The [reproduction guide](LOCAL_EXECUTION.md) defines the tested runtime,
required artifacts and replay limits. Large inputs and `work/` checkpoints
are excluded from Git: a clone alone does not restore the complete artifact
closure. [AWS instructions](AWS_DEPLOYMENT.md) remain unexecuted. Software
tests, code review and hashes are not biological or clinical validation.

## Follow-on analyses: what the null resolves

Four prespecified analyses were added after the v1 result, on 26–27 October 2026. Each froze its rules, inputs and thresholds in a hashed receipt before any new value was computed, reused the v1 estimator unchanged, refitted no disease model, and left the 47 protected v1 and positive-control files byte-identical. Each ran once on SRV01 in the v1 statistical runtime. None nominates a gene.

| Analysis | Design | Result | What it establishes |
|---|---|---|---|
| [Phase positive control](PHASE_POSITIVE_CONTROL_PROTOCOL.md) | v1 estimator, network, ranking and evaluation code unchanged; contrast secretory minus proliferative with disease, centered age, batch and adenomyosis adjusted; discovery 13 + 22 donors (225,207 cells), HECA 20 + 22 (94,150 nuclei); B0/B1/G1 lists frozen in discovery before the HECA fit | 7,367 and 3,986 genes at BH q < 0.05 (P1 pass); frozen B0 top-20 in HECA 20/20 same direction and 20/20 at q < 0.05 (P2 pass); marker panel 5/7 and 6/7 (explanatory S1 not met); G1 − B0 = −0.37, conditional 95% interval [−0.53, −0.21], 499/500 valid draws | The pipeline is sensitive to a large donor-level program and transfers it across platform and study, so the v1 disease null is not pipeline insensitivity. With a strong signal to propagate, G1 on this graph transfers worse than B0. [Results](../reports/phase_positive_control_results.md) |
| [Minimum detectable effect](MDE_ANALYSIS_PROTOCOL.md) | Exact non-central t solution per gene from the saved v1 moderated SE and total df, at α 0.05 and at the BH boundary for 1,000, 100, 10 and 1 discoveries; interval-exclusion counts; 100 donor-subsampling draws per level on the phase cohorts with the v1 estimator re-estimated per draw | Median gene MDE 1.5-fold at α 0.05, 2.2- / 2.3-fold at the boundary for 100 discoveries, 3.3- / 3.4-fold for one (discovery / HECA); per-gene 95% intervals within ±1 log2 for 72% / 78% of genes, family-wide Bonferroni intervals for 21% / 20%; at 11 + 11 donors the primary design recovers a median 45% / 13% of the full phase program (2.5th percentiles 18% / 0.1%) and is never estimable at 5 or 6 per arm | The v1 data exclude fourfold proliferative stromal disease differences for most genes under any reading and twofold differences for a fifth of genes family-wide; they say little about 1.5-fold differences family-wide. A 22-donor study of this design recovers even a several-thousand-gene program inconsistently. No gene's effect is thereby shown to be small. [Results](../reports/mde_results.md) |
| [Secretory program disease contrast](SECRETORY_PROGRAM_PROTOCOL.md) | Positive-control matrices; each cohort scored on the other cohort's frozen phase program (β-weighted projection on voom log-CPM; 1 = the signature cohort's full secretory shift); validity gate AUC ≥ 0.9; primary θ = secretory × ENDO interaction on all donors in the positive-control design; secondary δ in secretory donors with age and sub-phase; v1 gene-level interaction fit and limma `fry`; labels fixed before scoring | Gate passed at AUC 1.000 / 0.986 on 3,580 / 6,570 signature genes; θ = +0.037 [−0.179, +0.253], p 0.73 (discovery, 7 ENDO / 6 non-ENDO secretory donors) and −0.105 [−0.400, +0.189], p 0.47 (HECA, 16 / 4); 80%-power MDE 0.31 and 0.42 program units; δ −0.02 and +0.05; 0 interaction genes at q < 0.05; `fry` p 0.60 / 0.38; pair label `no_detectable_difference` | No detectable program-level disease difference in secretory stromal fibroblasts in two cohorts: discovery excludes an attenuation above 14% of the secretory shift, HECA admits anything short of a 71% loss. The one nominal sensitivity (discovery, adenomyosis donors excluded: +0.20 [+0.03, +0.37], p 0.026, 4 versus 6 donors) points to enhancement, is not reproduced in HECA and by rule changes no label. [Results](../reports/secretory_program_results.md) |
| Cohort extension audit | Metadata only: GEO SOFT files, file listings and the HECA authors' donor curation for GSE213216 and GSE179640; v1 eligibility rules (eutopic, single-cell, natural cycle, phase-resolved, age present); floor 5 donors per disease group within a phase, 5 pairs for a paired lesion design | Natural-cycle eutopic donors 3 ENDO / 0 non-ENDO proliferative and 2 / 0 secretory (GSE213216), 0 / 1 (GSE179640, a hormonally treated cohort by design); matched eutopic–lesion single-cell pairs 4 and 8; both studies lie inside the HECA cells atlas | Neither study can extend the eutopic donor-level contrasts or support a donor-level meta-analysis; nothing was downloaded. The remaining natural-cycle, phase-resolved public eutopic single-cell donors are already in the analyses above. [Report](../reports/cohort_extension_feasibility.md) |

Receipts: positive control [F0](../receipts/phase_positive_control_freeze.json), [F1](../receipts/phase_positive_control_f1_discovery_freeze.json) and [evaluation](../receipts/phase_positive_control_evaluation.json); MDE [freeze](../receipts/mde_analysis_freeze.json) and [results](../receipts/mde_analysis_results.json); secretory program [freeze](../receipts/secretory_program_freeze.json) and [results](../receipts/secretory_program_results.json); [extension audit](../receipts/cohort_extension_metadata_audit.json). Secretory-phase genes are not endometriosis genes: the positive-control lists are cycle biology, not a shortlist. The MDE values are conditional on the realized variance estimates of the original fits; the program-level bounds rest on 13 and 20 secretory donors, 6 and 4 of them non-ENDO; effects below the stated bounds are not shown to be absent.

## Figure guide

These panels cover six scientific topics. Panels within a topic are not additional biological experiments.

| Topic | Figure / panel | What it establishes |
|---|---|---|
| 1. Source-to-analysis flow | [Donor selection](../reports/figures/donor_flow_v1.svg) | Mutually exclusive losses and case/control donor denominators |
| 2. Donor composition/QC and balance | [Both cohorts, source states, primary counts, age and batch](../reports/figures/donor_composition_qc_v1.svg) | Descriptive donor heterogeneity; count depth is not purity |
| 3. Cross-study effects | [Effect scatter](../reports/figures/replication_v2_scatter.png) | Discovery/external associations, uncertainty and missingness |
| 4. Methods and sensitivity | [Transfer/coverage](../reports/figures/replication_v2_counts.png), [donor stability](../reports/figures/replication_v2_lodo.png), [paired bootstrap](../reports/figures/replication_v2_bootstrap.png), [matched sets](../reports/figures/matched_null_sensitivity_v1.svg) | Agreement, missingness, stability and distinct sampling/reference comparisons |
| 5. Evidence integration | [All 40 frozen genes](../reports/figures/evidence_matrix_v1.svg) | Expression, limited genetics, unknown conflicts/intervention; no supported shortlist |
| 6. Design resolution and the secretory program | [Follow-on summary](../reports/figures/follow_on_summary_v1.svg); [v1 detectable effects, discovery](../reports/mde_analysis/v1_discovery_mde_ecdf.png) and [HECA](../reports/mde_analysis/v1_heca_mde_ecdf.png); [subsampling curves, discovery](../reports/mde_analysis/curve_discovery.png) and [HECA](../reports/mde_analysis/curve_heca.png); [program scores, discovery](../reports/secretory_program/discovery_scores.png) and [HECA](../reports/secretory_program/heca_scores.png) | Per-gene detectable effects of the v1 fits by scenario; recovery of the phase program by donor number and design; donor-level program scores by phase and disease |

The SVGs render existing tables. Their [renderer](../scripts/render_summary_figures_v1.py) runs with the Python standard library; [receipt](../receipts/summary_figures_v1.json) preserves input/script/output hashes and exact denominator checks. No disease model was rerun. The topic-6 PNGs were drawn by the follow-on analysis scripts from their own outputs; the results receipts linked above record their hashes. The follow-on summary SVG is rendered from the published summaries by [its script](../scripts/render_summary_figure_follow_on_v1.py), whose [receipt](../receipts/summary_figure_follow_on_v1.json) records the input, script and output hashes and every number drawn.
