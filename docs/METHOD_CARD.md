# EndoTarget Evidence — method card

**Use:** auditable observational benchmark of cross-study stromal expression prioritization in endometriosis. **Result:** no supported shortlist; no evidence that network propagation outperforms expression ranking. Research follow-up only.

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

Use edgeR `filterByExpr` (min.count 10, min.total.count 15, large.n 10, min.prop 0.7; disease-group size, not full-design leverage), recalculate library totals, TMM normalization, then limma-voom (fixed span 0.5), least-squares fit and robust empirical Bayes moderation. Fixed numeric designs contain disease, centered age, capture batch in discovery / literal library-group category in HECA, and adenomyosis. Discovery design rank/residual df = 12/10; HECA = 10/12. No covariate is removed to rescue an unfavorable primary result.

BH adjustment uses each complete expression-pass family: 17,733 discovery and 17,643 external genes, both with zero q < 0.05. Gene-level 95% intervals are conditional model intervals, not simultaneous or causal confidence statements. [Full decisions](STATISTICAL_IMPLEMENTATION_DECISIONS.md) and [machine parameters](../config/statistical_analysis.json) specify all options and failure handling.

## Split, baselines and evaluation

This is a discovery/external-cohort split, not a random cell split. Freeze donor IDs, discovery universe, network, lists, directions, parameters and 200 expression/degree-matched reference sets per method before external disease-effect computation. External count detectability does not alter discovery lists. The retrospective protocol is not preregistered or literature-blind.

| Method | Ranking on identical discovery universe |
|---|---|
| B0 | Descending absolute moderated discovery disease t statistic |
| B1 | Descending weighted network degree; topology control |
| G1 | Restart random walk seeded by normalized absolute discovery t; restart 0.5, L1 tolerance 1e−10, maximum 10,000 iterations; no training/tuning |

All directions come from the saved discovery disease coefficient; unsigned propagation is not activation/inhibition. [Ranking configuration](../config/ranking_algorithm.json) and [freeze](../receipts/discovery_rank_freeze.json) retain the definitions.

Primary endpoint: top-20 sign agreement, reported as agree/evaluable/selected without refill. Missing-outcome bounds retain denominator 20; fewer than 16 evaluable genes is coverage-limited. B0 = 14/20/20, B1 = 3/13/20, G1 = 6/19/20. No same-sign selected gene has external q < 0.05. Top-50 results, effect scatter and all missing reasons remain in the [evaluation report](../reports/replication_v2_results.md).

External bootstrap: same disease-stratified donor draw for every method; 500 draws, fixed discovery lists/directions and original external masks; renormalize/refit without refiltering. Keep 13 singular failures; 487 valid draws yield a conditional G1−B0 95% percentile interval [−0.38947, +0.08684]. Discovery leave-one-donor-out reruns filtering/model/rankings: 22/22 valid. Matched sets are descriptive, share genes/donors and do not supply a calibrated p value. [Resampling definitions](RESAMPLING.md) distinguish these uncertainty questions.

## Decision, limits and reproducibility

The frozen 40-gene pool yields one limited SKIDA1 mapping and 39 unavailable mappings; all intervention directions are unknown. Only inconclusive PARP9 is selected by the deterministic dossier rule: hold. Unknown evidence is not a negative test. No target efficacy, antibody-design, GNN-learning, causal or clinical claim follows.

Limitations include source labels and hormone/phase ascertainment, scRNA/snRNA differences, residual confounding, heterogeneous author-state proportions, discovery estimated-count uncertainty and the unresolved 8,248-count component of the HECA H5 discrepancy. Three viable clinical sensitivities were completed; two were not estimable. Their implementation followed primary results using previously recorded subset definitions, not preregistration.

Executed science: native SRV01 project environment; R 4.4.3, edgeR 4.4.0, limma 3.62.1, statmod 1.5.2; Python 3.10.22, NumPy 2.2.6, SciPy 1.15.2 in the [explicit lock](../requirements/statistical-runtime.explicit.txt). Use the preserved execution checkpoint for completed commands and [local continuation](LOCAL_EXECUTION.md) for the current location. The local Docker build and [manifest](../receipts/local_container_manifest.json) are verified. [Synthetic checks](../receipts/local_container_nextflow.json) completed five fresh and five cached tasks; [real-primary reproduction](../receipts/local_primary_container_replay.json) completed five fresh tasks and matched all 34 scientific tables plus eight JSON result metadata files. Discrete selections matched exactly; continuous values met absolute tolerance 1e-8 (largest difference 1.10e-13). The original 500 bootstrap and 22 donor-exclusion checkpoints were checked, not rerun. [AWS instructions](AWS_DEPLOYMENT.md) are unexecuted. Historical frozen status text is preserved rather than rewritten.

Freeze receipts record [pre-effect implementation](../receipts/statistical_pre_effect_freeze.json), [discovery](../receipts/discovery_rank_freeze.json) and [pre-external implementation](../receipts/replication_implementation_freeze.json). The [migration integrity receipt](../receipts/local_migration_integrity.json) checks preserved evidence, not new biological validation. Large inputs and `work/` checkpoints are ignored: a Git clone alone cannot reproduce the complete artifact closure. This is an AI-assisted project; software tests, AI reviews and hashes are not empirical validation.

## Figure guide

The roadmap's five scientific topics are covered by these linked panels; panels within a topic are not counted as additional biological experiments.

| Topic | Figure / panel | What it establishes |
|---|---|---|
| 1. Source-to-analysis flow | [Donor selection](../reports/figures/donor_flow_v1.svg) | Mutually exclusive losses and case/control donor denominators |
| 2. Donor composition/QC and balance | [Both cohorts, source states, primary counts, age and batch](../reports/figures/donor_composition_qc_v1.svg) | Descriptive donor heterogeneity; count depth is not purity |
| 3. Cross-study effects | [Effect scatter](../reports/figures/replication_v2_scatter.png) | Discovery/external associations, uncertainty and missingness |
| 4. Methods and sensitivity | [Transfer/coverage](../reports/figures/replication_v2_counts.png), [donor stability](../reports/figures/replication_v2_lodo.png), [paired bootstrap](../reports/figures/replication_v2_bootstrap.png), [matched sets](../reports/figures/matched_null_sensitivity_v1.svg) | Agreement, missingness, stability and distinct sampling/reference comparisons |
| 5. Evidence integration | [All 40 frozen genes](../reports/figures/evidence_matrix_v1.svg) | Expression, limited genetics, unknown conflicts/intervention; no supported shortlist |

The new SVGs only render existing tables. Their [renderer](../scripts/render_summary_figures_v1.py) runs with the Python standard library; [receipt](../receipts/summary_figures_v1.json) preserves input/script/output hashes and exact denominator checks. No disease model was rerun.
