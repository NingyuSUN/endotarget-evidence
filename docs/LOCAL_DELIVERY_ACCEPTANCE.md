# Completed evidence index

This index connects the reported scientific conclusions to their source,
analysis and reproduction evidence. The study did not produce a supported target
shortlist or establish that network propagation improves external expression
replication. A positive target finding is not required to report the benchmark.

| Evidence area | Delivered artifacts | What the evidence establishes |
|---|---|---|
| Population and inputs | [Clinical eligibility](../reports/final_clinical_gate_report.md), [input provenance](INPUT_PROVENANCE.md) | Audited donor/lineage/phase selection within the documented source-label and count limitations |
| Donor summaries | [Composition audit](../reports/stromal_composition.md), [donor flow](../reports/figures/donor_flow_v1.svg), [composition and balance](../reports/figures/donor_composition_qc_v1.svg) | 22 donors per cohort; 150,788 discovery cells and 52,744 external nuclei; donors are the statistical units |
| Fixed method comparison | [Method card](METHOD_CARD.md), [pre-effect freeze](../receipts/statistical_pre_effect_freeze.json), [discovery freeze](../receipts/discovery_rank_freeze.json) | Expression, degree and unsigned propagation share one universe; no trained GNN or external-outcome tuning |
| External evaluation and sensitivity | [Replication report](../reports/replication_v2_results.md), [clinical sensitivities](../reports/clinical_sensitivity_report.md), [matched sets](../reports/figures/matched_null_sensitivity_v1.svg) | Complete lists and missingness; 487/500 valid external bootstrap draws with 13 failures retained; 22/22 discovery donor-exclusion fits; no established propagation advantage |
| Genetics and research decision | [Genetic evidence](../reports/genetic_evidence_report.md), [40-gene matrix](../reports/figures/evidence_matrix_v1.svg) | One inconclusive PARP9 dossier, no supported shortlist and unknown intervention directions; unavailable mappings are not negative genetic tests |
| Scientific interpretation | [Biology report](../reports/biology_report_zh.md), [method card and figures](METHOD_CARD.md), [English brief](../reports/executive_brief_en.pdf) | Integrated conclusions, denominators, uncertainty and limitations |
| Software reproduction | [Runtime and prerequisites](LOCAL_EXECUTION.md), [synthetic fresh/resume checks](../receipts/local_container_nextflow.json), [real-primary replay](../receipts/local_primary_container_replay.json), [continuous software checks](../.github/workflows/ci.yml) | Five fresh and five cached synthetic tasks; five fresh real-primary tasks; 34 scientific TSVs and eight JSON result metadata files reproduced; synthetic unit tests and a link check on every push |
| Follow-on: pipeline sensitivity (October 2026) | [Positive-control protocol](PHASE_POSITIVE_CONTROL_PROTOCOL.md), [results](../reports/phase_positive_control_results.md), [evaluation receipt](../receipts/phase_positive_control_evaluation.json) | Same pipeline detects 7,367 and 3,986 phase genes and transfers the frozen top-20 20/20; the v1 disease null is not pipeline insensitivity |
| Follow-on: design resolution | [MDE protocol](MDE_ANALYSIS_PROTOCOL.md), [results](../reports/mde_results.md), [results receipt](../receipts/mde_analysis_results.json) | Median gene detectable at 1.5-fold per test and 2.2–3.4-fold at the multiplicity faced; 22-donor studies recover a known program inconsistently |
| Follow-on: secretory program | [Protocol](SECRETORY_PROGRAM_PROTOCOL.md), [results](../reports/secretory_program_results.md), [results receipt](../receipts/secretory_program_results.json) | No detectable program-level disease difference in secretory stromal fibroblasts in either cohort; detectable attenuation 0.31 / 0.42 program units |
| Follow-on: public cohort extension | [Feasibility report](../reports/cohort_extension_feasibility.md), [audit receipt](../receipts/cohort_extension_metadata_audit.json) | No eligible natural-cycle eutopic disease stratum in GSE213216 or GSE179640; nothing downloaded |
| Follow-on summary figure | [Figure](../reports/figures/follow_on_summary_v1.svg), [receipt](../receipts/summary_figure_follow_on_v1.json) | Rendered from the published summaries above; no model rerun |

## Interpretation and execution boundaries

The three primary scientific freezes cover 95 distinct artifacts. The primary
replay retained exact discrete selections and used absolute tolerance 1e-8 for
continuous scientific values; the maximum observed difference was 1.10e-13.
It did not rerun raw import, 500 bootstrap fits, 22 donor-exclusion fits or
clinical sensitivities. Their original checkpoints remain evidence of the
original executions.

Engineering reproduction of the same data is distinct from cross-study
transcriptomic replication and from therapeutic validation. No target efficacy,
causal mechanism, intervention direction, wet-lab experiment or cloud execution
is established by these software checks.

A Git clone excludes the large inputs and `work/` products required to restore
the complete checkpoint. See the [reproduction guide](LOCAL_EXECUTION.md) for
artifact, runtime and path prerequisites. The [AWS guide](AWS_DEPLOYMENT.md)
remains an unexecuted plan. The [publication history](PUBLICATION_HISTORY.md)
defines the scope of historical verification receipts and documentation hashes.
