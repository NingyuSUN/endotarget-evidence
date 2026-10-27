# Documentation

Start with the [method card](METHOD_CARD.md) for the study population, model,
method comparison, results and limitations. The completed analysis found no
supported target shortlist and no established advantage from network propagation.

## Results and interpretation

- [English research brief](../reports/executive_brief_en.pdf)
- [Results summary, English](../reports/results_summary_en.md)
- [Integrated biology report, Chinese](../reports/biology_report_zh.md)
- [External evaluation and donor sensitivity, Chinese/English generated report](../reports/replication_v2_results.md)
- [Clinical sensitivity analyses, Chinese generated report](../reports/clinical_sensitivity_report.md)
- [Stromal composition](../reports/stromal_composition.md) and [genetic evidence](../reports/genetic_evidence_report.md)
- [Phase positive-control results](../reports/phase_positive_control_results.md) — the same pipeline on a proliferative-versus-secretory contrast: 7,367 and 3,986 genes at BH q<0.05, the frozen top-20 transfers 20/20; the v1 disease null is not pipeline insensitivity

## Methods and source evidence

- [Input provenance](INPUT_PROVENANCE.md) and [final clinical eligibility](../reports/final_clinical_gate_report.md)
- [Statistical implementation decisions](STATISTICAL_IMPLEMENTATION_DECISIONS.md)
- [Ranking implementation](RANKING_IMPLEMENTATION.md) and [network audit](../reports/string_gene_reference_audit.md)
- [External evaluation definitions](REPLICATION_EVALUATION.md) and [resampling](RESAMPLING.md)
- [Genetic evidence method](GENETIC_EVIDENCE_METHOD.md)
- [Phase positive-control protocol](PHASE_POSITIVE_CONTROL_PROTOCOL.md) — frozen 2026-10-26 before any secretory-phase data; not executed
- [Phase positive-control feasibility](../reports/phase_positive_control_feasibility.md) — metadata only; both designs estimable
- [Phase positive-control HECA counts](../reports/phase_positive_control_heca_counts.md) — 42-donor HECA matrix prepared 2026-10-26; v1 columns reproduced exactly; no model fitted
- [Discovery secretory counts acquisition plan](PHASE_POSITIVE_CONTROL_DISCOVERY_COUNTS_PLAN.md) — measured resource constraints and procedure for the 41.7 GB GEO object; executed 2026-10-26
- [Phase positive-control discovery counts](../reports/phase_positive_control_discovery_counts.md) — 35-donor discovery matrix streamed from the GEO object; v1 columns reproduced to 4e-15; no model fitted
- [Phase positive-control discovery fit and F1 freeze](../reports/phase_positive_control_f1_discovery.md) — 7,367 genes at BH q<0.05; B0/B1/G1 lists frozen; marker panel 5/7
- [Phase positive-control evaluation receipt](../receipts/phase_positive_control_evaluation.json) and [published products](../reports/phase_positive_control_f2/) — HECA fit, P1/P2/S1/S2 against the frozen lists, 500-draw paired bootstrap
- [Minimum detectable effect protocol](MDE_ANALYSIS_PROTOCOL.md) — frozen 2026-10-27 before computation; per-gene MDE from the saved v1 fits and a donor-subsampling curve on the phase signal; no model refitted on disease
- [Minimum detectable effect results](../reports/mde_results.md) — median gene: 1.5-fold at α 0.05, 2.2–3.4-fold at the multiplicity faced; twofold effects excluded family-wide for a fifth of genes; a 22-donor study recovers 45% / 13% of the phase program
- [Secretory program disease contrast protocol](SECRETORY_PROGRAM_PROTOCOL.md) — frozen 2026-10-27 before any donor was scored; program-level test of an attenuated secretory stromal program in endometriosis on the positive-control matrices, signature defined in the other cohort
- [Cohort extension feasibility](../reports/cohort_extension_feasibility.md) — metadata only, 2026-10-27: GSE213216 and GSE179640 cannot extend the eutopic donor-level contrasts (natural-cycle eutopic donors 5 ENDO / 0 non-ENDO and 0 / 1); lesion-context resources only; nothing downloaded
- [Secretory program disease contrast results](../reports/secretory_program_results.md) — no detectable difference in either cohort: θ +0.04 [−0.18, +0.25] (discovery) and −0.11 [−0.40, +0.19] (HECA) program units; detectable attenuation 0.31 and 0.42 units; gene-level set tests agree; signature transfers at AUC 1.00 / 0.99

## Reproduction

- [Tested runtime, required artifacts and reproduction scope](LOCAL_EXECUTION.md)
- [Completed evidence index](LOCAL_DELIVERY_ACCEPTANCE.md)
- [Scientific package lock](../requirements/statistical-runtime.explicit.txt)
- [Actual container primary replay](../receipts/local_primary_container_replay.json)
- [AWS portability plan — unexecuted](AWS_DEPLOYMENT.md)

A Git clone does not contain the large count matrices, runtime assets or complete
`work/` checkpoints required for the scientific replay. The reproduction guide
describes the preserved environment rather than a fresh-clone installer.

## Immutable historical method evidence

The [research protocol](RESEARCH_PROTOCOL.md),
[pre-effect decisions](PRE_EFFECT_METHOD_DECISIONS.md),
[statistical implementation](STATISTICAL_IMPLEMENTATION_DECISIONS.md),
[external execution contract](REPLICATION_EXECUTION.md),
[evaluation definitions](REPLICATION_EVALUATION.md) and [resampling definitions](RESAMPLING.md)
retain the bytes recorded by the scientific freezes. Their original status wording
describes the time of each record, not the current completion state.

The [scientific-work audit](REMAINING_SCIENTIFIC_WORK.md) is likewise a
provenance record rather than an instruction to resume completed analyses. See
[publication history](PUBLICATION_HISTORY.md) for documentation changes and
the scope of historical verification receipts.
