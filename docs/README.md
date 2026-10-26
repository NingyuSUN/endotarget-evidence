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
