# Local delivery acceptance — 14 October 2026

This is the current acceptance map for the original roadmap.
The user moved remaining execution to the local computer; the original SRV01
scientific runs and their frozen evidence remain preserved. Historical design/status
files are not rewritten to suggest they knew later results.

## Required evidence

| Milestone | Delivered evidence | Interpretation |
|---|---|---|
| M0 / ETE-01–04: population, inputs and feasibility | [Final clinical gate](../reports/final_clinical_gate_report.md), [input provenance](INPUT_PROVENANCE.md), execution checkpoint | Donor/lineage/phase eligibility and resource gates passed within the documented source-label and count limitations |
| M1 / ETE-05–06: counts and donor summaries | [Composition audit](../reports/stromal_composition.md), [donor flow](../reports/figures/donor_flow_v1.svg), [composition and balance](../reports/figures/donor_composition_qc_v1.svg) | 22 donors in each cohort; 150,788 discovery cells and 52,744 external nuclei; cells are not independent patient replicates |
| M2 / ETE-07–08: fixed baselines and discovery freeze | [Methods](METHOD_CARD.md), [pre-effect freeze](../receipts/statistical_pre_effect_freeze.json), [discovery freeze](../receipts/discovery_rank_freeze.json) | Expression, degree and unsigned network propagation on the same universe; no trained GNN or external-outcome tuning |
| M3 / ETE-09: external effects and sensitivity | [Replication report](../reports/replication_v2_results.md), [clinical sensitivity](../reports/clinical_sensitivity_report.md), [matched sets](../reports/figures/matched_null_sensitivity_v1.svg) | Complete lists, missingness, 487/500 valid bootstrap draws with 13 failures retained, and 22/22 donor-exclusion fits; no established network advantage |
| M4 / ETE-10: separate genetics and research decision | [Genetics report](../reports/genetic_evidence_report.md), [40-gene evidence matrix](../reports/figures/evidence_matrix_v1.svg), [English brief](../reports/executive_brief_en.pdf) | One inconclusive PARP9 dossier and no supported shortlist; unavailable genetics is not a negative test; intervention direction unknown |
| M4–5 / ETE-11: readable report and method card | [Biologist report](../reports/biology_report_zh.md), [method card and five figure topics](METHOD_CARD.md), [two-page brief source](../reports/executive_brief_en.md) | Reports existing results and their limits, with no new scientific fit for figure rendering |
| M5 / ETE-12: clean local reproduction | [Local runtime and reproduction](LOCAL_EXECUTION.md) | Same immutable image: [5 fresh + 5 cached synthetic tasks](../receipts/local_container_nextflow.json) and [5 fresh real-primary tasks](../receipts/local_primary_container_replay.json) passed; [final artifact verdict](../receipts/local_delivery_verification.json) verifies the delivered evidence |

Completion does not require a positive target or a winning graph method. It requires
an auditable answer under the fixed protocol. Engineering reproduction of these data
is distinct from external biological replication and from therapeutic validation.

## Resource and release boundaries

Only the discovery three-task chain and external two-task chain are designated for
the clean real-data replay. Raw import, 500 bootstrap fits, 22 donor-exclusion fits
and clinical sensitivities are not repeated for this delivery. Their original
checkpoints are preserved and checked against their recorded hashes.

The local runtime shares WSL's 8 GiB memory / 2 GiB swap limit. Nextflow runs one
container task at a time with a 2 GiB memory cap and a 512 MiB driver heap. A cap is
not a measured peak. The full migrated scientific/work/log tree contains about
3.48 GB; this excludes added container images, package downloads, build caches and
the WSL virtual disk. The real replay took 96.46 seconds and the synthetic fresh/resume check 87.93 seconds.
The [resource snapshot](../receipts/local_container_resources.json), taken between those checks,
recorded 7.70 GB used in the new guest ext4 filesystem and 442.8 GB free on G. Docker image
and cache accounting share layers and cannot be summed as independent physical storage.
No peak resource or AWS price/performance measurement is inferred.

The [AWS guide](AWS_DEPLOYMENT.md) remains an unexecuted plan. Its original
container-status wording is historical; current local evidence is recorded in
[LOCAL_EXECUTION.md](LOCAL_EXECUTION.md). No cloud resources or external image
publication are part of the local execution described above. GitHub distribution is a separate
publication step; see the current repository homepages. A Git clone
excludes large ignored inputs and work products; full checkpoint reproduction also
requires the hashed artifact bundle described in the execution documentation.

All 95 distinct scientific freeze artifacts remain unchanged. The previous two
homepages are byte-preserved in [the archive](archive/INDEX.md), while the current
homepages describe this local delivery. The final verifier handles only these two
explicit historical-path exceptions and independently hashes the new documents.
