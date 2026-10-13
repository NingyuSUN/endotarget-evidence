# EndoTarget Evidence

A public-data endometriosis target-evidence project using donor-level stromal RNA, a fixed association network, and planned independent genetic annotation. All execution currently runs on SRV01 under /home/analyst/te. The Git repository is local; it has not been published to GitHub.

**Current result (2026-10-13): real discovery analysis completed and frozen; external replication is pending.** Among 22 discovery donors (12 ENDO/10 non-ENDO), 17,733 genes passed the fixed expression filter. None met BH q<0.05 (minimum q=0.5101). This is insufficient single-gene evidence under this model, not proof of no expression difference.

## Question and comparison

Can discovery stromal expression priorities transfer to a source-separated public cohort, and does network propagation improve on simple baselines? The fixed primary scope is author-adjudicated proliferative-phase eutopic endometrial stromal fibroblasts. One donor is one statistical observation.

- B0: absolute moderated disease t statistic.
- B1: weighted network degree.
- G1: fixed-restart propagation of discovery expression evidence.

The models adjust for age, the adopted technical batch grouping, and adenomyosis. All methods compare the same discovery-expressed induced graph: 11,749 genes and 123,313 edges. Their expression directions come from the discovery coefficient. STRING associations do not establish physical binding, causal regulation or intervention direction.

## Actual input and execution evidence

| Cohort | Selected donors | Selected cells/nuclei | Stable input genes | Disease analysis |
|---|---:|---:|---:|---|
| GSE266265 discovery, audited secondary donor aggregation | 12 ENDO / 10 non-ENDO | 150,788 cells | 26,109 | Completed |
| HECA new E-MTAB-14039 nuclei, public raw matrices | 15 ENDO / 7 non-ENDO | 52,744 nuclei | 26,125 | Pending |

Discovery retains original fractional Alevin estimates. HECA uses 22 public integer-count archives. The source-independent recruitment evidence, inherited labels, phase/hormone qualifications and mapping losses remain documented; no genotype-based participant non-overlap or cell-purity claim is made.

The actual three-stage Nextflow25.10.2 workflow completed donor fitting, ranking and matched-set sampling in the isolated native environment. Full/top20/top50 rankings, all-U directions and 200 matched 20-gene sets per method are frozen before external effects. B0/G1 top20 overlap is 0/20; B1/G1 overlap is 9/20, making hub preference an explicit question for external evaluation.

- [Readable discovery results](reports/discovery_primary_results.md) and [actual execution receipt](receipts/discovery_primary_execution.json).
- [Pre-effect freeze](receipts/statistical_pre_effect_freeze.json) and [complete discovery freeze](receipts/discovery_rank_freeze.json).
- [Source provenance](docs/INPUT_PROVENANCE.md), [clinical audit](reports/final_clinical_gate_report.md), [fixed marker audit](reports/lineage_marker_audit.md).
- [Statistical decisions](docs/STATISTICAL_IMPLEMENTATION_DECISIONS.md), [ranking implementation](docs/RANKING_IMPLEMENTATION.md), [matched sets](docs/NULL_SAMPLING.md).
- [Runtime and locks](docs/STATISTICAL_RUNTIME.md), [Nextflow execution](docs/NEXTFLOW_EXECUTION.md), current continuation state.

## Validation and remaining work

24 synthetic R checks and 96 Python tests passed. The complete synthetic Nextflow chain ran all three stages, reused all three on resume with byte-identical outputs, and rejected changed or incomplete freezes. This establishes software behavior, not scientific replication or therapeutic utility.

External replication, paired donor bootstrap, discovery leave-one-donor-out stability, sensitivities, independent genetic evidence, dossiers and the complete biology report remain unfinished. The discovery result will not be rescued by changing the selected phase, covariates or significance threshold.

Docker client29.0.1 is installed but the project account cannot access the daemon. The Docker profile is an unexecuted design; no container build/run or AWS deployment is claimed. Docker access details. AWS documentation remains a planned portability deliverable.

The external raw counts exceed saved H5 totals by 9,607 of 360,765,431 UMIs (0.00266%); the complete cause is not established. Counts were preserved. Controls may have other gynecological conditions. Source limitations remain part of any final interpretation.

To inspect this populated SRV01 checkpoint without rerunning analyses, use the manifest verification commands in CURRENT_STATE. Fresh acquisition, container execution and cloud portability are not yet a completed end-to-end release.
