# EndoTarget Evidence

A public-data endometriosis target-evidence benchmark: donor-level stromal RNA, fixed expression, degree and network-propagation baselines, source-separated transcriptomic replication, and a separate genetics evidence layer. The original scientific analysis was executed on SRV01. Local continuation, container testing and bounded reproduction now use G:\projects\te in the project-owned EndoTargetRuntime WSL guest. The repository has not been published to GitHub.

**Result: the prespecified analysis did not produce a significantly supported cross-study target shortlist or evidence that network propagation improves on expression ranking.** Negative and not-estimable outcomes are retained.

## Start with the evidence

- [Two-page English research brief](reports/executive_brief_en.pdf) and [editable source](reports/executive_brief_en.md).
- [Method card and five figure topics](docs/METHOD_CARD.md).
- [Current local execution and reproduction](docs/LOCAL_EXECUTION.md) and [milestone acceptance map](docs/LOCAL_DELIVERY_ACCEPTANCE.md).

- [Integrated biology report, Chinese](reports/biology_report_zh.md).
- [Detailed bilingual replication report and figures](reports/replication_v2_results.md).
- [Clinical sensitivities and phase ascertainment](reports/clinical_sensitivity_report.md).
- [Donor subtype composition](reports/stromal_composition.md) and [genetic evidence](reports/genetic_evidence_report.md).
- Preserved SRV01 scientific checkpoint.

## Question, data and result

The fixed scope is author-adjudicated proliferative eutopic endometrial stromal fibroblasts.
Each donor contributes one aggregated expression profile. Separate models adjust for age,
technical batch and adenomyosis.

| Cohort | Case/control donors | Selected cells/nuclei | Tested genes | BH q<0.05 |
|---|---:|---:|---:|---:|
| GSE266265 discovery |12/10|150,788 cells|17,733|0|
| HECA E-MTAB-14039 replication |15/7|52,744 nuclei|17,643|0|

B0 ranks absolute moderated disease statistics; B1 ranks weighted network degree; G1 propagates
discovery evidence on a fixed STRING direct-evidence network. They share 11,749 genes/123,313 edges.
All lists, directions and 200 matched sets per method were frozen before external disease effects.

| Frozen top20 method | Agree/evaluable/selected | Same direction and external q<0.05 |
|---|---:|---:|
| Expression B0 |14/20/20|0|
| Degree B1 |3/13/20|0|
| Propagation G1 |6/19/20|0|

The paired G1-minus-B0 direction difference is −0.3842. Of 500 donor bootstrap draws, 487 were valid
and 13 failed the unchanged design checks; the conditional 95% percentile interval is
[−0.3895,+0.0868]. It includes zero. Direction agreement is not target accuracy or therapeutic utility.
All 22 discovery leave-one-donor-out folds completed. Network rankings were more stable, which
did not establish better external replication.

Three prespecified viable clinical sensitivities also had no BH-significant genes; two scenarios
failed donor minima and were not fitted. All 40 dossier-pool genes retain separate genetic annotations.
The deterministic rule selected one inconclusive PARP9 dossier and no supported target dossier.

## Reproducibility that can be inspected

Actual Nextflow 25.10.2 runs: three-stage discovery and two-stage replication, with pinned R/Python
runtime, input/output hashes, meaningful synthetic boundary tests, cached-resume checks, and
deterministic resampling. All 500 actual bootstrap checkpoints and 22 LODO output sets were checked.

- [Research contract](docs/RESEARCH_PROTOCOL.md) and [implementation decisions](docs/STATISTICAL_IMPLEMENTATION_DECISIONS.md).
- [Input provenance](docs/INPUT_PROVENANCE.md), [clinical gate](reports/final_clinical_gate_report.md), [marker audit](reports/lineage_marker_audit.md).
- [Pre-effect](receipts/statistical_pre_effect_freeze.json), [discovery](receipts/discovery_rank_freeze.json), and [replication implementation](receipts/replication_implementation_freeze.json) freezes.
- [HECA execution](receipts/heca_primary_execution.json) and [real resampling receipt](receipts/real_resampling_execution.json).
- [Replication instructions](docs/REPLICATION_EXECUTION.md), [runtime lock](requirements/statistical-runtime.explicit.txt), [AWS deployment guide](docs/AWS_DEPLOYMENT.md).

Large data/work artifacts were copied to local G:\projects\te; SRV01 originals are preserved. They are excluded from Git, so a clone alone does not restore the complete frozen checkpoint. Historical scientific documents retain their original status wording. The prior homepages are byte archives in docs/archive; current local status is described above.

## Limits and release status

Author phase/cell labels, source-qualified hormone eligibility, possible residual clinical and cell
mixture differences, secondary discovery aggregation and the unresolved small HECA raw/H5 count
difference are explicit in the reports. Recruitment sources are separated; genotype-based individual
non-overlap and cell purity are not established. Genetics annotations do not establish causal or
intervention direction.

The local Docker build, registry-manifest check and actual Nextflow container checks have passed. Synthetic tests completed five fresh and five cached tasks. A fresh real-primary replay completed all five tasks: all 34 scientific TSVs and eight JSON result metadata files agreed with the original science; discrete selections matched exactly and the maximum continuous difference was 1.10e-13 (declared absolute tolerance 1e-8). The 117 scientific packages remain locked. The additive local recipe pins 37 OS packages supporting task metrics and real Git provenance. AWS is a documented, unexecuted portability plan. No therapeutic validation, wet experiment, cloud deployment or successful target-discovery claim is made.

Actual evidence: [container/cache checks](receipts/local_container_nextflow.json), [real-primary reproduction](receipts/local_primary_container_replay.json), [final delivery integrity verdict](receipts/local_delivery_verification.json).
