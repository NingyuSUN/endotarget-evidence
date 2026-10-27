# EndoTarget Evidence — research decision brief

**Decision: hold the candidate pool for further evidence. This public-data benchmark did not produce a statistically supported cross-cohort target shortlist or evidence that network propagation improves on expression ranking.** It demonstrates a traceable research decision, including negative results, rather than a validated therapeutic discovery. Follow-on analyses frozen before execution show that the same pipeline detects a known cycle-phase program at these sample sizes, bound the null to effects below roughly 2.2- to 3.4-fold for the median gene at the multiplicity faced, and find no detectable disease difference in the secretory-phase stromal program.

Prepared 14 October 2026; updated 27 October 2026 with four prespecified follow-on analyses. Scientific results were executed on SRV01; bounded container reproduction passed locally. A [two-page PDF edition](executive_brief_en.pdf) accompanies this source.

## Question and study design

Does a fixed network-based prioritization improve the transfer of stromal expression associations in endometriosis across public studies? The comparison uses author-adjudicated proliferative-phase eutopic endometrial stromal fibroblasts. Non-endometriosis surgical controls can have other gynecological conditions.

| Cohort | Analyzed people | Primary cells / nuclei | Genes tested | BH q < 0.05 |
|---|---:|---:|---:|---:|
| Discovery, GSE266265 | 12 cases + 10 controls | 150,788 cells | 17,733 | 0 |
| Replication, HECA E-MTAB-14039 | 15 cases + 7 controls | 52,744 nuclei | 17,643 | 0 |

Counts were aggregated per donor, then modeled separately by cohort with edgeR normalization and limma-voom, adjusting for age, source batch and adenomyosis. The sample sizes are **22 people per cohort**, not the cell totals. BH q values account for the complete tested gene family; zero significant genes means insufficient evidence under this design, not proof of no disease effect. See the [method card](../docs/METHOD_CARD.md).

## Frozen comparison and external results

All methods used the same 11,749 discovery-expressed network genes and 123,313 STRING functional-association edges. B0 ranks absolute moderated disease statistics; B1 ranks weighted network degree; G1 propagates discovery evidence through the frozen network. G1 is a restart random walk, with no trained graph neural network or parameter search. Lists, expression directions, parameters and matched reference sets were frozen before calculating HECA disease effects.

| Fixed top-20 method | Same direction / evaluable / selected | Same direction among evaluable | Same direction and external q < 0.05 |
|---|---:|---:|---:|
| B0: expression | 14 / 20 / 20 | 70.0% | 0 |
| B1: network degree | 3 / 13 / 20 | 23.1% | 0 |
| G1: propagation | 6 / 19 / 20 | 31.6% | 0 |

Seven B1 genes and one G1 gene failed external expression filtering; none was replaced. B1 misses the prespecified 16/20 evaluability gate. Direction agreement measures estimated expression signs, including small noisy effects; it is not target accuracy or therapeutic success.

G1 minus B0 was **−38.4 percentage points** in the fixed sample. Paired resampling of external donors completed 500 draws: 487 valid and 13 singular-design failures, retained without replacement. The conditional 95% interval was **−38.9 to +8.7 points**, including zero. It covers external sampling sensitivity given the fixed discovery lists/directions, not full discovery-selection or future-experiment uncertainty. No superiority claim follows.

All 22 discovery leave-one-donor-out analyses completed. Median retention of the original 20 genes was B0 14.5, B1 20 and G1 19.5; stable network lists did not yield better external sign agreement. Three viable clinical subset models also yielded no BH-significant genes; two other scenarios failed fixed donor minimums and were not fitted. The 30/50/100-cell thresholds select the same donors, so the primary fits were reused.

## What the null means: four follow-on analyses

Added 26–27 October 2026, each with rules, inputs and thresholds frozen and hashed before any new value was computed. They reuse the v1 estimator unchanged, refit no disease model and leave every v1 artifact byte-identical. None nominates a gene.

| Analysis | Result | What it bounds |
|---|---|---|
| Phase positive control: the v1 pipeline on proliferative versus secretory donors (35 and 42) | 7,367 and 3,986 genes at BH q < 0.05; the frozen expression top-20 transfers 20 / 20 same direction, all at external q < 0.05; G1 minus B0 −37 points, interval −53 to −21 | The pipeline detects and transfers a large donor-level program at these sample sizes; the disease null is not pipeline insensitivity. Propagation transfers worse than expression ranking even with a strong signal. |
| Minimum detectable effect from the saved v1 fits | Median gene: 1.5-fold at α 0.05, 2.2- to 2.3-fold at the BH boundary for 100 discoveries, 3.3- to 3.4-fold for one; at 11 + 11 donors a study recovers a median 45% (discovery) and 13% (HECA) of the phase program | Fourfold proliferative stromal disease differences are excluded for most genes, twofold for a fifth of genes family-wide; 1.5-fold effects stay open. Twenty-two donors recover even a several-thousand-gene program inconsistently. |
| Secretory-phase program contrast: each cohort scored on the other cohort's frozen phase program; phase × disease interaction | Signature transfers at AUC 1.00 and 0.99; θ +0.04 (−0.18 to +0.25) and −0.11 (−0.40 to +0.19) program units; detectable attenuation 0.31 and 0.42; no interaction gene at q < 0.05 | No detectable program-level disease difference in secretory stromal fibroblasts: discovery excludes an attenuation above 14% of the secretory shift, HECA admits anything short of a 71% loss. One nominal sensitivity (adenomyosis donors excluded, p 0.026) points to enhancement and is not reproduced. |
| Cohort extension audit: GSE213216 and GSE179640, metadata only | Natural-cycle eutopic donors 5 ENDO / 0 non-ENDO and 0 / 1; matched eutopic–lesion pairs 4 and 8; both studies lie inside the HECA cells atlas | No public single-cell cohort extends the eutopic donor-level contrasts; nothing was downloaded. The remaining eligible public donors are already in the analyses above. |

Secretory-phase genes are not endometriosis genes; the positive-control lists are cycle biology, not a shortlist. The [positive control](phase_positive_control_results.md), [detectable effect](mde_results.md), [secretory program](secretory_program_results.md) and [cohort extension](cohort_extension_feasibility.md) reports give denominators and receipts.

## Evidence integration and recommendation

All 40 genes in the frozen B0/G1 top-20 union remain visible. SKIDA1 has limited positional/QTL mapping in the prespecified 2023 genetics publication. The other 39 have unavailable mapping in that source, not negative genetic tests. Neither locus assignment nor QTL annotation establishes a causal gene or whether to activate or inhibit it.

The deterministic dossier rule selected only the inconclusive example PARP9: discovery log2 difference +0.466 versus external −0.142, with external 95% interval [−0.413, +0.130] and q = 0.970. Its opposite point-estimate sign does not establish a precise conflict. Recommendation: **hold**, pending better expression precision, cell-state applicability and independent mechanistic evidence. No antibody design, clinical utility or intervention direction is established.

## Audit trail and interpretation limits

[Six figure topics and linked panels](../docs/METHOD_CARD.md#figure-guide) cover donor selection, composition/QC, effect transfer, method sensitivity, the evidence matrix, and design resolution with the secretory program. The [full replication report](replication_v2_results.md), [clinical sensitivities](clinical_sensitivity_report.md), [composition audit](stromal_composition.md) and [genetic audit](genetic_evidence_report.md) provide denominators and provenance.

Discovery reuses author secondary Alevin estimated counts; the full original RDS was not recomputed. HECA has a raw-versus-H5 total discrepancy of 9,607/360,765,431 counts (0.00266%): 1,359 explained and 8,248 not assigned feature by feature. Inherited labels, scRNA/snRNA differences, heterogeneous subtype mixtures and residual clinical confounding limit interpretation. Recruiting sources are separated; genotype-confirmed individual non-overlap is unavailable.

AI assisted implementation and review; methods, evidence and limitations are documented in the audit trail. The local Docker/Nextflow check completed five fresh and five cached synthetic tasks. A [real-primary replay](../receipts/local_primary_container_replay.json) reproduced 34 scientific tables and eight JSON result files across five fresh tasks; discrete selections matched exactly and continuous values met absolute tolerance 1e-8. This is engineering reproduction of the same data. The [AWS guide](../docs/AWS_DEPLOYMENT.md) remains unexecuted; [local scope and receipts](../docs/LOCAL_EXECUTION.md) document the limits.
