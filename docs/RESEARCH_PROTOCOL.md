# Research protocol — design v0.1

Date: 2026-10-13. Status: **PROPOSED, not preregistered or executed**.

Pre-effect implementation clarifications are adopted in [PRE_EFFECT_METHOD_DECISIONS.md](PRE_EFFECT_METHOD_DECISIONS.md): use author-adjudicated proliferative phase consistently across cohorts, retain source-qualified hormone evidence, declare pathology/batch handling and reversible compact-count ingestion. These supersede tentative clinical scenarios; the complete M0 and ranking freezes remain incomplete. No disease effects have been computed.
 This protocol was drafted after reading public study descriptions. It is a retrospective analysis with a planned freeze before computing external disease effects, not a blinded validation claim.

## 1. Decision and scope

Produce an evidence-based shortlist for follow-up research in endometriosis. The scientific question is whether stromal expression associations transfer across studies and remain plausible after considering independent evidence and alternative explanations.

Primary lineage: endometrial stromal fibroblasts, using a documented author-label/marker mapping. Exclude perivascular/mural, endothelial, epithelial and immune populations from this primary definition. Do not equate all mesenchymal cells with endometrial stroma. Do not choose the best-looking subtype after disease testing.

Primary tissue: eutopic endometrium in reproductive-age donors. Exclude lesions, organoids, postmenopausal tissue and donors on exogenous hormones from the primary comparison. Non-endometriosis surgical controls may have other conditions; preserve control indications rather than calling all controls healthy.

V1 includes one primary lineage and one disease contrast. Macrophages, perturbation-response prediction, spatial analysis, cell communication, full GWAS fine-mapping/MR, docking and foundation-model training are outside v1. scRNA-seq plus genetic annotations is integrative analysis, not a claim of matched multiomic profiling.

## 2. Stage 0: inputs must pass before analysis

Candidate discovery: GSE266265. Candidate replication: newly generated HECA nuclei from E-MTAB-14039. Resource existence is verified; file contents and eligibility are not.

Required audit outputs:

- A donor/sample/cell mapping with disease, source cohort, tissue, assay, cycle phase and its ascertainment, hormone exposure, age where available, batch and author cell labels.
- An independence ledger using study, recruiting site, recruitment period and public donor/sample provenance. Unresolved overlap means the second resource cannot be labeled independent. Different accession numbers alone do not prove independence.
- Nonnegative, finite raw/author-filtered RNA counts, stable gene IDs, and an explicit counts-layer verification. Integrated/SCT/log values are not counts. Noninteger estimated counts require a documented method decision before analysis, not blind integer rounding.
- A cross-tabulation of disease by phase, batch, age availability, tissue and control indication. Missing critical fields remain missing; no label inference from expression to manufacture eligibility.
- Import size, peak-memory feasibility, disk use, public access and licensing audit. Do not launch large downloads before checking storage and conversion requirements.

### Fixed primary phase and eligibility

Use a **proliferative-phase-only** comparison: at least 8 eligible case and 8 control donors in discovery, and at least 5 per group in replication, each with >=50 eligible stromal cells/nuclei. These are pragmatic minimums, not power guarantees. Public metadata supports this as a feasible candidate contrast; eligibility after cell QC is not yet established.

Exclude ambiguous/periovulatory, secretory, menstrual and exogenous-hormone categories from the primary comparison. Do not broaden the phase after inspecting results or to recover failed eligibility. If the fixed contrast fails, stop v1 inferential analysis and issue a feasibility report or write a revised protocol before proceeding.

Full-rank design, residual degrees of freedom >=5 and disease/batch overlap are also required in each cohort. Apply all donor minima and design checks to the FINAL model matrix after cell QC and age/other covariate missingness exclusions, not only to the initial metadata. A perfectly confounded disease/batch contrast fails; batch correction cannot reconstruct the missing comparison. Unknown cycle or hormone status fails primary eligibility.

Use >=50 cells/nuclei per donor-lineage as the primary measurement threshold and 30/100 as sensitivity thresholds. Recheck donor minima and design identifiability for every sensitivity fit; failed checks yield 'not estimable', not a relaxed model. Freeze all eligibility rules and exact cohort IDs before generating disease rankings. Publish all exclusions; threshold-related exclusions may be nonrandom.

## 3. Counts, annotations and the statistical unit

Reuse documented author QC as the starting point, keeping inherited filtering visible. Add only predeclared technical filters, not disease-dependent trimming. Audit broad lineage labels with a fixed marker panel; do not reuse published disease-associated gene lists as an annotation target or ranking seed.

Aggregate unintegrated RNA counts by donor and primary lineage. Technical libraries from the same biopsy may be summed; repeated biological samples are not independent donors and require a predeclared selection/repeated-measure rule. Multiple technical libraries do not increase n.

Use R edgeR library-size normalization followed by limma-voom donor-level modeling. Primary disease coefficient is ENDO minus non-ENDO. Use disease plus age and batch when available, identifiable and not redundant. Specify the exact design from metadata before effect calculation; exclude missing age from an age-adjusted analysis and report the unadjusted sensitivity, rather than silently dropping that covariate to improve results. If age is unavailable across an entire study, report residual age confounding explicitly. Record original cycle day and phase ascertainment; a broad proliferative label does not remove all hormonal variation.

Filter genes with a fixed count-based expression rule within each study; freeze the exact rule/software version at implementation. Discovery determines candidate eligibility. The external study may assess label-blind detectability, never feed its disease effects into discovery selection.

Use original counts for inference. Cohort integration or embeddings, if used for visualization, must not supply the differential-expression input. Model cohorts separately.

Report donor-level cell-state proportions and marker patterns to distinguish lineage-composition changes from within-state expression changes. Predefined stromal subtype sensitivity is descriptive unless each subtype independently passes sample thresholds; no extra confirmatory target rankings.

Method rationale: [Squair et al., 2021](https://www.nature.com/articles/s41467-021-25960-2) emphasizes accounting for biological replicate variation. This supports donor-level analysis, not automatic validity of every pseudobulk result.

## 4. Freeze and non-circular comparison

Before external disease effects are computed, commit: protocol/config, eligible donor IDs, gene universe, annotation crosswalk, design matrices, software/network snapshots, all three ranking tables, top-20/top-50 lists, and evaluation rules. Produce hashes and a timestamp.

Permitted before the freeze: external provenance, metadata and label-blind feature-coverage checks. Disallowed: external DE-guided thresholds, retuning graph parameters, reselecting cell types, or replacing unsuccessful candidates. Public papers have already been read; this is an auditable freeze, not ignorance of previous literature.

## 5. Three graph-ranking baselines

Comparative universe U: start with discovery-expressed genes with unambiguous identifiers, construct their induced subgraph using the frozen permitted human-network edges, then remove nodes with zero degree IN THAT induced subgraph and freeze U. An edge only to a gene outside this set does not prevent isolation. Report how many discovery genes are excluded by this choice. Keep a separate all-expressed-gene DE table for biological interpretation, but compare methods on identical U.

Use a frozen human STRING evidence-channel download. Retain experimental/database-supported edges, excluding text-mining-only and expression-coexpression-only support. At implementation, document whether the selected network is a functional association or physical network; never describe association edges as causal regulation. Pin the original file, channel definitions, score transformation, version and hash. If channel separation is unavailable, the graph comparison fails its gate rather than silently using a text-enriched combined score.

| Method | Frozen definition | Purpose |
|---|---|---|
| B0: expression only | Descending absolute moderated disease t statistic on U | Strong, simple transcriptomic ranking |
| B1: degree only | Descending weighted degree on the same U/network | Detect preference for heavily connected genes |
| G1: network propagation | Restart random walk seeded by normalized absolute discovery t statistics | Test whether neighborhood evidence improves transfer |

G1: symmetric nonnegative adjacency; keep edges with the maximum permitted channel score >=0.7 after scaling to [0,1], edge weight equal to that maximum. Exclude self loops and collapse duplicate undirected edges by maximum. This is an operational rule, not a calibrated interaction probability. Column-normalize adjacency into W. Seed s is abs(t) normalized to sum 1. Iterate p_next = 0.5*s + 0.5*W*p until L1 change <1e-10 or 10,000 iterations; nonconvergence is a failure. If all seed values vanish or U has fewer than 50 genes, no meaningful ranking benchmark is produced. Break ranking ties by stable gene ID.

No parameter search in v1. No external disease information, GWAS labels, druggability labels or published disease-gene lists enter these rankings. Network propagation gives an unsigned priority, not an activation/inhibition direction. For all methods, the predicted expression direction is the independently saved sign of the discovery disease coefficient.

## 6. Evaluation endpoints and uncertainty

Primary descriptive endpoint: **top-20 external direction concordance**, separately for B0, B1 and G1. Freeze the top 20 before evaluation. Report agreeing / evaluable / selected genes, including all missing/zero-sign/failed-fit cases; do not refill lists from lower ranks. For each method also report missing-outcome bounds: agree/20 through (agree + unevaluable)/20. If fewer than 16 of 20 are evaluable, mark the comparison coverage-limited and make no superiority claim. Even above that coverage threshold, do not discuss superiority if the direction of G1-minus-B0 can reverse across these missing-outcome bounds. A direction match can occur for tiny noisy effects and is not biological validation.

Secondary endpoints: top-50 concordance; discovery/external log-fold-change scatter; directional external effect strength sign(discovery coefficient)*external t; and the number with same direction plus external BH-adjusted q<0.05. External q values use the complete externally tested gene family, not just selected genes. Report effect sizes and uncertainty even when q is not small. A lack of significance does not establish no effect.

Primary method contrast is G1 minus B0 at k=20. B1 reveals topology-related behavior; it is not an independent therapeutic truth set. Treat these endpoints as transcriptomic replication, not target accuracy/precision/clinical utility. Do not calculate AUROC against an invented list of negative targets.

Prespecified sensitivity:

- Fix discovery lists AND discovery expression directions; bootstrap only external donors within disease groups, with 500 resamples. Refit the external model using the same resample for all methods to obtain a paired G1-minus-B0 concordance interval. This describes external sampling uncertainty conditional on the discovered lists/directions; it excludes discovery selection/direction uncertainty and future-laboratory variation. Keep the original external evaluability mask fixed; if a selected gene's refit fails, mark that resample invalid rather than changing the denominator. Log singular/failed fits; if >10% fail, suppress interval claims and report instability.
- Discovery leave-one-donor-out analysis recomputes DE and each ranking; report top-20 retention/Jaccard stability separately. Do not describe fixed-list bootstrap as full pipeline uncertainty.
- Generate 200 label-free null gene sets matched on discovery expression and network degree bins; compare concordance descriptively. Freeze bin rules and random seed 42 before external effects. Gene dependence and reuse mean this is a sensitivity analysis, not 200 independent biological replicates.
- Report cycle-day/phase-ascertainment sensitivity where available, age adjustment availability, 30/100 cell thresholds, and scRNA-versus-snRNA measurement differences. A reversely directed, precise effect is a conflict; a broad interval is inconclusive.

## 7. From reproducibility to research-target evidence

Freeze a dossier pool from the union of B0 and G1 top-20 genes before adding genetic annotations. All pool members remain in the final table, including failures.

Use the 2023 endometriosis GWAS and its publicly provided gene-mapping/QTL evidence as a separate evidence layer. Retain locus, ancestry, source analysis, tissue and mapping method. Nearest gene, positional overlap, eQTL, colocalization and perturbation evidence are distinct categories. A GWAS association is not proof that the nominated gene is causal. No de novo causal inference is claimed in v1. Missing public mapping evidence is 'unavailable', not 'negative'.

Do not validate a ranking against evidence used to construct it. Public annotations that reuse either expression cohort are flagged as dependent; multiple portals pointing to one publication count as one underlying source.

Create up to five dossiers by a frozen deterministic rule: select up to three externally direction-consistent genes with q<0.05, ordered by B0 rank within the pool; add one strongly conflicting case and one missing/inconclusive-evidence case if present. Break ties by stable gene ID. If the supported category is empty, write failure/inconclusive dossiers instead and state that no supported shortlist emerged. Do not force a positive result.

Every dossier records: expression evidence, cell-state specificity, genetic support, contradictions, whether intervention direction is known, tissue/assay applicability, and a research recommendation ('follow up', 'hold', or 'deprioritize'). Recommendations are reasoned judgments, not a single arbitrary weighted target score. No drug name or actionability claim is added without a cited public source.

## 8. What constitutes completion

Scientific delivery: external analysis and failure assessment are reported for all prespecified lists; candidate/negative dossiers are traceable; uncertainty and unresolved evidence are explicit. No requirement that G1 win or a drug target be discovered.

Engineering delivery: inputs, transformations, software and outputs are traceable; meaningful boundary tests and a real Nextflow/container run pass; execution on the user-selected SRV01 server has a receipt. AWS remains an optional deployment guide and has not been executed. These checks demonstrate reproducibility, not therapeutic validity.

A failed Stage 0 can complete a feasibility milestone but **does not complete the proposed cross-study scientific project**. Record remaining work honestly and revise before launching a substitute analysis.
