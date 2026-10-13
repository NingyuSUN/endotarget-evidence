# Pre-effect method decisions — 2026-10-13

Status: method decisions adopted before computing disease effects; the complete M0/input freeze has **not** passed. No expression rankings or external DE were used. This is a retrospective, auditable analysis, not a preregistered or blinded study.

## 1. One phase rule across both studies

Primary phase means **author-adjudicated proliferative phase**, after explicit exclusion of author ambiguous/periovulatory, secretory, menstrual, hormone-treatment and unknown categories. Retain the original histology, clinical day/progesterone and expression-assisted phase fields. Do not infer or relabel phase from our expression analysis.

This clarification uses the same evidence hierarchy in discovery and replication. Both sources combine or adjudicate clinical/histological and transcriptomic phase information. A primary rule using the final author label is reproducible across them; it does not establish independent histological confirmation. Histology conflicts remain visible and support a stricter predeclared sensitivity. They are not silently corrected.

Discovery: exclude the five explicitly flagged clinical-exclusion samples81,108,115,119,121 and quarantine sample73 because non-ENDO conflicts with Endometrioma=1. The author-proliferative clinical subset then has22 candidates (12 ENDO,10 non-ENDO). The two remaining histology-phase conflicts (samples76,123) are retained in this primary source-defined phase scenario and excluded in the stricter sensitivity, yielding20 (11,9). No disease label is changed.

HECA: use the H5 author genotype/Donor crosswalk and author clinical table, not SDRF individual values that disagree with actual demultiplexed annotations. The initial author-proliferative subset has24 donor IDs (16 ENDO,8 controls), before cell threshold and remaining clinical/design checks. Applying strict histology-conflict exclusion consistently would leave fewer than5 controls, so that external sensitivity fails the fixed minimum. It must be reported as not estimable; do not relax its threshold or call primary results robust to it.

The original fixed sample minima, >=50 primary-lineage cells/nuclei, identical-phase restriction and final design rank/residual-df rules remain in force. The former20-donor discovery set was a proposed scenario, not a frozen primary analysis. This decision was made after examining metadata eligibility but before any disease effects; disclose that timing.

## 2. Hormone evidence remains source-qualified

Discovery has explicit FALSE in the archived source clinical hormone/IUD field, with the three-month exclusion interval specified in the linked primary paper.

HECA donor treatment NA and H5 literal nan remain missing values. For an author-defined natural-cycle phase group, use `hormone_eligibility_basis=study_defined_natural_cycle_group`, supported by the paper definition requiring at least three months without hormones and its separate Hormones group. This is source-derived cohort eligibility, not an invented per-person negative or stop date. Explicit contradictory treatment, Hormones and Unknown categories fail primary eligibility. Preserve source fields and run the coding consistency audit.

## 3. Clinical covariates and designs

In discovery, start from age + capture batch + adenomyosis + disease. Use the actual complete-case final matrix; require full rank and residual degrees of freedom>=5. Sequencing batch is nested in capture batch and must not be added as a redundant extra set of columns. Capture batches without both disease groups do not create within-batch disease comparisons; report mixed capture batches12,14,16 and18 in the22-donor author-phase set, versus only14 and16 in the20-donor strict-phase sensitivity. Metadata-only checks cannot show power or eliminate residual confounding.

Include adenomyosis as a prespecified adjustment when it varies and is identifiable. The additive disease effect remains an association and does not erase all pathology differences. An adenomyosis-free sensitivity is a separate cohort; recheck donor minima and design, and omit a now-constant adenomyosis column rather than fit a singular model. Do not conflate the18-donor author-phase sensitivity (8 cases,10 controls) with the16-donor strict-phase sensitivity (7 cases,9 controls; fails minimum).

The external design still requires its own batch/repeated-library audit; it is not accepted merely because discovery passes algebraic checks. No model may be selected using DE significance.

## 4. Compact discovery counts are acceptable with explicit provenance

Accept the public, documented donor aggregation as a candidate counts source, not as a fresh reconstruction of the complete original RNA matrix. It retains unnormalized, unrounded fractional Alevin estimates. Per-donor totals and cell annotation sums agree in the inspected exports; archive CRC/SHA and aggregation code have been audited. Original RDS reconstruction has not been independently rerun here.

Exclude all50 irreversible multi-source canonical groups. Recover original IDs for the34,550 singleton groups before further annotation. These include9,349 versioned Ensembl IDs and25,201 source symbols/other IDs, not34,550 stable Ensembl IDs. Per-donor estimated-UMI loss is about0.093-0.162%; this small total fraction does not prove that removed genes lack biological importance. Record exclusions and normalize consistently after the declared feature filtering.

A frozen human annotation must still map accepted features without ambiguity to common stable gene IDs. Unmapped, ambiguous and duplicate-target mappings require an explicit pre-effect rule. This stable-ID gate, external counts-layer audit, final lineage mapping, final clinical/design ledger and all ranking/network freezes remain outstanding.

## 5. Primary lineage boundary

Do not treat all mesenchymal labels as endometrial stromal fibroblasts. HECA primary candidates are eStromal_MMPs, eStromal, eStromal_cycling, dStromal_early, dStromal_mid, dStromal_late and Fibroblast_basalis. Exclude ePV subtypes and uSMCs (mural/myometrial); exclude HOXA13. The primary HECA paper describes HOXA13+ mesenchymal cells as likely cervical contamination, supported by a cervix reference comparison. Inclusion must not be changed to rescue FX1215 or FX1141, which fail50 cells under the anatomically restricted labels. Discovery's P077 fibroblast definition still needs an explicit source-label crosswalk before final acceptance.

YELLOW: understanding not yet assessed. Defend author-adjudicated phase versus independent clinical confirmation, additive confounder adjustment versus causal inference, and reversible singleton mapping versus irreversible gene merging. No new user approval is required for these transparent pre-effect methodological decisions within the authorized project.
