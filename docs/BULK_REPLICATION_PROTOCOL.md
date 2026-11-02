# Program-level replication in public bulk eutopic endometrium — protocol

Status: rules frozen before any expression value is read (freeze receipt
[`bulk_replication_freeze.json`](../receipts/bulk_replication_freeze.json)). Written 2026-11-02
after the single-cell results (v1, positive control, MDE, secretory program, combined bound,
composition and genetic set tests) were known. Only GEO metadata (sample characteristics, file
listings) had been read for the bulk series when this protocol was frozen; the strata below
come from that [metadata audit](../manifests/bulk_replication_GSE51981_samples.tsv).

## 1. Purpose

The single-cell analyses bound proliferative-phase stromal disease differences and the
secretory stromal program's disease response in 22 + 22 and 35 + 42 donors. Public bulk
microarray cohorts of eutopic endometrium are several times larger. This analysis takes the
**frozen single-cell programs and lists** into those cohorts and asks, with fixed estimands and
labels, whether what the single-cell data could not resolve is visible at tissue level in
hundreds of donors. The estimand changes from stromal lineage to whole tissue; that change is
stated, not hidden.

## 2. Cohorts (metadata audited, expression not yet read)

| Series | Platform | Processed data | Eutopic samples by phase and group (GEO characteristics) | Role |
|---|---|---|---|---|
| GSE51981 (Tamaresis et al. 2014) | GPL570, Affymetrix HG-U133 Plus 2.0 | authors' series matrix (log2, GCRMA per the paper) | proliferative: 29 endometriosis (12 minimal/mild, 17 moderate/severe), 15 other uterine/pelvic pathology, 20 no pathology; early secretory 18 / 6 / 6; mid-secretory 28 / 14 / 8; 4 of indeterminate phase excluded | **primary** |
| GSE141549 (Gabriel et al. 2020) | GPL10558 and GPL13376, Illumina HT-12 v4 and WG-6 v2 | authors' batch-corrected combined normalized matrix (`GSE141549_batchCorrectednormalizedArrayscombined.xlsx`) with `GSE141549_Sample_link.xlsx`; the GEO series matrices carry no expression values | endometrium: controls 7 proliferative, 17 secretory (10 on medication, 8 unknown); patients 17 proliferative (stage 1: 2, 3: 6, 4: 9), 27 secretory, 43 on medication, 7 menstrual, 10 unknown; replicates flagged | audited and downloaded with the other files, **not analyzed in this version**: its sample-to-column linkage lives in workbooks whose structure was not inspected before the freeze, so it is reported as `not_analyzed` and left to a later prespecified extension |
| GSE6364 (Burney et al. 2007) | GPL570 | authors' series matrix | proliferative 6 / 5, early secretory 6 / 3, mid-secretory 9 / 8 (endometriosis / control) | secondary; the founding study, small |

GSE25628 (GPL571, 22 samples, paired eutopic and ectopic tissue) was audited and is not
used: it offers no disease-versus-control eutopic contrast.
Hormonal-treatment samples (GSE141549 "medication") and menstrual-phase samples are excluded.
Technical replicates in GSE141549 are averaged per donor before analysis. Age is used as a
covariate where present (GSE141549); GSE51981 and GSE6364 have no age in GEO.

## 3. Inputs pinned at the freeze and at download

Frozen before download: the series and file URLs, the SOFT metadata files (SHA-256), the
frozen single-cell artifacts used as inputs: the positive-control phase programs
(`reports/phase_positive_control_f1/de_all.tsv.gz`, `reports/phase_positive_control_f2/heca_de_all.tsv.gz`),
the frozen B0 top-20 and top-50 lists and directions (`reports/replication_v2_dossier_pool.tsv`
and the ranking tables of the discovery freeze), the pooled per-gene table
(`reports/combined_bound/combined_per_gene.tsv.gz`), the HGNC complete set for probe-to-gene
mapping, and the GPL570 / GPL10558 / GPL13376 annotation files (URLs pinned; hashes recorded at
download). At download, a data receipt records every file's SHA-256 before any value is read
by the analysis.

## 4. Preprocessing (fixed)

Expression as provided by the authors (log2 scale; verified to be in the expected range,
otherwise log2-transformed once). Probes mapped to Ensembl genes through the platform
annotation (Entrez ID, then symbol) and the pinned HGNC set; probes without a unique gene are
dropped; one value per gene per sample by the probe with the highest mean expression in the
cohort. Genes are centered per cohort on the mean of its natural-cycle proliferative
non-endometriosis samples. No batch correction across series; each series is analyzed
separately.

## 5. Estimands and models

**E1 — secretory program disease response (tissue level).** Program score per sample =
β-weighted projection on the frozen phase program (β from the other single-cell cohort's
positive-control fit, q < 0.05 genes, as in the secretory program protocol; both the
HECA-derived and the discovery-derived programs are scored, HECA-derived primary). Gate: in a
model score ~ phase + disease (+ age), the secretory-versus-proliferative coefficient must be
positive with AUC ≥ 0.9; otherwise `NOT_INTERPRETABLE`. Primary θ = secretory × endometriosis
interaction (all natural-cycle samples; secretory = early plus mid secretory); secondary δ =
disease difference within secretory samples; sub-phase (early versus mid) as a covariate in a
sensitivity. Ordinary least squares, t tests, 95% intervals, 80%-power MDE of θ.

**E2 — direction transfer of the frozen single-cell list.** In proliferative samples, per-gene
disease coefficient by ordinary least squares with empirical-Bayes moderation (limma-style,
implemented in Python on the gene-level log2 matrix; endometriosis versus all
non-endometriosis; sensitivity: versus other-pathology controls only, the v1-like contrast).
Same-direction count of the frozen B0 top-20 against its discovery signs (the published
dossier pool carries the list and signs; the top-50 list is not published and is not used),
with missing-outcome bounds, and the percentile of the count against 10,000
expression-matched random 20-gene sets drawn from the bulk family (the v1 matched sets are
not published and are not reused).

**E3 — tissue-level per-gene contrast (descriptive context).** The limma proliferative
contrast of E2: number of genes at BH q < 0.05, median MDE, interval-exclusion counts, and the
Spearman correlation and sign agreement between bulk t and the pooled single-cell t over
shared genes.

## 6. Fixed labels

| Estimand | Label | Rule |
|---|---|---|
| E1 per series | `attenuation_detected` / `enhancement_detected` | θ sign and two-sided p < 0.025 |
| | `no_detectable_difference` | otherwise, with interval and MDE |
| E1 across series | `replicated` | primary series and at least one secondary series with p < 0.05 and the same sign |
| E2 per list | `transfers` | ≥ 16 / 20 (≥ 40 / 50) same direction and matched percentile ≥ 95 |
| | `does_not_transfer` | otherwise |
| E3 | descriptive only | no label |

## 7. What this establishes and does not

A detected E1 effect would show that the secretory stromal program is disease-responsive at
tissue level in cohorts large enough to see it, and bound its size; a null bounds it in
hundreds of donors. E2 tells whether the single-cell lists' directions hold in whole tissue.
Whole-tissue expression mixes stromal, epithelial, immune and vascular compartments and
reflects composition as well as within-cell change; microarray platforms and author
normalizations differ from the count-based single-cell estimator; disease ascertainment and
controls differ by series (GSE51981 distinguishes other pelvic pathology from no pathology;
GSE141549 controls are healthy). Nothing here nominates a target.

## 8. Order of operations

Freeze (hash this protocol, the configuration, code and tests; pin metadata hashes and URLs;
refuse if outputs exist) → download expression matrices and platform annotations on the
analysis host and write the data receipt → one run of `scripts/bulk_replication_analysis.py`
(Python 3.10, NumPy 2.2.6, SciPy 1.15.2 of the v1 runtime) → published tables, summary,
results receipt and report. Defects found before any value is read are handled by a
superseding freeze; afterwards results stand and corrections are published separately.
