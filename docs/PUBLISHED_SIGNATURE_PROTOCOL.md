# Published bulk disease signatures in the single-cell stromal results — protocol

Status: rules frozen before any single-cell statistic is read for the signature genes (freeze
receipt [`published_signature_test_freeze.json`](../receipts/published_signature_test_freeze.json)).
Written 2026-11-02 after the v1 results, the follow-on analyses and the
[bulk program-level replication](BULK_REPLICATION_PROTOCOL.md) were known. The signature sets
are built from locally downloaded publisher supplements before the freeze and pinned by hash;
no set-level or donor-level statistic on the single-cell data had been computed when this
protocol was frozen.

## 1. Purpose

The bulk replication took the frozen single-cell programs and lists into public bulk cohorts.
This analysis runs the opposite direction: it takes the **published whole-tissue disease
signatures** of eutopic endometrium in the proliferative phase and asks whether they are
visible, with their published directions, in the two single-cell stromal cohorts of v1
(22 + 22 proliferative donors; surgical controls). A tissue-level signature that is
composition-driven or confined to other compartments should not transfer into stroma; one that
reflects within-stroma change should. The estimand is the published gene direction and the
donor-level signature score in the stromal lineage, under the v1 primary design.

## 2. Signature sets (built and published before the freeze)

| Set | Definition and direction convention | Source (local supplement, receipted) | Genes published / mapped |
|---|---|---|---|
| `polineto2020_proliferative_stage_I_II` | Dataset 2A: differentially expressed genes, eutopic endometrium of stage I–II endometriosis versus healthy controls, proliferative phase; limma meta-analysis of GSE4888, GSE6364, GSE7305 and GSE51981 at adjusted p < 0.05 and twofold; logFC = endometriosis minus healthy | Poli-Neto et al. 2020, *Sci Rep*, `41598_2019_57207_MOESM3_ESM.xlsx` (Europe PMC supplementary bundle, SHA-256 in the build receipt) | 4,783 symbols / 4,757 genes (2,097 up, 2,660 down in disease) |
| `polineto2020_proliferative_stage_III_IV` | Dataset 2D: as above for stage III–IV | same article, `..._MOESM6_ESM.xlsx` | 238 / 237 (64 up, 173 down) |
| `bunis2022_proliferative_disease_vs_control` | Supplementary Table 1, sheet `DvC_PEsamples`: proliferative-phase disease versus control in GSE51981 (limma on Brainarray Entrez probe sets, adjusted p < 0.05 and twofold); the published logFC is control over disease, so the direction is its negative; probe sets collapsed to genes, genes with conflicting probe directions dropped | Bunis et al. 2022, *Front Immunol*, `Table_1.xlsx` | 2,570 symbols / 2,547 genes (667 up, 1,880 down) |
| `consensus_proliferative` | genes in both the stage I–II and the Bunis set with the same published direction | derived | 1,630 of 1,630 shared genes agree in direction (390 up, 1,240 down) |

Symbols are mapped to Ensembl gene IDs with the pinned HGNC complete set
(`hgnc_complete_set_2026-10-06.tsv`, SHA-256 8bf6f686…): approved symbol first, then a unique
previous or alias symbol; 49 unmapped symbols are listed in
[`unmapped_or_conflicting.tsv`](../reports/published_signature_test/unmapped_or_conflicting.tsv).
The direction convention of the Poli-Neto tables is confirmed by their consistency with the
Bunis table after its documented flip (FOS, FOSB and JUNB higher in disease in both) and by
the perfect agreement of the 1,630 shared genes. Build:
[`published_signature_build.py`](../scripts/published_signature_build.py), receipt
[`published_signature_build.json`](../receipts/published_signature_build.json), table
[`signature_sets.tsv`](../reports/published_signature_test/signature_sets.tsv).

The three published sets are not independent of each other (Bunis reanalyzes GSE51981, which is
the largest input of the Poli-Neto meta-analysis) and all contrast disease with **healthy**
controls, the contrast that the bulk replication showed to dominate the GSE51981 proliferative
signal. The single-cell cohorts use surgical controls. This difference is part of what the
analysis measures.

Not available at the freeze and therefore not in this version: the Tamaresis et al. 2014
disease-versus-other-pathology list (Supplemental Table 5, publisher site only) and the Burney
et al. 2007 proliferative-phase list (supplemental tables, not open access). Either may be
added later only by a separately frozen, prespecified extension.

## 3. Family, statistics and models

- **Family for T1 and T2**: the 15,435 genes of the [two-cohort combined bound](COMBINED_BOUND_PROTOCOL.md)
  with the pooled inverse-variance t (primary) and each cohort's t (secondary), read from the
  published per-gene table; expression quintile per gene from the discovery MDE table. Each set
  is intersected with the family; sets with fewer than 10 family genes are reported, not labelled.
- **T1 — unsigned enrichment** (per set): one-sided Wilcoxon rank-sum test that |t| of set genes
  exceeds |t| of the other family genes, and the percentile of the set's mean |t| among 10,000
  expression-quintile-matched random sets (seed 20261102). Identical to the genetic set test.
- **T2 — directional concordance** (per set; the primary transfer statistic): signed mean
  statistic S = mean over set genes of (published direction × pooled t), compared with 10,000
  matched random sets that receive the same direction labels bin by bin (percentile and
  two-sided matched p); sign agreement share between the published direction and the pooled
  estimate with an exact binomial test against 0.5 and the matched reference; Spearman
  correlation of the published logFC with the pooled t (descriptive).
- **T3 — donor-level signature score** (per cohort): pseudobulk log-CPM per donor,
  log2((count + 0.5) / (effective library size + 1) × 10⁶), from the positive-control count
  matrices and the published normalization tables, restricted to the 22 proliferative donors of
  the v1 primary design and to signature genes passing the cohort's published expression
  filter; score per donor = mean over signature genes of direction × (log-CPM − gene mean),
  in log2 units (1 = every signature gene shifted by one log2 unit in the published disease
  direction); ordinary least squares on the v1 primary design (intercept, disease, centered age,
  batch columns, adenomyosis; constant columns dropped), disease coefficient with 95% interval,
  two-sided p and 80%-power MDE (exact non-central t). Sensitivities: logFC-weighted projection
  score; leave-one-donor-out range of the coefficient; AUC of the score for disease. None
  changes a label.

## 4. Fixed labels

| Test | Label | Rule |
|---|---|---|
| T1 per set | `enriched` / `no_enrichment` | Wilcoxon p < 0.01 and matched percentile ≥ 99 |
| T2 per set | `concordant` | S percentile ≥ 99, agreement share > 0.5 and binomial p < 0.01 |
| | `discordant` | S percentile ≤ 1, agreement share < 0.5 and binomial p < 0.01 |
| | `no_detectable_concordance` | otherwise |
| T3 per cohort | `difference_detected_toward_published_direction` / `..._against_published_direction` | two-sided p < 0.025, sign of the coefficient |
| | `no_detectable_difference` | otherwise, with interval and MDE |
| T3 pair | `replicated` | both cohorts p < 0.05 with the same sign |
| | `not_replicated` | otherwise |
| any | `not_labelled` | fewer than 10 family genes, or a model not estimable |

## 5. What this establishes and does not

A concordant and replicated signature would show that part of the published whole-tissue
disease signal lives within the stromal lineage and survives the change from healthy to
surgical controls; a null bounds the stromal share of that signal in score units with an MDE.
Whole-tissue signatures mix compartments and composition; the published sets were selected at
tissue level in cohorts that overlap each other; the single-cell cohorts are small; nothing here
nominates a target or validates a signature clinically.

## 6. Order of operations

Download supplements on the analysis host with receipts → build and publish the sets
(`published_signature_build.py`, build receipt) → freeze (this protocol, the configuration with
LF-normalized hashes of every repository input and raw hashes of the two host count matrices,
code, tests; refuse if outputs exist) → one run of `published_signature_test_analysis.py` on
the analysis host in the v1 runtime (Python 3.10, NumPy 2.2.6, SciPy 1.15.2) → published tables,
summary, results receipt and report. Defects found afterwards are corrected in a separately
frozen run; earlier outputs stand as reported.
