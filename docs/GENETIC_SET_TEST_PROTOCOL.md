# Genetic support and the stromal disease contrast — set-level test protocol

Status: rules frozen before any set statistic is computed (freeze receipt
[`genetic_set_test_freeze.json`](../receipts/genetic_set_test_freeze.json)). Written
2026-11-02. The v1 per-gene results, the pooled two-cohort estimates and the published
genetics tables were known; no set-level statistic had been computed. The gene sets are
built from external sources before the freeze and pinned by hash; the test statistics are
read only after it.

## 1. Purpose

The v1 genetics layer annotated the 40 frozen candidate genes and found one limited mapping.
It answered a per-gene question and could not say whether the stromal disease signal, though
below the per-gene threshold everywhere, is **concentrated among genes implicated by
endometriosis genetics**. This analysis asks that set-level question with fixed gene sets
and fixed statistics.

## 2. Gene sets (defined before any statistic is read)

| Set | Definition | Source |
|---|---|---|
| `rahmioglu_nearest` | unique `Nearest gene` symbols of the 2,733 credible-set variants (Supplementary Table 11), excluding "." | Rahmioglu et al. 2023, publisher supplement, SHA-256 01cc0de4… (pinned in `manifests/genetic_sources.tsv`) |
| `rahmioglu_smr` | union of genes in the SMR tables: Supp14 (`Gene`, uterus GTEx), Supp15 (`Gene`, endometrium expression), Supp16 (`Gene`, blood expression), Supp17 (`Closest Gene`, endometrium methylation) | same |
| `opentargets_genetic_any` | Open Targets Platform targets associated with endometriosis (MONDO_0005133) with a `genetic_association` datatype score > 0, direct associations only | Open Targets Platform API, data version recorded at fetch; fetched JSON saved and hashed, derived table published |
| `opentargets_genetic_0.5` | as above with `genetic_association` score ≥ 0.5 | same |

Symbols are mapped to Ensembl gene IDs with the pinned HGNC complete set
(`hgnc_complete_set_2026-10-06.tsv`, SHA-256 8bf6f686…): approved symbol first, then a
unique previous or alias symbol; ambiguous or unmapped symbols are listed and dropped. Each
set is then intersected with the test family; set sizes before and after are reported.
Sets with fewer than 10 family genes are reported but not labelled.

## 3. Family and statistics

- Family: the 15,435 genes of the [two-cohort combined bound](../reports/combined_bound_results.md)
  (estimable in both v1 fits), with the pooled |t| as the primary per-gene statistic and
  each cohort's |t| as secondary statistics, read from the published per-gene tables
  pinned by hash. Expression quintile per gene is taken from the discovery MDE table.
- **Primary test per set**: one-sided Wilcoxon rank-sum test that |t| of set genes exceeds
  |t| of the other family genes (competitive), and an **expression-matched reference**:
  10,000 random gene sets of the same size drawn with the same expression-quintile
  composition (fixed seed 20261102), giving the percentile of the observed mean |t| among
  matched sets.
- Secondary: the same two tests on each cohort's |t|; the share of set genes with the same
  sign of the disease effect in both cohorts against the matched sets (descriptive); the set
  genes listed with their pooled estimates and intervals.

## 4. Fixed labels

| Label | Rule |
|---|---|
| `enriched` | Wilcoxon one-sided p < 0.01 **and** matched-reference percentile ≥ 99 (the 0.01 level accounts for the four sets) |
| `no_enrichment` | otherwise |
| `not_labelled` | fewer than 10 family genes in the set |

An `enriched` label would mean that genetically implicated genes carry larger stromal
disease statistics than expression-matched genes; it would not make any gene a target and
would not identify which genes drive it.

## 5. What this does and does not establish

It tests whether the sub-threshold stromal disease signal is concentrated among genes with
genetic support, under rules fixed in advance. It inherits the sources' definitions
(nearest gene is a positional proxy; Open Targets scores aggregate heterogeneous evidence
and may include the 2023 study itself), uses expression matching as the only matching
variable, and cannot separate composition from within-state effects. No causal, directional
or therapeutic claim follows.

## 6. Order of operations

Build and publish the sets (fetch Open Targets, parse the supplement, map symbols; record
hashes) → freeze (hash protocol, configuration, code, tests; pin inputs and sets; refuse if
outputs exist) → one run of `scripts/genetic_set_test_analysis.py` (Python 3.10, NumPy
2.2.6, SciPy 1.15.2, local) → published tables, summary, receipt and report. Defects found
before any value is read are handled by a superseding freeze; afterwards the results stand
and a correction is published separately.
