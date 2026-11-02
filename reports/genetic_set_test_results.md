# Genetic support and the stromal disease contrast — set-level results

**Labels under the frozen rules: `no_enrichment` for all four gene sets.** Genes implicated
by endometriosis genetics do not carry larger stromal disease statistics than
expression-matched genes in the two v1 cohorts, pooled or separately. For the 83 credible-set
nearest genes of Rahmioglu et al. 2023 present in the family, the mean pooled |t| is 0.86
against 0.96 for expression-matched random sets (percentile 11; one-sided Wilcoxon p 0.95);
for the 289 Open Targets genetic-association targets it is 0.90 against 0.96 (percentile 6;
p 0.98). The sets' genes agree in sign across cohorts no more often than matched genes
(43–50% versus 53%).

Completed 2026-11-02 under the [protocol](../docs/GENETIC_SET_TEST_PROTOCOL.md) frozen the
same day after the sets were built and before any set statistic was computed (freeze
[`genetic_set_test_freeze.json`](../receipts/genetic_set_test_freeze.json), git head 3b6f84f;
results receipt [`genetic_set_test_results.json`](../receipts/genetic_set_test_results.json)).
Family and statistics: the 15,435 genes of the [two-cohort combined bound](combined_bound_results.md)
with the pooled |t| (primary) and each cohort's |t| (secondary); expression quintile per gene
from the discovery MDE table. Run once locally (Python 3.10.21, NumPy 2.2.6, SciPy 1.15.2).

## Gene sets

| Set | Source | Genes in source | Mapped Ensembl genes | In the family |
|---|---|---:|---:|---:|
| `rahmioglu_nearest` | nearest genes of the 2,733 credible-set variants, Supplementary Table 11 of Rahmioglu et al. 2023 (pinned publisher supplement) | 102 symbols | 100 | 83 |
| `rahmioglu_smr` | SMR genes, Supplementary Tables 14–17 (uterus GTEx, endometrium and blood expression, endometrium methylation) | 20 symbols | 18 | 14 |
| `opentargets_genetic_any` | Open Targets Platform 26.09, endometriosis (MONDO_0005133), direct associations with a genetic-association score > 0 | 356 | 356 | 289 |
| `opentargets_genetic_0.5` | as above, score ≥ 0.5 | 55 | 55 | 48 |

Symbols were mapped with the pinned HGNC complete set (approved, then unique previous or
alias symbol); four symbols did not map (`LOC100506236`, `QARS`, `HCG4P7`, `Metazoa_SRP`) and are
listed. Build receipt: [`genetic_set_build.json`](../receipts/genetic_set_build.json) (fetch
time, Open Targets data version, raw-response hash); published tables:
[gene sets](genetic_set_test/gene_sets.tsv), [Open Targets targets](genetic_set_test/opentargets_endometriosis_targets.tsv).

## Tests

Primary statistic pooled |t|; one-sided Wilcoxon rank-sum (set > rest) and the percentile of
the set's mean |t| among 10,000 random sets with the same expression-quintile composition
(seed 20261102). Label `enriched` required p < 0.01 and percentile ≥ 99.

| Set | n | Mean \|t\|, set | Mean \|t\|, rest | Matched reference mean | Matched percentile | Wilcoxon p (one-sided) | Label |
|---|---:|---:|---:|---:|---:|---:|---|
| `rahmioglu_nearest` | 83 | 0.86 | 0.96 | 0.96 | 11.2 | 0.95 | `no_enrichment` |
| `rahmioglu_smr` | 14 | 0.96 | 0.96 | 0.96 | 51.1 | 0.42 | `no_enrichment` |
| `opentargets_genetic_any` | 289 | 0.90 | 0.96 | 0.96 | 5.9 | 0.98 | `no_enrichment` |
| `opentargets_genetic_0.5` | 48 | 0.92 | 0.96 | 0.95 | 38.8 | 0.70 | `no_enrichment` |

Per cohort (secondary): discovery matched percentiles 15, 36, 29 and 19; HECA 14, 64, 39 and
44; no Wilcoxon p below 0.30. Sign agreement across cohorts among set genes: 47%, 43%, 50%
and 48%, against 53–54% for matched sets (percentiles 10, 13, 7 and 18). Full table:
[`set_results.tsv`](genetic_set_test/set_results.tsv); set genes with their pooled estimates:
[`set_members.tsv`](genetic_set_test/set_members.tsv).

The set gene with the smallest pooled q is SKIDA1 (ENSG00000180592; in the nearest-gene and
Open Targets sets): pooled β −0.85 log2, q 0.21, same sign in both cohorts, the seventh largest
|t| in the combined family. It is the one frozen-pool gene for which the v1 genetics layer
found a limited mapping. This is a descriptive coincidence of two sub-threshold observations,
reported because the protocol lists the set genes; it nominates nothing.

## Readings

- **The sub-threshold stromal disease signal is not concentrated among genetically
  implicated genes.** Their statistics are, if anything, slightly smaller than those of
  expression-matched genes, and their cross-cohort sign agreement is at chance.
- **This does not contradict the genetics.** Risk loci act in tissues, cell types, phases and
  on scales (regulatory, developmental, lesion) that a proliferative-phase eutopic stromal
  pseudobulk of 22 + 22 donors may not capture; the test says only that in this measurement,
  under this design, the genetic sets look like matched genes.
- **The genetics layer of the project now has a set-level statement** to go with its
  per-gene annotation: no supported shortlist, and no enrichment of the disease signal among
  GWAS-implicated genes.

## What this establishes and does not

Under rules fixed before computation, four genetically defined gene sets show no enrichment
of stromal disease statistics in the pooled or single-cohort v1 fits. Nearest-gene and Open
Targets sets are positional and aggregate proxies, Open Targets evidence may include the
2023 study itself, expression is the only matching variable, and the statistics cannot
separate composition from within-state change. No causal, directional or therapeutic claim
follows.

## Files

[`reports/genetic_set_test/`](genetic_set_test/): set results, set members, gene sets,
Open Targets targets, unmapped symbols, summary. Code:
[`scripts/genetic_set_build.py`](../scripts/genetic_set_build.py),
[`scripts/genetic_set_test_analysis.py`](../scripts/genetic_set_test_analysis.py),
[`scripts/followon_freeze.py`](../scripts/followon_freeze.py), [tests](../tests/test_genetic_set_test.py).
