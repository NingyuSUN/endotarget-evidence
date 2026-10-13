# EndoTarget Evidence

**Public single-cell evidence for prioritizing endometriosis research targets.**

Status: **TWO-COHORT COUNT/DESIGN AUDIT PASSED; DISEASE ANALYSIS NOT YET RUN — 2026-10-13.**
All project downloads and runs are on SRV01 under `/home/analyst/te`. Real donor count matrices, clinical designs, a fixed marker audit and a gene-network reference are available. No target rankings, external disease-effect replication, Nextflow/container run or AWS deployment is claimed. This is a server-local Git repository; it has not been published to GitHub.

[中文](README.zh-CN.md) · [Protocol](docs/RESEARCH_PROTOCOL.md) · [Actual input provenance](docs/INPUT_PROVENANCE.md)

## Research question

Which expression changes in **eutopic endometrial stromal fibroblasts** reproduce across independently recruited endometriosis studies, and which genes merit further research after genetic evidence, contradictions and uncertain intervention direction are considered?

The planned comparison is expression-only ranking, degree-only ranking and network propagation on the same discovery-defined graph universe. The endpoint is cross-study transcriptomic replication. It is not therapeutic target accuracy, a trained GNN or a foundation model, and it does not require the graph method to win.

## Real inputs prepared

| Input | People and primary cells/nuclei | Stable named genes | Evidence |
|---|---|---:|---|
| [GSE266265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265), using documented [secondary OSF donor sums](https://doi.org/10.17605/OSF.IO/YJM7H) | 22 people: 12 ENDO, 10 non-ENDO; 150,788 cells | 26,109 | [Discovery count audit](reports/discovery_counts_preparation.md); original fractional Alevin estimates retained |
| Newly generated [HECA E-MTAB-14039 nuclei](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039) | 22 people: 15 ENDO, 7 non-ENDO; 52,744 selected nuclei | 26,125 | [HECA count audit](reports/heca_counts_preparation.md); all selected barcodes matched to 22 public integer-count archives |
| Human STRING 12.5 direct experimental/database evidence | 159,312 gene-level edges; 14,554 nodes | Background only | [Gene mapping audit](reports/string_gene_reference_audit.md); final expression-induced universe not constructed |

Both cohorts use the author-adjudicated proliferative phase and source-qualified hormone eligibility. The smaller histology-conflict-excluded HECA sensitivity has only three controls and is not estimable under the fixed minimum. Controls may have other gynecological conditions; they are not all healthy volunteers.

There are 21,585 shared stable gene IDs before expression filtering. Missing IDs remain missing; external availability does not filter discovery candidates. The [clinical/design audit](reports/final_clinical_gate_report.md), [source independence ledger](reports/cohort_independence_audit.md) and [fixed marker audit](reports/lineage_marker_audit.md) retain source and interpretation limits. Source-level independence is supported by distinct recruitment frameworks; no genotype identity check or independent cell reannotation is claimed.

## What is verified, and what remains

- Count values, donor order, gene uniqueness, source hashes and model rank were checked on the actual prepared matrices. Discovery has residual degrees of freedom 10; HECA has 12. This is mathematical estimability, not evidence of power or absence of confounding.
- The fixed 23-marker audit is descriptive and does not change eligibility. Markers are not exclusive, and pseudobulk abundance is not cell purity.
- The HECA H5 expression layer was rejected because it is normalized/log-transformed. Public unnormalized matrices provide the accepted counts. A 9,607-UMI excess over the author's saved H5 totals (0.00266% of raw totals) remains recorded; gene filtering explains part, with the full cause not proven.
- All 63 integrated software tests passed on the populated SRV01 project. [Receipt](receipts/m0_integrated_tests.json). This is not a clean container execution or scientific validation.
- Next: freeze the statistical software/filtering implementation, run discovery, freeze the three rankings and directions, then run external replication, sensitivities, genetic evidence dossiers and a biology report.

The original full discovery RDS was not re-imported. The graph maps protein preferred names to HGNC genes; sequence/genome-release equivalence is not established. The analysis is retrospective, not blinded or preregistered.

## Executed checks

These commands apply to the **existing populated SRV01 project**; a clean-install Nextflow entry point is still being implemented.

```bash
cd /home/analyst/te
envs/m0-h5/bin/python -m unittest discover -s tests -v
envs/m0-h5/bin/python scripts/audit_m0_contracts.py
```

The contract command writes a hash snapshot and validates the prepared count/design interfaces. It does not compute disease effects or independently decide all scientific gates. Data, environments, large outputs and caches are ignored by Git.

## Intended delivery

A donor/exclusion flowchart, complete gene-level effect tables, three-method comparison, uncertainty and failure analysis, up to five traceable evidence dossiers, and an English biology brief. Public genetics from [GCST90205183](https://www.ebi.ac.uk/gwas/studies/GCST90205183) remains a separate evidence layer after ranking freeze; no new wet experiment or de novo causal claim is required.

Nextflow and actual container execution remain engineering deliverables. Docker is installed but the current account cannot access its daemon; see the verified access explanation. AWS is only an unexecuted portability guide. No cloud resource has been created.

Source terms are recorded separately; this project does not relicense third-party material or redistribute patient-level matrices by default.
