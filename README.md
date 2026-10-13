# EndoTarget Evidence

**Public single-cell evidence for prioritizing endometriosis research targets.**

Status: **M0 IMPLEMENTATION IN PROGRESS — 2026-10-13.** The project runs on SRV01 under /home/analyst/te. Bounded public-metadata retrieval and audit implementation have started. Scientific readiness remains NOT PASSED: no expression matrices analyzed, target rankings computed, container workflow executed, or AWS job run. This server has a local Git repository; no GitHub repository has been published.

[中文说明](README.zh-CN.md) · [Research protocol](docs/RESEARCH_PROTOCOL.md) · [Data sources and gates](docs/DATA_SOURCES.md)

## Research question

Which disease-associated expression changes in **eutopic endometrial stromal cells** reproduce in an independent study, and which genes warrant further research when public genetic evidence, contradictory findings, and intervention uncertainty are considered?

The intended decision is a short, justified research shortlist. Expression differences alone do not establish a causal target, an inhibition/activation strategy, or therapeutic efficacy.

## Why this project

This project connects reproductive biology and omics research with computational target assessment. Its planned contribution is a transparent cross-study analysis, a controlled graph-method comparison, and evidence dossiers that a biological research team can inspect. It is a retrospective public-data reanalysis, not a claim of a newly invented algorithm or a blinded prospective study.

## Bounded first release

- One disease, one primary stromal lineage, one discovery study and one independent replication study.
- Donor-level differential expression and explicit menstrual-phase/treatment controls.
- Three ranking baselines: expression evidence, network degree, and network propagation.
- Public human-genetic evidence as a separate evidence layer, followed by 3–5 candidate dossiers; report fewer if the evidence does not justify that many.
- Nextflow, containers, a documented SRV01 execution, optional AWS portability instructions, and an English report for biologists are implementation deliverables. Actual execution was moved to SRV01 by the user; it supplies no AWS execution evidence.

```mermaid
flowchart LR
    A[Public counts and donor metadata] --> B[Data and independence gates]
    B --> C[Discovery: donor-level stromal analysis]
    C --> D[Freeze gene rankings and evaluation]
    D --> E[Independent cohort replication]
    E --> F[Genetic evidence and contradictions]
    F --> G[Research target dossiers]
```

## Candidate inputs

| Role | Source | Current evidence |
|---|---|---|
| Discovery | [GSE266265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265) | Public endometrial scRNA-seq; catalog reports 60 donors. Large RDS import and count-layer checks remain open. |
| Replication | [HECA / E-MTAB-14039](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039), newly generated nuclei | Public author portal exists. The eligible donor subset, raw-count layer and provenance must pass audit. |
| Genetic context | [GCST90205183](https://www.ebi.ac.uk/gwas/studies/GCST90205183) and the originating study's public gene-mapping evidence | Association-to-gene mapping must retain its evidence type. No de novo MR/colocalization in v1. |

These are candidate inputs, not an analysis-ready manifest. The complete audit and known limitations are in [DATA_SOURCES.md](docs/DATA_SOURCES.md).

## What the graph comparison measures

Can network propagation rank genes whose disease-associated expression changes transfer to another cohort better than expression-only and connectivity-only rankings? The endpoint measures **cross-study transcriptomic reproducibility**, not therapeutic target accuracy. Network propagation is a graph baseline; it is not a trained GNN or a foundation model.

## Results and reproducibility

**No scientific results are available yet.** The SRV01 implementation plan, [resource gate](docs/RESOURCE_FEASIBILITY.md), runtime status, and [preliminary scientific readiness](reports/m0_scientific_readiness.md) distinguish actual checks from pending analysis. Metadata/software test results do not validate a therapeutic target.

Planned evidence: donor/exclusion flowchart; cross-cohort effect plots; three-method comparison; sensitivity/failure analysis; 3–5 evidence dossiers; a two-page biology brief; exact input/software/run receipts.

## Executed metadata-audit entry point

On SRV01 only:

```bash
cd /home/analyst/te
python3 -m unittest discover -s tests -v
python3 scripts/m0_metadata.py
```

The real audit and 22 boundary tests passed. Validated cached files are reused with original retrieval provenance. See [audit summary](reports/m0_metadata_summary.md) and [readiness JSON](reports/m0_metadata_readiness.json). Successful metadata execution leaves Stage 0 OPEN until clinical, cell, count and resource gates pass.

## Repository contents today

```text
README.md / README.zh-CN.md
docs/{RESEARCH_PROTOCOL,DATA_SOURCES}.md
manifests/source_inventory.csv
templates/{TARGET_DOSSIER,BIOLOGY_BRIEF}.md
receipts/metadata_audit_2026-10-13.json   # actual small metadata audit only
```

This is the working project root on SRV01. The listing above preserves the original design package; implementation records additionally live in scripts/, tests/, reports/ and receipts/. Scientific modules remain gated; no empty success placeholders are provided.

## Data use and attribution

All scientific inputs must be publicly obtainable without private-data agreements or new experiments. Download from original providers; do not commit patient-level matrices or redistribute third-party resources by default. Record source terms individually. A publication's open-access license does not automatically establish the license of every dataset or dependency. A code license must be selected before software publication; this design does not relicense source material.

See [the source register](docs/DATA_SOURCES.md) for original studies.
