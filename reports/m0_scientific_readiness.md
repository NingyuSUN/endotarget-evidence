# M0 scientific readiness — preliminary, 2026-10-13

Status: NOT PASSED. Software and catalog checks cannot establish the analysis cohort. No disease effects or target rankings have been computed.

## Verified source-method facts

Public primary-paper XMLs were downloaded to SRV01 and hashed in receipts/reference_metadata.json. The discovery study describes Salmon/AlevinQC v1.6.0 processing and exclusion of exogenous hormones within three months before surgery. This confirms provenance wording, not the numeric contents or integer status of the RNA count layer.

The HECA paper states that naturally cycling donors had no hormonal treatment for at least three months and that recruitment was in reproductive years. It also describes genotype-based donor demultiplexing for mixed libraries, with doublet/unassigned barcodes discarded. These facts support a common study-level hormone window but do not replace per-donor label and cell-crosswalk checks. Some donors had other benign pathology identified later; control indications must be retained.

Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC13385355/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC11387200/ .

## Remaining gate checklist from Astra scientific review

1. Namespace each cell barcode by library, then join author-demultiplexed donor, biopsy and SDRF individual without multiplying rows. File names must not substitute for the declared donor.
2. Resolve clinical phase/hormone labels, source reproductive-age ascertainment, duplicate biopsies and conflicting annotations before disease effects. The reported 6 replication proliferative controls are candidates; losing two fails the >=5 rule.
3. Verify raw/author-filtered count layer, finite/nonnegative values, axes and stable gene IDs. Floating dtype alone does not identify normalized or estimated counts; no blind rounding.
4. Establish selected-donor provenance and source independence. Recheck final donor minima, cell thresholds, full-rank model, residual df>=5 and batch overlap after exclusions.
5. Prove an import path within SRV01 memory/disk constraints. Bounded metadata acquisition is allowed while M0 is open, but disease modeling is not.

## Interpretation limits

No external replication, target validation, Docker execution, Nextflow run or AWS deployment has been completed. A failed M0 may finish a feasibility milestone but does not fulfill the scientific goal. No phase, donor threshold or disease scope has been loosened.

Understanding label: YELLOW — understanding not yet assessed. The key distinction is between independent patients and their many cells, and between reproducible expression associations and causal therapeutic targets.
