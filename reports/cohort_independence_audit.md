# Public cohort provenance audit

Status: source-level separation SUPPORTED; selected-donor audit PENDING. This audit does not assert genotype-confirmed participant non-overlap.

The selected discovery source is the Swiss-approved GSE266265 study (Swiss Ethics Committee KEK-BE 01780, 2019). Its paper describes consent before planned laparoscopic procedures, endometrial pipelle biopsy and surgical/histopathological disease ascertainment. Bern institutional involvement is documented, but the inspected sample-methods paragraph does not itself name every recruiting clinic or provide all recruitment dates. Do not invent dates from the ethics identifier.

The candidate external source is the **newly generated HECA nuclei** associated with E-MTAB-14039. The primary paper assigns its superficial endometrial samples to ENDOX, FENOX, Sanger Human Cell Atlasing, and Immunology and Subfertility. It provides separate UK ethics approvals and describes Oxford surgery/outpatient recruitment. The selected-donor subset still requires the actual author donor/library/cell crosswalk. The complete integrated HECA cell atlas must not be substituted for this source.

These publicly documented recruitment frameworks support separate study provenance. Different accession numbers or distinct local donor-ID strings by themselves would not establish independent participants. No shared recruiting study has been identified between these selected sources in the inspected methods. This is evidence supporting independence, not proof based on cross-study genotyping. HECA cells and nuclei can share individuals within HECA and must not be counted as separate external cohorts.

Before final eligibility is frozen, record for every selected external donor the author donor ID, original library and barcode join, source study and assay, clinical phase/disease, and repeated-library handling. Resolve contradictory identity fields with the author-demultiplexed annotation and explicit source tables. Do not silently choose the clinical label that preserves the planned sample size.

## Sources and reproducibility

The exact source XML hashes and inspected recruitment paragraphs are in `receipts/cohort_independence_sources.json`. Primary publications: [discovery](https://doi.org/10.1038/s41467-026-73020-4) and [HECA](https://doi.org/10.1038/s41588-024-01873-w). These XMLs are existing downloaded inputs, not newly retrieved clinical or expression data.

No external disease effects, target rankings, or genetic target labels were used in this audit. YELLOW: understanding not yet assessed; explain the difference between independent recruiting studies, unique local labels, and actual evidence of participant overlap.
