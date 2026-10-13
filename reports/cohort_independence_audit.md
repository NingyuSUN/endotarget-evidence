# Public cohort provenance audit

Status: source-supported independence for the **final 22 discovery donors and 22 HECA new-nuclei donors**. Selected donor/sample provenance joins are executed and recorded. This is not genotype-confirmed participant non-overlap.

The discovery cohort is GSE266265, the Swiss-approved study (Swiss Ethics Committee KEK-BE 01780, 2019). Its paper describes consent before planned laparoscopic procedures, endometrial pipelle biopsy and surgical/histopathological disease ascertainment. Bern institutional involvement is documented, but the inspected sample-methods paragraph does not name every recruiting clinic or provide all recruitment dates. No dates are inferred from the ethics identifier.

The external cohort is specifically the **newly generated HECA nuclei** associated with E-MTAB-14039, not the complete integrated atlas. The paper assigns newly generated superficial endometrial samples to ENDOX, FENOX, Sanger Human Cell Atlasing, and Immunology and Subfertility, with separate UK ethics approvals and Oxford surgery/outpatient recruitment. Every selected sample is linked through the actual H5 author-deconvolution donor/sample fields and official author clinical/library CSVs; selected study_original is Marecckovaa. These fields identify the new study component but do not distinguish which recruiting sub-study supplied each person. Do not infer sub-study from EX/FX prefixes.

## Selected provenance and repeated specimens

- Discovery: manifests/final_discovery_clinical_selected.tsv contains 22 unique donors, original sample IDs, 22 GSM accessions, original archived clinical workbook rows/hash, and 150,788 included fibroblast cells. Sample73's disease conflict is excluded by the pre-effect rule. The author-clinical supplement copy is tied to the public OSF release manifest/ZIP; primary-publisher byte equivalence is not independently verified.
- HECA: manifests/final_heca_clinical_selected.tsv contains 22 unique author donors, selected library/sample IDs, clinical/phase/pathology and source fields, and 52,744 included nuclei. All 22 author_deconvolution_Donor values match donor_id. SDRF declared individual identifiers are retained but do not override the author cell-demultiplex join.
- FX1101 and FX1236 each contribute only the eligible library with most frozen primary-lineage nuclei (lexicographic tie rule). Other eligible libraries are excluded; 2,741 nuclei are lost. Unknown technical-replicate versus repeat-biopsy relationships remain unknown, and specimens are not counted as independent donors.
- Clinical eligibility, hormone evidence hierarchy and histology-conflict sensitivity are documented in reports/final_clinical_gate_report.md and docs/PRE_EFFECT_METHOD_DECISIONS.md; this report does not revise those rules.

The public recruitment frameworks support different study provenance; no shared recruiting study has been identified in the inspected methods. Different accessions or disjoint local ID strings alone would not establish independent people. No cross-study genomic/person linkage is available. These limitations do not require private clinical reconfirmation before a transparently source-attributed independent-cohort analysis.

Both studies inherit related Tan/reference annotation methodology. That methodological dependency is explicit: it does not imply shared people, but broad label agreement cannot be presented as independent cell-identity validation. The fixed marker check is a limited descriptive population plausibility audit. HECA cells and nuclei can share participants within HECA; those two assays are not separate replication cohorts here.

## Sources and reproducibility

receipts/cohort_independence_sources.json retains the exact primary XML hashes/recruitment paragraphs and adds both final ledger hashes and original selected IDs. Source primary publications: [discovery](https://doi.org/10.1038/s41467-026-73020-4) and [HECA](https://doi.org/10.1038/s41588-024-01873-w). All inputs were cached; this audit retrieved no new data and used no disease effects or target rankings.

YELLOW: understanding not yet assessed. The user should distinguish independent recruiting studies, unique local labels, shared annotation methods, and actual genotype evidence of participant overlap.
