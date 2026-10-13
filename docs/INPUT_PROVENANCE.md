# Active input provenance and limits

This is a public-data reanalysis. The usable inputs are the compact discovery donor sums and the selected public HECA integer matrices below. The original large discovery RDS was not re-imported; the normalized HECA H5 expression layer is rejected for counts. A count preparation receipt is an engineering record, not an independent biological validation.

| Resource | Actual use | Source/use evidence and boundary |
|---|---|---|
| GSE266265 original study | Discovery study identity, clinical and annotation methods | GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265 and primary DOI https://doi.org/10.1038/s41467-026-73020-4. Original code is not copied from the unlicensed author repository. Original full RDS reconstruction is not claimed. |
| OSF yjm7h | Secondary public P077 donor sums, original-feature mapping, cell annotation membership, archived clinical workbook; P005 comparator for marker audit | https://doi.org/10.17605/OSF.IO/YJM7H ; selected archive members checked against ZIP CRC and public SHA manifest. Original OSF materials declare CC BY 4.0; third-party material retains its source terms. This is a secondary export of original author data, not a new recruitment cohort. Publisher byte equivalence of the archived clinical workbook is not verified. |
| E-MTAB-14039 processed archives | Exactly 22 selected donor/sample integer CellRanger matrices; 52,744 author-selected nuclei | https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039 . Exact public URLs and byte sizes are in reports/heca_counts_raw_inventory.json; final download and count receipts record hashes. These are processed unnormalized counts, not raw FASTQ. Controlled-access reads are not used. No blanket redistribution license is inferred from download access. |
| Official HECA metadata and H5 observation fields | Donor, sample, library, clinical phase and cell annotation linkage | Author repository commit 244c282c20c3efd53707e3eb9f58bbd1cbd56ad7; exact metadata URLs/hashes in receipts/discovery_clinical_gate_heca_downloads.json. H5 source validators and extracted observation hashes remain recorded. The H5 .X normalized/log expression matrix is not an accepted count input. |
| HGNC 2026-10-06 | Named-gene stable identifier crosswalk | Official versioned object and generation in config/gene_mapping_reference.json; CC0 1.0 at https://hgnc.genenames.org/about/license/ . Approved symbols, unique official previous symbols and verified version-stripped ENSG IDs are distinguished. This is not sequence or genome-coordinate equivalence across reference releases. |
| STRING 12.5 human | Direct experimental/database-supported functional association background | Official versioned download URLs in config/network_reference.json and source receipts. CC BY 4.0 plus applicable original-provider terms: https://string-db.org/cgi/access?footer_active_subpage=licensing . Association edges do not prove direct physical interaction or causal regulation. Protein preferred-name to HGNC mapping is explicit; excluded ambiguities are retained. |

All downloads, processing, environments, temporary files and outputs for this round are under `SRV01:/home/analyst/te`. Public raw matrices and large generated data are ignored by Git; committed source URLs and retrieval/processing scripts provide the reproduction route. No new data-use restriction, scientific approval or legal conclusion is invented here.

## Denominators and non-equivalence

Discovery uses 22 people (12 ENDO, 10 non-ENDO), 150,788 inherited-label primary cells and 26,109 unambiguous named genes. Its original Alevin fractional estimates are preserved. HECA selection uses 22 people (15 ENDO, 7 non-ENDO) and 52,744 inherited-label primary nuclei; those counts are original integers from the selected public matrices. Two libraries from one person never become two independent people.

The assays and original gene reference/filtering differ. Gene-ID overlap, marker plausibility and source-supported recruitment independence are separately audited. They do not establish identical cell composition, pure anatomical origin, independent reannotation, equivalent gene models or absence of residual clinical confounding. Missing gene measurements are not evidence of no disease effect.

## Distinct freezes

1. Source/clinical/count preparation establishes eligibility and traceability.
2. The count/design contract snapshot establishes exact column identity, finite nonnegative values and model estimability. It cannot by itself pass scientific M0.
3. Expression filtering and the final discovery-induced graph universe are fixed during statistical implementation.
4. All rankings, selected lists, directions and evaluation rules are frozen before external disease effects.

No disease effects, rankings or external replication results are produced by this provenance document.
