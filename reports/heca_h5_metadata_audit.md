# HECA author H5AD obs identity audit

Status: VERIFIED author-file correspondence; final clinical eligibility, genotype truth, independent cohort status and inference readiness NOT VERIFIED.

This audit answers whether deposited cell/nucleus annotations connect to the public author donor records, and how that connection compares with the SDRF. It reads only obs metadata, with no expression dataset opened, no differential expression, no target ranking, and no genotype reanalysis.

## Verified source and bounded retrieval

Official anonymous source: https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_nuclei.h5ad

Exact size 8,739,186,814 bytes, checked against every HTTP 206 Content-Range. ETag `"3b257ce77e9a967b574b55810dfa189d-556"`; Last-Modified Wed, 08 Nov 2023 10:30:00 GMT. The multipart ETag is not a whole-file SHA256. Source byte ranges are hashed individually, with a full local range manifest and original response logs retained in ignored data/metadata/heca_h5/. Tracked receipts/heca_h5_source_validation.json gives the manifest hash and source validators. No whole-file hash is claimed.

Fetched 28,901,376 bytes (27.56 MiB; approximately 0.33% of the full object), within the 150 MiB ceiling. Final retrieval finished 370.393 seconds after the first source-cache write, within the 600-second cumulative exploration ceiling. Offline cross-table processing ran afterward without source requests. Blocks are aligned to 64 KiB; edge blocks may contain neighboring bytes, but expression objects were never opened or interpreted.

Runtime is isolated at envs/m0-h5, Python 3.10.16, numpy 2.2.6, h5py 3.14.0, requests 2.32.5. Exact installed dependencies are requirements/metadata-probe.txt. System Python lacked ensurepip; the newly created incomplete venv was rebuilt with the existing Miniconda Python, leaving existing environments unchanged. One source connection dropped; a bounded two-retry adapter and hashed cache recovered the export.

## What the original obs establishes

312,246 author-retained nuclei, 68 combined library/donor samples, 41 raw libraries, 63 author genotype IDs. All 312,246 library-plus-barcode keys are unique; every key has its decoded library prefix and a 16-base A/C/G/T barcode. Every sample equals Library_genotype. Raw barcodes alone are not asserted globally unique.

This older AnnData file stores categorical codes as integer datasets plus shared obs/__categories dictionaries. The probe decodes that original scheme; integer codes must not be mistaken for donor identifiers or disease labels. Complete original fields, category dictionaries and bounded heads are in reports/heca_h5_schema.json. Cell-cycle `phase` has G1/G2M/S and is distinct from menstrual-cycle `Binary Stage` and `Stage`.

Original obs includes sample, library, genotype, Library_genotype, Group, Age, Hormonal treatment, Binary Stage, Stage, Biopsy_type, Tissue_sampled, author QC labels, LT_celltype, LT_mesen_celltype, celltype, lineage and label_long. It has no explicit biopsy date, hormone washout duration or donor-specific independent phase-ascertainment provenance.

## Identity and clinical label crosswalk

Join performed: unique library-plus-barcode -> H5 sample/Library_genotype/genotype -> author metadata-libraries_deconvoluted.csv Library_genotype/Donor -> metadata-donors.csv Donor. Author CSVs were retrieved by the separate clinical branch at pinned repository commit 244c282c20c3efd53707e3eb9f58bbd1cbd56ad7; their source receipts are receipts/discovery_clinical_gate_heca_downloads.json. Their file hashes are also recorded in receipts/heca_h5_audit.json.

| Comparison | Matches / 68 samples | Meaning |
|---|---:|---|
| Unique author deconvolution record | 68/68 | All observed combined sample IDs have one public author record. |
| H5 genotype versus author Donor | 68/68 | Published cell assignment agrees with published donor crosswalk. |
| H5 disease Group versus author Group | 68/68 | Same author annotation correspondence. |
| H5 exact Age versus author Age | 68/68 | Actual age values can replace SDRF recruitment-range placeholders. |
| H5 genotype versus SDRF declared individual | 40/68 | 28 deposited identity annotations conflict. |
| H5 disease versus SDRF disease | 52/68 | 16 disease labels differ. |
| H5 broad cycle phase versus SDRF broad phase | 54/68 | 14 phase labels differ. |
| H5 broad cycle phase versus author Binary Stage | 68/68 | Author final clinical category agrees with H5, distinct from histology-only evidence. |

Among the original 22 declared SDRF proliferative candidate IDs, 12 are involved in the identity disagreement. The original source-only report remains a correct description of SDRF, but its 16-case/6-control candidate count must not be used to join cells. Actual H5 author proliferative annotation has 26 samples from 24 genotype IDs: 16 Endo_Superficial and 8 Ctrl_Superficial. No final eligibility is claimed.

Example: UA_Endo11282193_EX0540 has author genotype/Donor EX0540 and Endo_Superficial, whereas SDRF declares FX1146/normal. The paired sample ending FX1146 has author donor FX1146/normal. Deconvolution field `genotyped` is FX1146 for both rows, so it is not a unique donor assignment key. Similar behavior occurs for the FX1104/FX9005 pair. Use the explicit author Donor/genotype crosswalk, not the `genotyped` field or filename suffix alone.

PROPOSED precedence before effects: use H5 assignment supported by the matching author deconvolution and donor tables for the cell join, retaining all SDRF originals and discordances. This establishes documented author provenance, not independent genotype validation. A root clinical adjudication should record that rule before freezing any cohort. Do not silently rewrite SDRF.

## Hormones, ages and phase disagreement

Across 68 sample records, 54 carry the literal H5 category `nan` for Hormonal treatment and 14 carry explicit treatment labels. All 26 H5 proliferative sample records have literal `nan` and corresponding author donor CSV `NA`; zero have explicit-treatment/proliferative conflict. Missing treatment coding is not an explicit negative exposure record. Paper-defined natural-cycle eligibility (no hormone treatment for at least three months) may provide linked study-definition evidence if the author final category is documented to encode that group; it must remain separate from a measured negative hormone field. That interpretation is pending root adjudication.

The 24 author proliferative IDs have ages including EX0795=18 and FX0052=46. These differ from the universal SDRF “20 to 45” recruitment range; never treat that range as individual age or reject the actual author ages merely for disagreeing with the range.

H5 final proliferative labels do not uniformly agree with histology. Nine unique genotype IDs have a histology field starting with Secretory: EX0540, FX0011, FX1097, FX1215 (Secretory?), FX1236, FX1104, EX0795, FX1119, EX0288. FX0052 is mixed “Proliferative disordered; Secretory”. FX1101 and FX9024 have no histology phase. Corresponding expression-derived author phase fields are preserved, including Proliferative Late and Early Proliferative; final labels must not be represented as uniformly independent histology evidence. This audit does not generate phase labels from expression.

Repeated libraries for FX1101 and FX1236 remain explicit. Shared donor ID does not establish the same biopsy or a technical replicate; summing for inference awaits the protocol repeat rule.

## Primary stromal counts and annotation decision

PROPOSED annotation-only main count uses eStromal, eStromal_MMPs, eStromal_cycling, dStromal_early, dStromal_mid, dStromal_late and Fibroblast_basalis. It excludes perivascular ePV populations, uSMCs, epithelial, endothelial and immune cells. HOXA13 is kept separately and excluded from endometrial primary stroma: the original HECA article describes mesenchymal HOXA13+ cells as likely cervical contamination, supported by comparison with a public cervix dataset (Extended Data Fig. 1e–i). The same paragraph describes uSMCs/myometrial PV cells as myometrial populations. Source: [Mareckova et al., original article](https://doi.org/10.1038/s41588-024-01873-w), local data/reference/heca_article.xml, paragraph beginning 'As expected, most of the cell populations'. This is source-supported anatomy, not a choice made to meet donor minima. This mapping is not frozen or marker-validated, and no expression markers were inspected.

Across all 312,246 nuclei: 122,229 match the proposed primary labels; 5,133 are HOXA13; all Mesenchymal contains 129,370. Thus all Mesenchymal is not the same biological object as the proposed primary stroma.

Under that proposed mapping, among 24 author proliferative IDs, FX1215 has 0 primary nuclei and FX1141 has 45. Other candidates have at least 203; the descriptive >=50 count is therefore 15 case and 7 control IDs before clinical, repeated-biopsy and design adjudication. These are observed integer annotation counts, not final qualifying donor totals.

The HOXA13 choice is material: FX1215 has 57 HOXA13 nuclei and FX1141 has 16. Including them would move both above 50. Do not choose annotation boundaries to recover minimum sample counts. The author anatomical evidence supports excluding HOXA13 regardless of whether inclusion would recover donor minima. The remaining label crosswalk still needs a formal freeze before effects.

## Validation, outputs and open gates

The probe verifies partial-response size, range, stable response validators, category-code bounds, equal obs lengths, complete crosswalk coverage, unique library/barcode keys and per-sample clinical consistency. Negative checks reject unbounded reads, reads above 8 MiB, out-of-object seeks, budget overflow and wall-time overflow before network transfer; exclusion checks keep mural/perivascular labels outside the proposed primary mapping. These are engineering/source-consistency checks, not external biological validation.

- reports/heca_h5_sample_crosswalk.tsv: all 68 original sample mappings, identity conflicts, clinical fields and proposed lineage counts.
- reports/heca_h5_donor_counts.tsv: all 63 genotype IDs with observed counts; repeated samples are visible and inferential aggregation is not authorized by this table.
- reports/heca_h5_author_label_counts.tsv: per-sample counts for every original author cell type.
- receipts/heca_h5_audit.json, receipts/heca_h5_ranges.json and receipts/heca_h5_source_validation.json: observed outcome, hashes, validators and range budget.
- Ignored data/metadata/heca_h5/obs_selected.tsv: original decoded per-nucleus metadata; 79 MiB local output, not a full expression download.

Commands executed on SRV01:
```
envs/m0-h5/bin/python scripts/heca_h5_metadata.py
envs/m0-h5/bin/python scripts/heca_h5_metadata.py --export sample library genotype Library_genotype Group Age "Hormonal treatment" "Binary Stage" Stage Biopsy_type Tissue_sampled LT_celltype LT_mesen_celltype celltype lineage label_long Endometrial_pathology Endometriosis_stage discard_lowQC _index
envs/m0-h5/bin/python scripts/heca_h5_metadata.py --audit
envs/m0-h5/bin/python scripts/heca_h5_metadata.py --self-test
```

Remaining: root adjudication of author-vs-SDRF precedence, natural-cycle/hormone evidence, histology/expression phase differences, the remaining anatomical label crosswalk freeze (HOXA13 exclusion now supported by the primary source), repeated biopsies, independence, count-layer validity and final batch/design checks. No inference runtime, matrix download, DE or disease-ranking result was produced. No user decision is required for this read-only branch; the parent task owns the already-authorized evidence adjudication.

YELLOW — understanding not yet assessed: independently explain donor versus library, why an authoritative author crosswalk is correspondence rather than biological identity revalidation, why missing hormones are not measured absence, why cycle phase may be expression-derived, and why cervical-origin HOXA13 inclusion cannot be chosen to rescue a threshold.

