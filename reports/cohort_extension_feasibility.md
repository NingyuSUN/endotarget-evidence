# Cohort extension feasibility — metadata only

Date: 2026-10-27. Scope: whether the two public single-cell endometriosis studies registered
in [DATA_SOURCES.md](../docs/DATA_SOURCES.md) as possible future resources, **GSE213216**
(Fonseca et al. 2023) and **GSE179640** (Tan et al. 2022), can extend this project's
donor-level stromal contrasts: the v1 proliferative-phase disease contrast, the secretory
program-level disease contrast, or a donor-level meta-analysis across cohorts. Script:
[`scripts/cohort_extension_metadata_audit.py`](../scripts/cohort_extension_metadata_audit.py);
receipt: [`cohort_extension_metadata_audit.json`](../receipts/cohort_extension_metadata_audit.json);
manifests: [GSE213216 samples](../manifests/cohort_extension_GSE213216_samples.tsv),
[GSE179640 samples](../manifests/cohort_extension_GSE179640_samples.tsv),
[eutopic donor ledger](../manifests/cohort_extension_eutopic_donor_ledger.tsv),
[lesion pairs](../manifests/cohort_extension_lesion_pairs.tsv).

**Result: neither study can extend the eutopic donor-level contrasts, and no expression
download is justified for that purpose.** Under the project's natural-cycle rule, GSE213216
retains 5 eutopic single-cell donors, all with endometriosis (3 proliferative, 2 secretory)
and no eligible control; GSE179640 retains 1 (a proliferative control). GSE179640 was
designed as a hormonally treated cohort. The planned multi-cohort meta-analysis is therefore
not possible with these accessions, and the v1 bound and the secretory-program result stand
as this project's statements for public single-cell eutopic data. Both studies are
lesion-context resources: GSE179640 offers 8 matched eutopic–lesion single-cell pairs under
hormonal treatment, GSE213216 4 natural-cycle pairs. Nothing was downloaded beyond public
metadata (SOFT family files, file listings and the NCBI series summary, 31 KB in total), and
no expression value was read.

## What was read

| Source | Content | Bytes |
|---|---|---|
| GEO SOFT family files, GSE213216 and GSE179640 | sample titles, source names, characteristics, library strategy, supplementary file names, series summary and design | 6,115 and 8,690 (gzip) |
| GEO supplementary file listings | per-file names, sizes and the RAW archive size | 4,354 and 5,314 |
| NCBI series summary | PubMed links, sample counts | 1 request |
| HECA authors' public donor and library metadata (already local from the v1 clinical gate; repository commit 244c282c) | per-donor cycle stage, hormonal treatment, age and endometrial pathology for the Fonseca and Tan donors, and per-library tissue, QC and donor assignment | 30,823 + 45,010 + 76,356 |

The GEO records themselves carry **no cycle phase, no hormonal treatment and no age** for
GSE213216 and **no phase or hormone fields** for GSE179640. Those fields come from the HECA
authors' curation of the two studies, which they integrated into the HECA cells atlas. This
audit uses that curation as published and does not re-derive it from the primary papers'
supplementary tables; any later use of either study would have to confirm the fields
against the primary sources first.

## Rules applied

The v1 rule family, metadata level only: eutopic endometrium, single-cell gene-expression
library (bulk tissue, organoid and hashing libraries excluded); natural cycle, meaning no
hormonal treatment and a curated stage of proliferative or secretory (menstrual, hormone-
treated and menopausal donors excluded); disease by the study's own assignment or the
curated endometrial pathology; age present. A stratum is usable for a disease contrast if
both disease groups reach **5 donors** within one phase (the HECA floor; 8 is the discovery
minimum). A paired lesion-versus-eutopic design needs **5 matched pairs**.

## GSE213216 (Fonseca et al. 2023; 21 patients, 49 scRNA libraries, 2 Visium)

| Eutopic sample | Patient | Curated donor | Disease | Stage | Hormonal treatment | Age | Verdict |
|---|---|---|---|---|---|---|---|
| GSM6634274 | 3 | Fonseca_3 | non-ENDO (fibroids, polyp) | menopause, hormones | Vivelle | 62 | excluded |
| GSM6574512 | 6 | Fonseca_6 | non-ENDO (adenomyosis, fibroids) | menopause, hormones | estradiol patch + progesterone | 58 | excluded |
| GSM6574528 | 10 | Fonseca_10 | ENDO, stage 4 | proliferative | none | 40 | **eligible** |
| GSM6574534 | 11 | Fonseca_11 | ENDO + adenomyosis, stage 2 | secretory | none | 42 | **eligible** |
| GSM6574535 | 12 | Fonseca_12 | non-ENDO (adenomyosis, fibroids) | hormones | norethindrone acetate | 48 | excluded |
| GSM6574536 | 13 | Fonseca_13 | ENDO + adenomyosis | hormones | estrogen/testosterone creams + progesterone | 38 | excluded |
| GSM6574539 | 14 | Fonseca_14 | ENDO, stage 2 | secretory | none | 28 | **eligible** |
| GSM6574542 | 16 | Fonseca_16 | ENDO, stage 2 | proliferative | none | 21 | **eligible** |
| GSM6574543 | 17 | Fonseca_17 | ENDO + adenomyosis + fibroids, stage 4 | proliferative | none | 34 | **eligible** |
| GSM6574544 | 18 | Fonseca_18 | ENDO + fibroids, stage 4 | perimenopause, hormones | progesterone and estrogen pellets | 52 | excluded |

Eligible natural-cycle eutopic donors: **proliferative 3 ENDO / 0 non-ENDO; secretory 2 ENDO
/ 0 non-ENDO.** The study's non-endometriosis endometrium comes from older, hormone-treated
or post-menopausal surgical patients; the "no endometriosis detected" peritoneal and bladder
biopsies belong to endometriosis patients and are not controls. Lesion material is rich: 25
peritoneal or deep lesions and 8 endometriomas across 16 patients, with matched eutopic
endometrium in 4 natural-cycle patients (10 and 16 proliferative, 11 and 14 secretory).
Public files: one `tar.gz` of Cell Ranger output per library, 50 MB to 1.1 GB each, 16.8 GB
in total; samples 2, 3 and 4 (GSM6634272–4) list no supplementary matrix. Stromal and
cell-type labels are not in the GEO record.

## GSE179640 (Tan et al. 2022; 27 patients; 31 tissue scRNA libraries, 24 bulk, 4 organoid)

| Eutopic sample | Donor | Disease (study) | Curated stage | Hormonal treatment | Age | Verdict |
|---|---|---|---|---|---|---|
| GSM6102532 | C01 | control | hormones | norethindrone/estradiol 1/20 | 36 | excluded |
| GSM6102533 | C02 | control | hormones (early proliferative) | norethindrone/estradiol 1/20 | 22 | excluded |
| GSM6102534 | C03 | control | proliferative | none | 42 | **eligible** |
| GSM6102537 | E01 | endometriosis, stage 3 | menstrual | none | 40 | excluded (menstrual) |
| GSM6102540 | E02 | endometriosis, stage 2 | hormones (early proliferative) | norethindrone/estradiol 1/20 | 40 | excluded |
| GSM6102543 | E03 | endometriosis, stage 4 | hormones (inactive) | norethindrone/estradiol 1/20 | 33 | excluded |
| GSM6102546 | E04 | endometriosis, stage 4 | hormones | norethindrone/estradiol 1/20 | 27 | excluded |
| GSM6102549 | E05 | endometriosis, stage 4 | hormones (early proliferative) | norethindrone/estradiol 1/20 | 33 | excluded |
| GSM6102551 | E06 | endometriosis, stage 4 | hormones | norethindrone/estradiol 1/10 | 30 | excluded |
| GSM6102554 | E07 | endometriosis, stage 4 | hormones (inactive) | norethindrone/estradiol 1/10 | 34 | excluded |
| GSM6102555 | E08 | endometriosis, stage 4 | hormones (inactive) | norethindrone/estradiol 1/20 | 40 | excluded |
| GSM6102560 | E09 | endometriosis, stage 4 | hormones (early secretory) | norethindrone/estradiol 1/20 | 45 | excluded |

Eligible natural-cycle eutopic donors: **1 proliferative control, 0 endometriosis.** The
primary paper describes its cohort as individuals undergoing hormonal treatment; the
curated fields agree (10 of 12 eutopic donors on norethindrone/estradiol-type regimens).
Matched eutopic–lesion single-cell pairs: **8** (E01–E07, E09; peritoneal lesions, with
adjacent peritoneum for six and ovarian lesions for two), 7 of them under treatment and E01
in the menstrual phase. Public files: one filtered Cell Ranger `h5` per library, 10 to 36
MB, 795 MB in total; the bulk libraries (5 control and 7 endometriosis eutopic, plus
lesions) are whole tissue and outside the lineage-level design, and their donors' hormone
status is not curated.

## Overlap with the project's cohorts

- The v1 and positive-control HECA cohorts use only the newly generated nuclei of
  E-MTAB-14039 (Oxford recruitment). The Fonseca (Los Angeles) and Tan (Connecticut) donors
  are not among them; the discovery cohort is a Bern recruitment. No donor overlap is
  expected by source; none of this is genotype-confirmed, as for v1.
- Both studies are external datasets **inside the HECA cells atlas** (48 Fonseca and 31 Tan
  libraries in the authors' library table; the Tan libraries are flagged as analysed, the
  Fonseca libraries as not). The HECA label vocabulary and marker statements were derived on
  that atlas. Any analysis of these studies with HECA-derived labels is label transfer
  within the same resource and would not count as independent of HECA.

## Decisions

1. **Close the extension item.** No stratum of either study reaches the floor for a
   eutopic disease contrast in either phase; a donor-level meta-analysis with these
   accessions cannot be run. No expression download.
2. **Register both as lesion-context resources** in the data-sources register, with the
   counts above. A paired lesion-versus-eutopic stromal contrast is numerically feasible only
   in GSE179640 (8 pairs) and only under hormonal treatment; it is a different question from
   this project's eutopic contrasts, would need its own prespecified protocol, its own
   lineage definition without HECA labels, and confirmation of the curated clinical fields
   against the primary publication, and it is not adopted here.
3. **What this means for the project's conclusions.** The v1 proliferative null with its
   minimum-detectable-effect bound, the phase positive control and the secretory program
   result are the project's statements for public single-cell eutopic endometrium; the
   remaining public eutopic single-cell donors with natural-cycle, phase-resolved metadata are
   already in those analyses.

## Sources

According to PubMed: Fonseca et al. 2023, *Nat Genet*
([DOI 10.1038/s41588-022-01254-1](https://doi.org/10.1038/s41588-022-01254-1); GSE213216);
Tan et al. 2022, *Nat Cell Biol*
([DOI 10.1038/s41556-022-00961-5](https://doi.org/10.1038/s41556-022-00961-5); GSE179640;
author correction [DOI 10.1038/s41556-022-01023-6](https://doi.org/10.1038/s41556-022-01023-6));
Marečková et al. 2024, *Nat Genet*
([DOI 10.1038/s41588-024-01873-w](https://doi.org/10.1038/s41588-024-01873-w); HECA
metadata curation).
