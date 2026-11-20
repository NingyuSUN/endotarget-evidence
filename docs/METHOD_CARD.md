# EndoTarget v2 — Method Card

**Purpose:** Auditable multi-cohort benchmark testing whether dual-contrast host-tissue decontamination (`M1`) and human-genetics + niche-crosstalk integration (`M2`) overcome host-background confounding (`B0`) and global PPI hub bias (`G1`) in ectopic endometriosis target prioritization.

**Headline verdict:**
1. **Host-controlled external tissue validation (`GSE141549`, 408 arrays / 162 women) confirms the core confounding hypothesis:** Naive paired lesion ranking (`B0`) achieves an apparent `13 / 20` (`65%`) `BH q < 0.05` replication when comparing peritoneal lesions to eutopic endometrium (`PeL vs PE`, $n=41$ pairs), yet **100% (`13 / 13`) of those hits fail when compared against matched unaffected peritoneum from the same women (`PeL vs PP`, $n=20$ pairs; `10 / 13` invert sign)**. By contrast, host-decontaminated `M1` and `M2` preserve positive transfer against host peritoneum (`Paired PeL vs PP` direction concordance `63.2%` and `72.2%` vs `42.1%` for `B0`; paired bootstrap `M2 − B0 = +31.4` percentage points `[+8.2, +46.5]`; `All PeL vs PP+CP` program `AUC = 0.651, p = 0.035` for `M1` and `AUC = 0.684, p = 0.0045` for `M2`, vs `AUC = 0.430, p = 0.115` for `B0`; `Paired DIE vs PP` concordance `88.9%` for `M2` vs `36.8%` for `B0`).
2. **Uncontrolled external single-cell `EcP vs EuE` (`GSE213216`, 29 libraries / 18 patients) exposes two real cross-cohort limitations of `M1`/`M2`:** On `GSE213216` `EcP vs EuE` ($n=12$ vs $9$), `B0` appears to replicate `9 / 20` (`45%`, `AUC = 0.991`) — driven by the exact peritoneal resident fibroblast markers (`SFRP2`, `GPC3`, `TBX18`, `CDON`, `PDLIM3`, `SKAP2`, `S100A10`) that invert against `PP` in `GSE141549` — whereas `M1` and `M2` replicate `0 / 20` at `BH q < 0.05` (`M1 − B0 = −42.7` pp `[−55.0, −30.0]`; `M2 − B0 = −43.9` pp `[−55.0, −30.0]`). Diagnosis identifies (a) **hormonal regimen mismatch** (`7 / 8` discovery `GSE179640` pairs were on exogenous `norethindrone/estradiol`, suppressing acute cytokines `CXCL8`, `IL1B`, `SPHK1`, `KCNK3` in uterine `EuE`, whereas `GSE213216` `EuE` donors had untreated cycling endometrium with high physiological expression of those genes) and (b) **small `AdjP` sample size ($n=6$ triplets, `0` individual `BH q < 0.05` genes) plus minor lymphoid doublet bleed in one 120-cell discovery cluster**.
3. **Human genetics and multi-evidence target nomination:** Unlike eutopic stroma in v1 (null enrichment), human endometriosis GWAS/L2G/SMR genes are significantly enriched in ectopic lesion remodeling across discovery perivascular mural (`p = 7.5e-4`) and endothelial (`p = 0.0045`) cells, external `GSE213216` lesion stroma (`p = 1.9e-4`), and host-controlled `GSE141549` `PeL vs PP` (`p = 2.9e-5`). Exactly **one target (`RUNX1`, `ENSG00000159216`)** satisfies all four pre-locked Tier-1 criteria simultaneously (`OT genetic_association = 0.653`, discovery `EcP vs EuE +2.51` and `EcP vs AdjP +0.45`, external scRNA `GSE213216 +1.61 [q = 0.0082]`, external host-controlled `GSE141549 PeL vs PP +0.84 [q = 0.0492]`, small-molecule tractable).

---

## 1. Input Provenance and Design Matrix Summary

| Layer | Source & Identifier | Scope Retained | Statistical Unit & Modeling |
|---|---|---|---|
| **Discovery scRNA** | `GSE179640` (Tan et al. 2022) | `31` tissue `.h5` libraries, `115,402` QC cells (`14` donors; `37,691` stromal fibroblasts) | Donor × site pseudobulk (`min_cells >= 20`); paired `~ donor + site` `limma-voom` (`8` `EcP–EuE` pairs, `df_res=7`; `6` `EcP–AdjP` triplets, `df_res=5`) |
| **External scRNA** | `GSE213216` (Fonseca et al. 2023) | `29` locked `.tar.gz` libraries, `242,430` QC cells (`18` patients; `106,424` stromal fibroblasts) | Donor × site pseudobulk; primary `EcP (12) vs EuE (9)` across `17` patients via `duplicateCorrelation` GLS (`rho = 0.159`, `df_res = 15`) + `4` sensitivities |
| **External Bulk** | `GSE141549` (Gabriel et al. 2020) | `408` BeadChip arrays (`28,247` probes $\rightarrow$ `19,262` unique `ENSG` genes), `162` subjects (`SL 001–162`) | Subject-deduplicated per tissue group; paired `PeL vs PE` (`n=41`, `df=40`), host-controlled `PeL vs PP` (`n=20`, `df=19`), `DIE vs PP` (`n=20`), triplet (`n=18`), and covariate-adjusted GLS |
| **Human Genetics & Druggability** | Open Targets v26.09 (`MONDO_0005133` / `EFO_0001065`); Rahmioglu et al. 2023 (*Nat Genet*); GTEx v10 | `3,513` Open Targets genes (`356` `genetic_association > 0`, `55 >= 0.5`); `404` combined genetically supported genes | Pre-locked bounded score $G_{\text{prior}} \in [0,1]$ excluding literature text-mining and RNA expression channels |
| **PPI & Niche L–R** | STRING v12.5 (`9606`, direct `experiments`/`database >= 700`); curated L–R panel | `159,312` reference edges (`97,462` induced on `10,972`-gene shared universe $\mathcal{U}$); `168` L–R pairs | Fixed `r = 0.5` RWR for `G1`; discovery endothelial/myeloid partner activity $N_{\text{LR}} \in [0,1]$ for `M2` |

---

## 2. Why Two Apparent Contradictions Are Actually One Coherent Geometric Fact

| Observation | Naive Interpretation (Wrong) | Mechanistic Reality Demonstrated by Dual-Contrast Audit |
|---|---|---|
| `B0` replicates `9/20` (`45%`, `AUC 0.99`) in `GSE213216` scRNA `EcP vs EuE` and `13/20` (`65%`) in `GSE141549` `PeL vs PE`, but `1/20` (`5%`, `0/20` joint) in `GSE141549` `PeL vs PP` | "Naive DE transfers well across cohorts" | Comparing pelvic peritoneal biopsies to uterine biopsies reproducibly contrasts **normal peritoneal serosal mesenchyme (`SFRP2`, `GPC3`, `TBX18`, `CDON`, `PDLIM3`, `S100A10`)** against Mullerian endometrium. Because `GSE213216` `EcP vs EuE` repeats the same tissue-site contrast without an unaffected peritoneum control, it reproduces the **site-of-biopsy confound**, not lesion pathology. |
| `M1`/`M2` outperform `B0` and `G1` across every host-peritoneum-controlled contrast in `GSE141549` (`PeL vs PP`, `PeL vs PP+CP`, `DIE vs PP`, `OMA vs PP+CP`), yet score `0/20` in `GSE213216` `EcP vs EuE` | "Dual-contrast filtering overfit discovery" | `M1`/`M2` successfully removed peritoneal background genes, leaving inflammatory/hypoxic remodeling genes (`CXCL8`, `IL1B`, `CXCR4`, `KCNK3`, `SPHK1`). Those genes ARE elevated over host peritoneum `PP` in `GSE141549`, but in `GSE213216` `EcP vs EuE` their sign flips because `GSE213216` `EuE` is **untreated natural-cycle endometrium** (high baseline cytokines) whereas `GSE179640` `EuE` was **norethindrone/estradiol-suppressed**. |

---

## 3. Scope, Limitations, and Licensing Notices

1. **Observational reanalysis, not clinical validation**: All cohorts are retrospective surgical samples. Differential expression and genetic colocalization/L2G prioritize candidate mechanisms (`RUNX1`, `CXCR4–CXCL12`, `ESR1` loss, vascular/mural remodeling) but do not establish clinical efficacy, safety, or therapeutic modality direction (agonist vs inhibitor) without functional perturbation.
2. **Hormonal treatment heterogeneity**: Public lesion scRNA cohorts differ systematically in hormone use (`GSE179640` predominantly norethindrone/estradiol-treated; `GSE213216` mixed natural-cycle eutopic and surgical lesions; `GSE141549` `177 / 408` arrays on hormonal medication). Future single-cell lesion benchmarks require hormone-matched `(EuE, EcP, AdjP)` triplets from `>= 15–20` patients.
3. **Whole-tissue vs single-cell estimands**: `GSE141549` validates host-peritoneum double-differences at whole-biopsy tissue level across `162` women; cellular composition differences between fibrotic lesion cores and thin normal serosa still contribute alongside cell-intrinsic transcriptional remodeling.
4. **Third-party database notices**:
   - **Open Targets Platform**: Queried via the official GraphQL wrapper; users must review terms at https://platform-docs.opentargets.org/licence (recorded in `.licenses/opentargets_database_LICENSE.txt`).
   - **GTEx Portal API v10**: Queried via the official GTEx CLI wrapper; users must review terms at https://gtexportal.org/home/license (recorded in `.licenses/gtex_database_LICENSE.txt`).
   - **STRING v12.5**: Used under CC BY 4.0 (https://string-db.org/cgi/access).
