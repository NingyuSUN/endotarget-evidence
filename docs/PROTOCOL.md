# EndoTarget v2 — Prespecified Protocol (`configs/protocol_v2.yaml`)

## 1. Motivation and Continuity from EndoTarget v1

In EndoTarget v1 (`endotarget-evidence`), cross-sectional donor-level comparisons of proliferative-phase **eutopic** endometrial stromal fibroblasts (`GSE266265`, 22 donors; `E-MTAB-14039`, 22 donors; pooled 44 donors) detected `0` differentially expressed genes at `BH q < 0.05`, bounded median stromal disease differences below `1.5–1.8`-fold, found no enrichment of endometriosis GWAS genes in eutopic stroma (`Wilcoxon p = 0.42–0.98`), and showed that global STRING random walk with restart (`G1`) degraded cross-study sign agreement relative to raw expression ranking (`B0`) by `−38.4` percentage points (`−37` points even on a 7,367-gene cycle-phase positive control).

EndoTarget v2 moves from cross-sectional eutopic endometrium to **ectopic endometriosis lesions**, testing two explicit methodological hypotheses under a cryptographically sealed 4-commit zero-leakage chain:

1. **Host-tissue confounding of naive lesion contrasts (`B0` vs `M1`)**: Comparing peritoneal lesions (`EcP`) against eutopic endometrium (`EuE`) confounds lesion-intrinsic pathological remodeling with normal resident peritoneal fibroblast identity and ambient serosal mesothelial RNA (`MSLN`, `UPK3B`, `PRG4`, `S100A10`). Requiring concordant divergence from **both** matched uterine origin (`EcP vs EuE`) and matched unaffected adjacent peritoneum (`EcP vs AdjP`), penalized for mesothelial enrichment (`M1`), isolates host-decontaminated lesion remodeling.
2. **Causal target anchoring via niche crosstalk and human genetics (`M2` vs `G1` / `B0`)**: Rather than diffusing discovery statistics globally across high-degree housekeeping PPI hubs (`G1`), `M2` gates strictly on positive dual-contrast lesion remodeling (`M1_core > 0`) and multiplicatively integrates active vascular/myeloid-to-stromal ligand–receptor (`L–R`) partner expression with independent human germline endometriosis genetics.

---

## 2. Zero-Leakage Four-Commit Execution Chain

| Stage / Commit | Scope | Cryptographic Receipt | External Expression State |
|---|---|---|---|
| **Stage 1 (`Commit 1`, `33798ca`)** | Metadata-only audit of `GSE179640`, `GSE213216`, `GSE141549`; Open Targets + Rahmioglu 2023 human genetics table; de novo 7-lineage canonical markers; 168 curated L–R pairs; STRING v12.5 direct gene graph; freeze `configs/protocol_v2.yaml` | `receipts/m0_inputs_and_priors.json`, `receipts/protocol_freeze.json` (`d05e7a20...`) | Not downloaded; `96,810` platform lines in `GSE141549` SOFT skipped |
| **Stage 2 (`Commit 2`, `3635800`)** | Download 31 `GSE179640` `.h5` files; de novo Leiden + canonical marker annotation; exact donor × site × lineage pseudobulk; paired `EcP vs EuE` (`n=8`) and `EcP vs AdjP` (`n=6`) `limma-voom` models; freeze `B0`, `G1`, `M1`, `M2` rankings and `K=20, 50` shortlists | `receipts/discovery_GSE179640_data.json`, `receipts/prediction_seal.json` (`37d98e8c...`) | Not downloaded; verified absent before seal |
| **Stage 3 (`Commit 3`, `9f021d8`)** | Download 29 locked `GSE213216` scRNA archives and `GSE141549` 408-array workbook; verify `protocol_freeze.json` and `prediction_seal.json` hashes; run single-pass locked external evaluation, 2,000-draw donor/subject bootstraps, genetic enrichment, target dossiers, and SVG/PNG figures | `receipts/external_cohorts_data.json`, `receipts/locked_evaluation_v2.json`, `receipts/figures_v2.json` | Read once under locked Stage 1 & Stage 2 hashes |
| **Stage 4 (`Commit 4`)** | Bilingual documentation (`README`, `RESULTS`, `PROTOCOL`, `METHOD_CARD`) and automated `pytest` contract verification | `tests/test_endotarget_v2.py` | Verified byte-identical to Stage 3 |

---

## 3. Cohorts, Eligibility, and De Novo Lineage Pseudobulk Rules

### 3.1 Discovery Single-Cell Cohort (`GSE179640`, Tan et al. *Nat Cell Biol* 2022)
- **Included libraries**: All `31` tissue single-cell gene-expression `.h5` libraries across `14` donors (`E01–E11`, `C01–C03`); `4` organoid/hashing libraries and `24` whole-tissue bulk libraries excluded at M0.
- **Within-patient pairing**:
  - `8` matched `(EuE, EcP)` pairs (`E01, E02, E03, E04, E05, E06, E07, E09`).
  - `6` matched `(EuE, EcP, AdjP)` triplets (`E01, E02, E03, E04, E05, E09`; `E09`'s two adjacent peritoneum libraries `EcPA` and `EcPA2` are summed at raw integer count level prior to normalization).

### 3.2 External Single-Cell Replication Cohort (`GSE213216`, Fonseca et al. *Nat Genet* 2023)
- **Deterministic M0 selection rule**: Exclude `2` Visium spatial libraries, `3` GSMs without public matrices on GEO (`GSM6634272–4`), `4` normal whole-ovary libraries, and `1` contradictory metadata record (`GSM6574531`, labeled `No endometriosis detected` with anatomic location `Bladder endometriosis`). Select exactly `1` library per `(patient_id, site_class)` by ascending GSM accession to prevent multi-biopsy pseudoreplication and avoid arbitrary library selection.
- **Locked selection**: `29` scRNA libraries across `18` unique patients (`12` `EcP` peritoneal/pelvic lesion donors, `9` `EuE` eutopic donors [`7` endometriosis, `2` surgical non-endometriosis], `6` `EcO` endometrioma donors, `2` unaffected peritoneum donors; `4` matched natural-cycle `(EuE, EcP)` pairs `P10, P11, P14, P16`).

### 3.3 External Multi-Subtype & Host-Controlled Bulk Cohort (`GSE141549`, Gabriel et al. 2020)
- **Complete subject linkage**: All `408` Illumina BeadChip arrays (`GPL10558`, `GPL13376`) map with zero unlinked samples to `162` unique subjects (`SL 001`–`SL 162`) via `GSE141549_Sample_link.xlsx`.
- **Subject-level deduplication**: Collapse technical replicates and multiple biopsies of the same tissue group within each subject by arithmetic mean on the batch-corrected $\log_2$ scale so every subject contributes at most `1` observation per group (`CE=38`, `CP=24`, `PE=92`, `PP=37`, `Peritoneal_Lesion=54`, `DIE=57`, `OMA=25`; paired within-subject `(PE, PeL)` `n=41`, paired `(PP, PeL)` `n=20`, triplet `(PE, PP, PeL)` `n=18`, paired `(PP, DIE)` `n=20`).

### 3.4 Single-Cell QC and De Novo Lineage Assignment (Zero HECA Labels)
Because both `GSE179640` and `GSE213216` were integrated into the HECA atlas, inheriting HECA-derived annotations would constitute cross-cohort label leakage. Both cohorts are annotated **de novo** from `56` pre-locked canonical markers across `7` lineages (`manifests/canonical_lineage_markers.tsv`: `stromal_fibroblast`, `perivascular_mural`, `epithelial`, `endothelial`, `myeloid_macrophage`, `lymphoid`, `mesothelial`):
1. **Cell QC**: `n_genes_by_counts >= 250`, `total_counts >= 500`, `pct_counts_mt <= 25.0%`.
2. **Clustering and lineage scoring**: Normalize to `CP10K`, `log1p`, select `2,000` Seurat-flavor HVGs, compute `30` PCs (`seed = 20261120`), `15`-NN graph, and `igraph` Leiden clustering (`resolution = 1.2`). Z-score each canonical marker across cells (clipped to `[-3, +6]`), assign each cluster to its highest mean lineage score (`top_z >= 0.15`, else `unassigned_low_score`), and flag individual cells whose top single-cell lineage score exceeds their cluster-assigned lineage score by `> 0.50` as `ambiguous_doublet`.
3. **Pseudobulk**: Sum exact raw integer UMI counts per `(donor_id, site_code, lineage)` bucket; retain pseudobulks with `>= 20` cells; assert exact integer UMI conservation across retained, `<20`-cell dropped, and excluded doublet/unassigned buckets.

---

## 4. Shared Gene Universe and Four Frozen Ranking Methods

All four methods rank the identical shared universe $\mathcal{U}$ (`10,972` genes): canonical unversioned `^ENSG[0-9]{11}$` genes reaching `CPM >= 1.0` in at least `25%` of paired `GSE179640` `(EuE, EcP)` stromal pseudobulk samples (`12,188` expressed genes) and having degree $\ge 1$ in the induced human STRING v12.5 direct `experiments`/`database >= 700` functional association subgraph (`97,462` induced undirected edges).

For each paired discovery contrast (`EcP vs EuE`, $n=8$ pairs; `EcP vs AdjP`, $n=6$ pairs), fit `~ donor + site` via `limma-voom` (TMM factors, tricube local linear mean–variance precision weights, empirical Bayes `squeeze_var` moderation) and convert moderated $t_g$ with total degrees of freedom $df_{\text{res}} + d_0$ to exact two-sided signed normal $z$-equivalents $z_g = \operatorname{sign}(t_g)\Phi^{-1}(1 - p_g/2)$. Let $s_g = \operatorname{sign}(\text{logFC}_g^{(\text{EcP}-\text{EuE})}) \in \{-1, +1\}$ denote the frozen discovery sign.

1. **`B0` (Naive Paired Lesion DE)**:
   $$B_0(g) = \left|z_g^{(\text{EcP}-\text{EuE})}\right|$$
2. **`G1` (Global STRING RWR — v1 Baseline)**:
   $$p^{(t+1)} = 0.5\,p^{(0)} + 0.5\,W p^{(t)}, \quad p^{(0)}_g = \frac{|t_g^{(\text{EcP}-\text{EuE})}|}{\sum_{u \in \mathcal{U}} |t_u^{(\text{EcP}-\text{EuE})}|}$$
   iterated on symmetric column-normalized STRING adjacency $W$ until $\|p^{(t+1)} - p^{(t)}\|_1 < 10^{-10}$.
3. **`M1` (Host-Decontaminated Dual-Contrast Lesion Score)**:
   $$M_{1,\text{core}}(g) = w_{\text{ambient}}(g)\,\min\!\Big(\left|z_g^{(\text{EcP}-\text{EuE})}\right|,\;\max\!\big(0,\,s_g \cdot z_g^{(\text{EcP}-\text{AdjP})}\big)\Big)$$
   where $w_{\text{ambient}}(g) = \big[1 + \exp\big((\overline{\log_2\text{CPM}}_{g,\text{meso}} - \overline{\log_2\text{CPM}}_{g,\text{EcP\_stroma}}) - 2.0\big)\big]^{-1}$, with deterministic tie-breaker $+10^{-4}\big(|z_g^{(\text{EcP}-\text{EuE})}| + s_g z_g^{(\text{EcP}-\text{AdjP})}\big)$.
4. **`M2` (Genetics & Niche Crosstalk-Anchored Target Score)**:
   $$M_2(g) = M_{1,\text{core}}(g)\,\big(1.0 + 1.0\,G_{\text{prior}}(g) + 0.5\,N_{\text{LR}}(g)\big) + 10^{-4}\big(|z_g^{(\text{EcP}-\text{EuE})}| + s_g z_g^{(\text{EcP}-\text{AdjP})}\big)$$
   where $G_{\text{prior}}(g) \in [0, 1]$ combines Open Targets (`MONDO_0005133` / `EFO_0001065`) `genetic_association` and Rahmioglu 2023 GWAS/coding/SMR indicators (locked at Commit 1), and $N_{\text{LR}}(g) \in [0, 1]$ measures active cognate partner expression ($\overline{\log_2(\text{CPM}+1)} \ge 2.5$) and lesion regulation across discovery `endothelial` and `myeloid_macrophage` pseudobulks over `168` curated L–R pairs.

---

## 5. Prespecified External Endpoints and Uncertainty Rules

- **External Single-Cell (`GSE213216`)**: Primary donor-blocked `limma::duplicateCorrelation` GLS model on TMM `log2CPM` comparing `EcP` ($n=12$) vs `EuE` ($n=9$) across `17` unique patients (`df_res = 15`), plus `4` sensitivities (`ENDO-only EuE`, `4` matched natural-cycle pairs, `8 vs 5` donor-disjoint unpaired, and `EcO vs EuE` endometrioma). Uncertainty via `2,000` patient-cluster bootstrap draws (`seed = 20261120`).
- **External Host-Controlled Bulk (`GSE141549`)**: Co-primary within-subject paired `Peritoneal_Lesion vs PE` ($n=41$ pairs) and **host-peritoneum-controlled `Peritoneal_Lesion vs PP`** ($n=20$ pairs), quantifying (i) `BH q < 0.05` concordant replication, (ii) joint dual-contrast replication (`PE` and `PP` both `BH q < 0.05` same sign), (iii) host false-positive fraction among `PE` hits, (iv) `All Peritoneal_Lesion vs PP+CP` ($n=54$ vs $61$) program projection $\theta$ and `AUC`, and (v) subtype transfer to `DIE` and `OMA` plus `PE vs CE` eutopic control, with `2,000` subject-cluster bootstrap draws.
