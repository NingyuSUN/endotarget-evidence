# Multi-Evidence Endometriosis Lesion Target Dossiers (EndoTarget v2)

Generated strictly under locked `configs/protocol_v2.yaml` and `receipts/prediction_seal.json`.
Data sources include GEO (`GSE179640`, `GSE213216`, `GSE141549`), Open Targets Platform (`MONDO_0005133` / `EFO_0001065`; https://platform-docs.opentargets.org/licence), Rahmioglu et al. 2023 (*Nat Genet*), and GTEx Portal API v10 (https://gtexportal.org/home/license).

- **Tier 1 (`Tier_1_Genetics_Anchored_MultiModal`)**: **1** targets backed simultaneously by human endometriosis germline genetics (`G_prior > 0`), discovery dual-contrast host-decontaminated remodeling (`GSE179640`), external scRNA stromal concordance (`GSE213216`, `p < 0.05`), and external host-peritoneum-controlled tissue validation (`GSE141549` paired `Peritoneal_Lesion vs PP` `BH q < 0.05`).
- **Tier 2 (`Tier_2_Replicated_Niche_Tractable`)**: **0** frozen `M1`/`M2` Top-50 targets achieving strict `BH q < 0.05` replication in **both** external scRNA (`GSE213216`) and external host-controlled bulk (`GSE141549` `PeL vs PP`) with active vascular/myeloid L–R crosstalk or small-molecule/antibody tractability.

| Tier | Gene | Lesion Dir | Rank (M2/M1/B0) | Disc logFC (EuE / AdjP) | GSE213216 scRNA logFC (q) | GSE141549 PeL vs PP logFC (q) | DIE vs PP logFC (q) | Genetic Prior (OT / GWAS-SMR) | Niche L-R Axis | Tractability (SM / AB) | GTEx TPM (Ut/Fib/Ov) |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| Tier 1 | **RUNX1** (`ENSG00000159216`) | `UP_in_Lesion` | 121 / 406 / 55 | +2.51 / +0.45 | +1.61 (8.24e-03) | +0.84 (4.92e-02) | +0.53 (6.66e-02) | 0.39 (0.65 / False) | `-` (-) | SM=True, AB=False | 1.36/20.22/4.86 |
