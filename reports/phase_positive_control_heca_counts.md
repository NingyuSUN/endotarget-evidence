# Phase positive control — HECA counts prepared (protocol step 4, executed early)

Date: 2026-10-26. Executed on the native analysis host (SRV01) in the project's
preserved Python environment (Python 3.10.16, NumPy 2.2.6). **No phase or disease model
was fitted; no expression effect exists.** Downloading and aggregating counts is permitted
before the F1 list freeze under the protocol; fitting the HECA model is not.
Script: [`scripts/phase_positive_control_prepare_heca_counts.py`](../scripts/phase_positive_control_prepare_heca_counts.py).
Receipts: [downloads](../receipts/phase_positive_control_heca_downloads.json),
[counts](../receipts/phase_positive_control_heca_counts.json).

## Download

Twenty public E-MTAB-14039 processed-matrix archives for the frozen secretory-arm samples,
1,378,729,112 bytes in total, each checked against the server-reported size during
transfer and hashed on arrival; ETag, Last-Modified and wall time are in the receipt.
Three concurrent streams, 732 s. Archives are stored one per directory and are not
committed.

## Aggregation and checks

The script reuses the v1 MatrixMarket parser, barcode join, feature check and checkpoint
validation unchanged, and writes only to `data/processed/phase_positive_control/heca/`.
The 22 v1 proliferative donors were re-derived from their frozen archives in the same run.

| Check | Result |
|---|---|
| Donors | 42 = 22 proliferative (15 ENDO / 7 non-ENDO) + 20 secretory (16 / 4) |
| Primary stromal nuclei | 94,150 = 52,744 (v1) + 41,406 (secretory); every sample equals its ledger count |
| v1 lineage positions in the H5 observation table | identical to the frozen v1 index file |
| v1 donor columns in the new stable matrix | **identical, value for value**, to the frozen v1 matrix |
| v1 gene rows and feature mapping | identical (26,125 stable IDs; 10,474 unmapped and 2 ambiguous excluded, as in v1) |
| Frozen v1 matrix hash before / after | `77121542…29af0` both times; nothing in the v1 tree was written |
| Source features across all 42 archives | identical (36,601 CellRanger features) |
| Count integrity | integer, finite, non-negative; exact mass conservation per archive and through mapping |
| Raw totals versus author H5 per nucleus | 79,880 of 94,150 exact; 14,270 higher in the raw matrix by at most 9 UMIs; none lower — the same feature-filtering pattern v1 recorded |
| Shared stable features with discovery | 21,585 |
| Runtime | 606 s; peak resident memory about 0.7 GiB |

Outputs: `source_counts.tsv.gz` (36,601 × 42), `counts_stable.tsv.gz` (26,125 × 42),
`selected_donor_count_ledger.tsv`, `source_features.tsv`, `heca_features_mapping.tsv`,
`shared_stable_features.tsv`, `raw_vs_author_H5_cell_totals.tsv.gz`, `provenance.json`
and per-sample checkpoints, all under the new tree. The output hashes are in the receipt.

## What this does and does not establish

The HECA secretory arm is now analysis-ready in the exact format the v1 model code
consumes, and the pipeline reproduces the frozen v1 columns exactly from raw archives.
That is an engineering identity check on the same public data. It says nothing about
phase biology, disease, or the positive control's endpoints, which await the F1 freeze
of the discovery lists. The discovery secretory arm is still not acquired; see the
[acquisition plan](../docs/PHASE_POSITIVE_CONTROL_DISCOVERY_COUNTS_PLAN.md).
