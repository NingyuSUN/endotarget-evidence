# Phase positive control — discovery counts prepared (protocol step 2)

Date: 2026-10-26. Executed on the native analysis host (SRV01) per the
[acquisition plan](../docs/PHASE_POSITIVE_CONTROL_DISCOVERY_COUNTS_PLAN.md). **No phase or
disease model was fitted; no expression effect exists.** Receipts:
[fixture check](../receipts/phase_positive_control_rds_fixture.json),
[discovery counts](../receipts/phase_positive_control_discovery_counts.json) (which embeds
the download receipt). Scripts: [parser port](../scripts/phase_positive_control_rds_stream.py),
[fixture](../scripts/phase_positive_control_rds_fixture.R) and
[fixture check](../scripts/phase_positive_control_rds_fixture_check.py),
[downloader](../scripts/phase_positive_control_geo_download.sh),
[aggregation](../scripts/phase_positive_control_prepare_discovery_counts.py).

## Pre-flight and fixture gate

The ported parser keeps the archived author parser's logic byte for byte and changes only
the resident-memory unit, the memory ceiling (24 GiB) and the disk reserve (75 GiB,
configurable). A base-R Seurat-like fixture (dgCMatrix slot layout, two assays, factor,
logical-with-NA and numeric metadata, character row names, gzip, serialization v3)
round-tripped exactly on all 29 checks, including the retained-only-counts policy and a
forced memory-ceiling abort. The project R lacks the Matrix package, so the fixture uses
plain S4 classes rather than the authors' original fixture script. The first fixture
attempt, on the root volume, was stopped by the 75 GiB reserve guard — the intended
behavior — and was repeated on the secondary volume.

## Download, extraction, aggregation

| Step | Measured |
|---|---|
| Download | `GSE266265_EndoAtlas.rds`, 41,729,517,596 bytes; four ranged resumable streams, 2,329 s (17.9 MB/s aggregate); concatenation 662 s; SHA-256 604 s; **hash equals the author-route record** `44a832cd…cf94f` |
| Extraction | gzip, serialization v3; 92,176,487,074 uncompressed bytes streamed in 1,218 s; 12,327,363 nodes; peak resident memory 1.83 GiB; 14 GB retained (counts `i`, `p`, `x` and `meta.data` numeric columns); everything else discarded |
| Structure | 34,727 features × 466,371 cells, 60 samples, 1,188,298,946 nonzero entries, 104,078,188 fractional, no negative or non-finite value, maximum 3,919 — all equal to the author record |
| Conservation | per-cell totals equal metadata `nCount_RNA` for **466,371 of 466,371 cells** (maximum absolute difference 0) |
| Aggregation | one sequential pass, 257 s, peak 1.9 GiB; 35 donors; 225,207 primary-lineage cells |

## Donors and lineage

Lineage rule as frozen in [`config/discovery_lineage.json`](../config/discovery_lineage.json):
`AnnotationMain` = stromal **and** `AnnotationRefined` in the 12 frozen labels; missing or
conflicting labels are excluded, never inferred. Exactly one cell in the 35 donors had a
lineage refined label with a missing main label (sample136, `eF1`); it is excluded and
counted in the receipt. The 22 v1 donors had no such cell, and their membership equals the
frozen v1 ledger exactly. For every donor the RDS metadata agreed with the frozen clinical
ledger on disease, cycle-phase label and capture batch.

| Secretory donor | Sub-phase | Disease | Primary-lineage cells |
|---|---|---|---:|
| sample75 | mid | ENDO | 5,732 |
| sample85 | mid | ENDO | 4,398 |
| sample88 | mid | non-ENDO | 5,439 |
| sample100 | mid | ENDO | 6,433 |
| sample104 | late | non-ENDO | 5,223 |
| sample113 | early | non-ENDO | 6,370 |
| sample117 | mid | non-ENDO | 5,860 |
| sample127 | early | ENDO | 4,241 |
| sample128 | mid | ENDO | 5,904 |
| sample131 | late | non-ENDO | 6,054 |
| sample132 | early | ENDO | 6,063 |
| sample133 | early | non-ENDO | 5,269 |
| sample136 | mid | ENDO | 7,433 |

All 13 secretory donors exceed the 50-cell minimum (74,419 cells in total); the
proliferative arm is the frozen v1 set (150,788 cells). The final design check of the
protocol is unchanged from the feasibility step: 13 versus 22 donors, primary design
full rank, residual df 22.

## Features and the v1 identity check

| Check | Result |
|---|---|
| Irreversible multi-source groups excluded | 50 groups, 177 features; 34,550 singletons, as in v1 |
| Pinned HGNC mapping | 26,109 stable IDs; decisions identical to the v1 mapping table, feature by feature |
| v1 donor columns in the new stable matrix | all 22 within **4.4 × 10⁻¹⁵ relative** of the frozen v1 matrix (float64 summation order; the protocol tolerance was 10⁻⁶) |
| v1 all-feature donor totals | agree with the v1 ledger |
| Frozen v1 matrix hash before / after | `c286af7c…e067` both times; nothing in the v1 tree was written |

Outputs under `data/processed/phase_positive_control/discovery/`:
`source_counts_all_features.tsv.gz` (34,727 × 35), `counts_stable.tsv.gz` (26,109 × 35),
`selected_donor_count_ledger.tsv`, `discovery_features_mapping.tsv`,
`irreversible_group_exclusions.tsv`, `selected_cell_totals_vs_nCount_RNA.tsv.gz`,
`provenance.json`. Fractional Alevin estimates are kept unrounded; nothing is normalized.

## Deviations from the plan, all operational

- The downloader's free-space pre-check compared a number printed in scientific notation
  and errored instead of deciding; the transfer proceeded on a volume with 1.3 TB free,
  and the committed script formats the value correctly.
- The first aggregation run stopped on a missing copy of the discovery ledger on the host;
  the second on the one-cell label conflict above, before the rule in the config was made
  explicit in code; the third on the project's 12 GiB resident-memory guard, which counts
  memory-mapped file pages — the streaming reads were changed to explicit offset reads
  (peak 1.9 GiB). No numbers changed between attempts; only the completed run wrote outputs.
- Clean-up of the 39 GiB object and the 14 GB extract on the secondary volume (plan step
  8) has not been done; both remain, hash-recorded, pending a decision.

## What this does and does not establish

Both arms of the discovery cohort are now in the exact format the v1 model code consumes,
derived from the original public object rather than the authors' compact export, and the
pipeline reproduces the frozen v1 proliferative columns from raw counts. That is an
engineering identity check. It says nothing yet about phase biology or about the positive
control's endpoints: the next protocol step is the discovery `SEC_minus_PRO` fit and the
F1 list freeze, after which — and only after which — the HECA phase model may be fitted.
