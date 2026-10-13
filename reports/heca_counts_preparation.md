# HECA original-count preparation — public raw route completed

Status: **VERIFIED raw-source intake, fixed-cell correspondence and integer donor aggregation; no disease effects.** This milestone solves the transformed-X input failure using public raw processed files. It does not establish scientific project completion, independent replication, target validity or clinical utility.

## What was completed

Used the 22 exact author-selected samples from manifests/final_heca_clinical_selected.tsv, retaining 15 ENDO and 7 control donors, one chosen library each. Every one of the 52,744 original H5 author-selected primary-stroma nuclei matched a raw source library-plus-barcode. No selected cell was dropped, replaced or reassigned. The smallest selected donor has203 qualifying nuclei.

The previous official H5AD X is standardized log expression and remains ineligible as raw counts; it was not inverted or used for count aggregation. The replacement inputs are the 22 public processed .tar.gz files explicitly listed in cached E-MTAB-14039 BioStudies/SDRF and matching the fixed H5 sample crosswalk. Official source: [E-MTAB-14039](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14039), with count-generation methods in the [original HECA paper](https://doi.org/10.1038/s41588-024-01873-w).

All22 archived matrices have36601original gene features with stable Ensembl IDs, original symbols and Gene Expression modality. Their feature IDs/order are identical. Each matrix has the explicit header “MatrixMarket matrix coordinate integer general”; the complete values, sparse geometry and source identifiers were checked, not just this header.

## Method and source checks

1. Downloaded only the22 exact inventoried files, totaling1,582,044,098bytes; recorded source size, SHA256, HTTP redirect, ETag and Last-Modified. The earlier1MiB route probe agrees with its full archive ETag and prefix hash. No controlled reads, FASTQ, EGA data or complete8.739GB log H5AD was downloaded.
2. Validated complete outer gzip CRC, tar member/sample identity, unique matrix/features/barcode members, inner gzip CRC, all MatrixMarket coordinate bounds, nonnegative integer values, declared nonzero counts and coordinate uniqueness. Read archives without full extraction.
3. Connected raw16-base barcodes through their actual library to the exact author-selected identities. Ambiguous normalization, duplicate or missing selected keys stop processing. The own fixed-cell selection independently matches the clinical index-array hash.
4. Parsed2MiB text blocks and summed counts as exactly represented float64 integers below2^53. Per-gene selected-cell totals equal independent per-cell selected totals; all-source totals equal their per-cell totals. These are count sums, not normalized expression, estimated counts or inverse-rounded log values.
5. Applied the same pinned HGNC2026-10-06 approved-gene policy as discovery, giving priority to original source Ensembl IDs. Retained original symbols and source-to-HGNC disagreement fields. Excluded all unresolved mappings; no alias guess, gene replacement, sum or maximum selection resolves ambiguity.

Author-retained cells are inherited QC. This audit did not rerun Cell Ranger, genotyping, marker-based annotation or disease tests. Raw-source correspondence and checksum checks are engineering/provenance evidence; they do not independently certify biological donor identity.

## Actual results and denominators

| Item | Observed result | Meaning |
|---|---:|---|
| Final donors |22:15ENDO +7controls | Exact clinical-selected biological units; extra libraries do not add donors. |
| Qualifying nuclei |52,744/52,744 matched | Complete fixed author-cell coverage, with no silent loss. |
| Original source genes |36,601 | Original raw archive features, preserved before mapping. |
| HGNC retained stable genes |26,125 | Unique official stable targets in the named-gene annotation universe. |
| Mapping exclusions |10,474unmapped +2ambiguous | Missing official mapping is unavailable evidence, not proof of no biology. |
| Raw selected integer UMIs |360,765,431 | Original counted molecules summed over the fixed selected nuclei. |
| Retained stable-gene UMIs |352,307,494 |8,457,937UMIs excluded by the declared mapping rule. |
| Mapping loss |2.34444% overall;1.76777–3.43949% per donor | Small total mass loss does not establish unimportance of omitted genes. |
| Stable-ID overlap with discovery |21,585 | Label-blind identifier overlap only; this is not the frozen expression/graph comparison universe. |

The annotation covers approved named genes, not the complete human genome feature catalog. Counts are preserved in source_counts.tsv.gz even for excluded features; the mapping decision ledger is manifests/heca_features_mapping.tsv. No donor/gene selection used disease expression effects.

## Raw-source versus H5 original total differences

The raw archives preserve36601features, whereas H5 var has33578after author filtering. Original H5 obs/n_counts was retrieved under the already approved small probe budget and preserved unchanged.

Among52,744matched nuclei,44,434have exactly equal raw/H5-original totals;8,310have raw totals higher by1–8UMIs; none has a lower raw total. Histogram:

| Raw minus H5-original UMIs per nucleus | Nuclei |
|---:|---:|
|0|44,434|
|1|7,230|
|2|916|
|3|134|
|4|19|
|5|5|
|6|2|
|7|2|
|8|2|

Raw total is360,765,431versus H5-original360,755,824: difference9,607UMIs,0.00266295% of raw total. Source author loading code applies filter_genes(min_cells=3), which provides a plausible explanation for library-specific differences; that attribution is **INFERRED, not completely proved**.

An optional zero-network, zero-MTX-reparse check used only original source symbols unique within the raw feature table and exactly absent from H5 var. Those3022features explain1359UMIs, not all9607. Twenty raw feature rows have nonunique symbols and were excluded from that optional correspondence check rather than guessed. Global var absence alone is therefore insufficient; precise library-specific filtering attribution remains unresolved.

No equality was forced, cells were not deleted to make totals agree, counts were not changed, and no inverse reconstruction was performed. All exact differences remain in raw_vs_author_H5_cell_totals.tsv.gz. A small discrepancy is not an independent identity proof or an assurance about individual biological genes.

## Verification and runtime

Final hardened run on SRV01 completed in332.213seconds with peakRSS541576KiB (approximately529MiB). Free disk afterward was67,666,747,392bytes (approximately63.02GiB), above the20GiB reserve and far below the12GiB memory ceiling. The two-stream initial download completed in1238.528seconds, with exact1.582GB new source bytes; original source files remain in ignored project data directories.

13meaningful boundary tests passed. They cover missing/duplicate/ambiguous barcodes, negative or noninteger source values, duplicate coordinates, wrong source feature identities, changed selected-barcode pins, mutated NPZ contents, fractional cached arrays and deterministic gzip bytes. The original eight-test baseline failed before implementation and is retained in logs/heca_counts_tests_before.log. Independent final validation read all emitted source/stable counts and confirmed hashes, integer/nonnegative values, mapping correspondence, complete cell counts and per-donor/raw/retained/excluded mass conservation.

Gzip outputs use fixedmtime=0 and emptyfilename. Checkpoints bind archive hash, clinical hash, selected-barcode list, feature ordering, count-parser implementation, NPZ hash and schema version; reused arrays must pass geometry, finite/nonnegative/integer and count-mass checks. Old unpinned checkpoints received a one-time local numerical reparse; no source download was repeated. The execution script hash was captured before the final run and checked against disk before publishing its receipt. Final outputs are quiescent.

Model request: gpt-6.1-sol, high. The parent model_dispatches.json records explicit delegation configuration; this is not independent backend model attestation. Actual execution host, Python3.10.16, numpy2.2.6, requests2.32.5, wall time, memory and file hashes are in the runtime/source receipts. No fallback model or external model command was used.

## Deliverables and remaining limits

Primary interface (firstcolumn gene_id, then22exact original clinical donor IDs in fixed clinical order):
- data/processed/heca/source_counts.tsv.gz —36,601original features.
- data/processed/heca/counts_stable.tsv.gz —26,125mapped stable IDs; SHA25677121542de07e346ea2112f80bfee6248993c42f90c393ec3154b0a9de529af0.
- data/processed/heca/selected_donor_count_ledger.tsv —cells, source UMIs, mapping loss and original clinical fields.
- data/processed/heca/provenance.json; receipts/heca_counts_preparation.json; receipts/heca_counts_validation.json.
- manifests/heca_features_mapping.tsv; data/processed/heca/source_features.tsv; shared_stable_features.tsv; raw_vs_author_H5_cell_totals.tsv.gz.
- reports/heca_counts_raw_inventory.json; receipts/heca_counts_raw_downloads.json; receipts/heca_counts_raw_route_probe.json; receipts/heca_counts_original_totals.json; receipts/heca_counts_filter_difference_probe.json.

Executed:
```
envs/m0-h5/bin/python scripts/prepare_heca_counts.py --download-only
envs/m0-h5/bin/python scripts/prepare_heca_counts.py --wait-for-downloads
envs/m0-h5/bin/python scripts/prepare_heca_counts.py
envs/m0-h5/bin/python -m unittest discover -s tests -p test_prepare_heca_counts.py
envs/m0-h5/bin/python scripts/validate_heca_counts.py
```

The earlier transformed-X failure and initial metadata-budget accounting error remain documented in reports/heca_counts_layer_audit.md; they were not turned into passes. The public raw route now supplies a genuine original-count input, with exact correspondence to author-selected cells. Source filtering attribution is still incomplete; fixed-marker plausibility, complete M0/design/independence adjudication, gene-expression/network freezes and scientific inference remain parent-task gates. No disease ranking, target validation or clinical result exists from this round.

No new user decision is required for this already-authorized preparation. YELLOW — understanding not yet assessed: distinguish original counts from log expression, source-gene filtering from cell filtering, exact integer summation from scientific validation, and stable-ID overlap from a final comparison universe.

