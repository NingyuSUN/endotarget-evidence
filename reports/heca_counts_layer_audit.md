# HECA counts-layer input audit — .X fails raw-count gate

Status: **VERIFIED transformed input; raw-count readiness FAILED for this public H5AD .X.** Pseudobulk preparation has not run. No disease coefficient, gene ranking, expression target list or DE table was computed or read.

## Actual input and layer evidence

Source is the same author H5AD already audited for cell identity: https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_nuclei.h5ad ; size 8,739,186,814 bytes; ETag `"3b257ce77e9a967b574b55810dfa189d-556"`. Count-probe range responses retain exact hashes and Last-Modified in receipts/heca_counts_schema.json and receipts/heca_counts_numeric_probe.json.

Root keys are X, obs, obsm, obsp, uns and var. There is **no raw and no layers**. X is sparse CSR, shape 312,246 nuclei × 33,578 features, with 1,051,326,046 nonzeros. Data are float32 and column indices int32. Both arrays have chunks of 64,168 values and no compression.

| Test | Actual observed result | Interpretation |
|---|---|---|
| Three samples of 1,024 nonzero values, beginning/middle/end | All finite and nonnegative; all fail integer-only check | Finite/nonnegative does not establish count provenance. |
| Five selected primary-stroma rows | All sparse indices in bounds and unique | Engineering consistency only. |
| Five row totals versus original author obs/n_counts | All differ; row 1 has X sum 3,244.9323 versus original 4,967 | X cannot be summed as unnormalized UMIs. |
| Author pipeline code | Clustering notebook executes normalize_per_cell(counts_per_cell_after=1e4) and log1p | A concrete provenance explanation for the transformed values. |

The primary paper states that Cell Ranger 6.0.2 filtered count matrices were used, with GRCh38-2020-A. The pinned author code loads these filtered 10x matrices, writes an internal nuclei_rawcounts.h5ad, reads that file for integration, and later normalizes/log-transforms during clustering. The public V2 export-writing step was not observed. Thus the Cell Ranger statement supports upstream raw counts, not a claim that the downloadable X is still raw.

Author source: [HECA paper](https://doi.org/10.1038/s41588-024-01873-w), local data/reference/heca_article.xml; fixed public GitHub tree 244c282c20c3efd53707e3eb9f58bbd1cbd56ad7. Original notebooks and hashes are retained in ignored data/metadata/heca_h5/heca_count_authorcode/ and receipts/heca_counts_authorcode.json. Only code cells concerning matrix input/transformation were inspected; notebook outputs and disease-DE tables were not inspected.

## Fixed selected-cell input and measured resource plan

Clinical branch supplied manifests/final_heca_clinical_selected.tsv and work/final_clinical/heca_primary_h5_obs_indices.npy: 22 chosen samples, 15 ENDO and 7 controls, one library per donor, 52,744 selected original-row indices. Manifest SHA256 1a0ba5813046f9cb6e18862e17ca7fe64c31a39605d048309f1e3522447046cc; index NPY SHA256 d02ed9b8cee2c789c4f7aec7a3792a8f9444307163531933b632845cecf864ef. The clinical receipt owns selection rationale; this audit does not change it.

The selected rows contain 171,669,473 stored nonzeros and touch 15,955 of 16,384 chunks in each of data and indices. Exact selected-chunk storage is 8,190,403,520 bytes across both arrays, **97.39%** of their complete storage. Sparse row selection would save little network traffic for this file. At the measured host snapshot, approximately 65 GB disk was free, and 29 GiB RAM available; retaining a full 8.739 GB H5AD would leave more than the 20 GiB disk reserve, and chunked local aggregation could remain below the 12 GiB memory guard. This is a resource estimate, **not authorization to download transformed X**. No full object was downloaded.

## Feature identity recovered without additional requests

Var index contains 33,578 unique source symbols, not Ensembl IDs. There are 51 per-input categorical columns gene_ids-0 through gene_ids-50, introduced by author concatenation. Decoding every duplicated dictionary inflated metadata retrieval; that broad scan stopped at its source-byte ceiling.

Using only cached bytes, 22 of these 51 gene-ID columns were decoded fully. Of 33,578 rows, 33,450 have one observed stable Ensembl ID and 128 have no observed ID; no conflicts occurred among the 22 read columns. The remaining 29 dictionaries were not asserted concordant. Source symbols and observed IDs are in data/metadata/heca_h5/heca_count_probe/features_cached.tsv, SHA256 8437837dd60a0f9a273ee80ef869f2aaeee51913b7aa07896261355c6d959056. Its explicit column-coverage fields prevent a partial correspondence from appearing complete.

Do not select an arbitrary gene_ids column or strip numeric symbol suffixes by guess: suffixes may originate from unique-name handling or be part of an actual symbol. Preserve original index, all observed original IDs, and mapping coverage. The shared pinned HGNC rules must retain source-ID priority and expose source-ID versus approved-symbol disagreement. Gene mapping is a separate branch; this report does not establish a final common gene universe.

## Inverse-transform feasibility is not approved counts

For five already-cached rows only, a PROPOSED inverse expm1(X) × obs/n_counts / 10,000 produced estimates close to integers: maximum nearest-integer deviation across the five rows was 0.00010448. Nearest-integer row sums matched all five original totals. This is limited mathematical feasibility, conditional on the exact normalization scale, unmodified total-count metadata, complete feature retention and floating-point behavior.

It does not establish original raw-count recovery across the dataset. No counts have been accepted, rounded for analysis, or aggregated. The protocol prohibits silently interpreting transformed values as counts. The parent task explicitly directs searching public unnormalized processed files before considering a revised reconstruction method.

## Problems, budget accounting and verification

1. Broad var-category decoding reached exactly 200 MiB of H5 range bytes and stopped after 218.125 seconds. The error text inherited from the prior probe says “150 MiB”; its actual configured ceiling and observed stop were 200 MiB. No byte beyond that H5 ceiling was fetched.
2. Three author notebooks added 17,575,153 bytes (16.76 MiB), which were not included in that shared initial budget. The initial round therefore exceeded its intended total network allowance by that amount. This is an accounting error, explicitly retained; subsequent source requests stopped until a measured adjustment was approved.
3. The parent authorized a separate <=10 MiB numerical/row-pointer probe. Its exact new transfer and per-range hashes are in receipts/heca_counts_numeric_probe.json; that sub-budget passed. Cached feature recovery and inverse feasibility used zero new network bytes.
4. The sparse indptr is nondecreasing, starts at zero and ends at the declared nonzero count. Negative test verifies normalized-log values are rejected by the raw-integer gate; a synthetic inverse fixture is labeled feasibility only. Syntax and executed self-test passed on SRV01. These checks do not override the real source failure.

Commands:
```
envs/m0-h5/bin/python scripts/heca_count_layer_audit.py --mode self-test
envs/m0-h5/bin/python scripts/heca_count_layer_audit.py --mode cached-features
envs/m0-h5/bin/python scripts/heca_count_layer_audit.py --mode numeric-probe
envs/m0-h5/bin/python scripts/heca_count_layer_audit.py --mode inverse-feasibility
```
Default mode is offline self-test. Numeric probe has an explicit additional 10 MiB source ceiling; it is not a general large-transfer command.

## Current conclusion and next authorized action

M0 raw-count input remains OPEN. The clinical sample selection may pass its own checks while this matrix input fails; these are different gates. Investigate exact public unnormalized processed files already named in cached BioStudies/SDRF and author source inventories. A genuine library matrix can be restricted by the existing library-plus-barcode and selected-row crosswalk without retrieving controlled FASTQ/EGA inputs. No inferential analysis may use this X directly.

This branch has no new user decision. The parent task owns the public-raw-route investigation and any explicit future method revision. YELLOW — understanding not yet assessed: explain why logarithmic standardized values cannot be added as original molecule counts, why five successful inverse rows do not validate a complete reconstruction, and why selected sparse rows can still require nearly the whole source file because chunks mix selected and unselected cells.

