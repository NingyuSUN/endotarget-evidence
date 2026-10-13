# Task 1 bounded metadata retrieval and donor audit
Status: VERIFIED engineering metadata milestone; Stage 0 OPEN; scientific eligibility NOT VERIFIED.

## Scope and method
All downloads, code, tests and audit execution occurred on SRV01 in /home/analyst/te, Python 3.10.12 stdlib. SSH alone was launched from the controlling workspace. No expression matrices, supplements, DE, cloud/container run, sudo or external model CLI were used.

Sources: SDRF from source_inventory.csv, BioStudies JSON from docs/DATA_SOURCES.md; official GEO canonical record initially rejected as a CAPTCHA response. The parent authorized the standard accession-derived official HTTPS FTP family SOFT alternate, which was then actually fetched and validated. Exact original/resolved URLs, UTC, bytes and SHA256 are in receipts and manifests/inputs.lock.tsv.

Downloads use HTTPS certificate verification, HTTPS-only redirects, a 30-second socket timeout, 120-second streaming deadline, <=10 MiB compressed transfers and <=10 MiB decompressed SOFT. Exclusive .part files are validated then atomically renamed. Existing metadata requires matching hash, bytes, type and URL to be reused. --refresh preserves original pins and saves a changed payload separately, returning failure. This source-change behavior was unit-tested; no actual upstream change is claimed. Rejected diagnostic payloads are preserved remotely with receipts.

## VERIFIED actual results

SDRF: 103 library rows ({'scRNAseq': 35, 'snRNAseq': 68}). Exact deduplication key: Characteristics[10x technology] + Characteristics[individual]. Nuclei assay: 63 donor units; cell assay: 16 donor units. The sum 79 is NOT 79 independent people; shared IDs across assays are deliberately not merged.

| Assay | Exact source disease label | Proliferative candidates |
|---|---|---:|
| snRNAseq | endometriosis | 16 |
| snRNAseq | normal | 6 |
| scRNAseq | endometriosis | 5 |
| scRNAseq | normal | 5 |

Replication uses nuclei: 22/63 donor units match the exact disease and phase labels, before final eligibility exclusions. This reproduces the historical SDRF hash/count evidence, not clinical or biological validity. The normal label does not certify health.

Within-donor clinical conflicts: 0/79 assay-specific donor units. Technical library_id values are preserved separately and not treated as clinical contradictions. One nuclei donor has Unknown coarse/fine phase. This does not establish completeness of other clinical fields: hormone history and other protocol exclusions are not independently established by this table.

Filename screening flagged 34/103 library rows ({'scRNAseq': 6, 'snRNAseq': 28}): the declared donor ID is absent as a literal filename token. Original labels, filename, field and line number are retained as UNRESOLVED_NAMING_EVIDENCE. No donor was reassigned. This broad screen is not a diagnosis of sample swaps.

GEO: 60 GSM records, 8,459 compressed bytes, 362,489 decompressed bytes. Sample records are not verified eligible donors.

| Source | Bytes | SHA256 |
|---|---:|---|
| sdrf | 61604 | cdcd25f5519b3ddac7b68a41140ae577b929a856c35732e598fc95fcb0bbf36b |
| biostudies | 174654 | 0d91686fe91ff2223c85821e5ac511fd4d254762c953ed40e0dd594d0f57abdb |
| geo_soft | 8459 | 2b73147d402f758b00fe74c578600ed6a678448e546b3f65e178ca36852c41b7 |

## Validation and problems
VERIFIED: python3 -m unittest discover -s tests -v ran 22 tests on SRV01, exit 0. Coverage includes technical-library deduplication, assay-specific identity, clinical conflicts/missingness, required columns, ragged/non-table payloads, HTTPS policy, bounded streaming, gzip expansion limits, accession validation, cache tampering/reuse and changed-source preservation. Final real replay returned reused_validated for all three sources.

Pre-implementation collection initially failed because the module did not exist; this is not 17 separately executed assertion failures. Two subsequent source regression tests actually failed before fixes: exact cycle header mismatch and treating library_id differences as clinical conflicts. SOFT positive-format support also had a failing test before implementation. All 22 pass now. Raw outputs, hashes and failure history are in receipts/m0_metadata_test_evidence.json and run receipts.

Problems and resolution:
1. Assumed phase header did not match the real Characteristics[menstrual_cycle_stage]. Strict validation failed; original header was inspected and corrected.
2. Technical library_id differences were initially treated as clinical conflicts. A failing regression demonstrated the defect; library IDs now remain separate.
3. GEO canonical HTML returned reCAPTCHA. It was rejected, not solved or accepted as metadata; official family SOFT succeeded.
4. One SDRF request closed without response; one bounded retry succeeded. Both runs are preserved.
5. Server rg is unavailable; read-only grep was used instead.

## Exact execution and outputs
All commands ran at /home/analyst/te:
    python3 -m unittest discover -s tests -v
    python3 scripts/m0_metadata.py

Test stdout/stderr was redirected to logs/m0_metadata_tests_{red,regression_red,soft_red,final}.log; final script replay to logs/m0_metadata_cache_replay.log. The test-evidence receipt embeds these logs. --refresh is documented optional usage and unit-tested, but was not run as a real remote refresh this round.

Main outputs:
- scripts/m0_metadata.py and tests/test_m0_metadata.py
- manifests/inputs.lock.tsv
- reports/m0_metadata_audit.json: exact fields, labels, strata, conflicts and naming evidence
- reports/m0_metadata_readiness.json: stage0_passed=false and final_eligibility=NOT_VERIFIED
- reports/m0_metadata_summary.md and remote ignored work/task-1-report.md
- receipts/m0_metadata_{sdrf,biostudies,geo_soft}.json, run/rejection receipts, test evidence
- ignored raw metadata under data/metadata and logs under logs/

## Open gates and user understanding
NOT VERIFIED: donor-cell joins; >=50 eligible stromal cells/nuclei per donor; hormone/clinical exclusion; raw counts; age/batch design; independent donor provenance; resource feasibility; source terms. Discovery >=8/group and replication >=5/group apply after all exclusions and final design checks. No threshold changed. Stage 0 cannot pass from this audit alone.

Task 2 must resolve the clinical and cell mapping and safe import path before scientific inference. Engineering test success is not scientific/clinical/target validation. No new user decision is needed for this bounded milestone. YELLOW: understanding not yet assessed; explain why 103 libraries, 63 nuclei donor units and 22 phase/disease candidates answer different questions.
