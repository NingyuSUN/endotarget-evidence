# Discovery secretory counts — acquisition plan (protocol step 2)

Date: 2026-10-26. Status: **EXECUTED 2026-10-26** — see the
[discovery counts report](../reports/phase_positive_control_discovery_counts.md) and
[receipt](../receipts/phase_positive_control_discovery_counts.json). No model was fitted.
The plan below is kept as written before execution; it implements section 10 of the
[phase positive-control protocol](PHASE_POSITIVE_CONTROL_PROTOCOL.md) for the discovery
cohort and records the measured resource constraints before the large transfer, as the
repository rules require. Deviations from the plan during execution are listed in the
report.

## Why a new download is unavoidable

The v1 discovery matrix was built from the public OSF compact export, which aggregates
only the 23 author-selected proliferative donors. GEO lists no per-sample count files for
GSE266265; the only public object containing the 13 secretory candidate donors is the
full atlas `GSE266265_EndoAtlas.rds`. The
[feasibility step](../reports/phase_positive_control_feasibility.md) confirmed by a HEAD
request that the object is served at its recorded size.

## Source object, as recorded by the author route

| Item | Value | Provenance |
|---|---|---|
| URL | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266265/suppl/GSE266265_EndoAtlas.rds` | GEO series supplementary file |
| Size | 41,729,517,596 bytes | author route record; HEAD 2026-10-26 agrees |
| SHA-256 | `44a832cdbb7b43ebd3b4a2b6e375c52b08a8059229cc60228ab20a29249cf94f` | author route record; to be verified after download, never assumed |
| RNA counts | 34,727 source features × 466,371 cells, 60 samples; 1,188,298,946 nonzero estimates, of which 104,078,188 fractional | author route record ([route report](../reports/osf_route_counts_and_stream_reader.md)) |
| Metadata used | `sample`, `AnnotationMain`, `AnnotationRefined`, `X10x.capturing.batch`, `nCount_RNA` | author aggregation code and crosswalk |

The author stream parser (`rds_stream.py`, `aggregate_author_counts.py`) is archived under
the OSF route inputs. It reads the serialization in 8 MiB chunks, keeps the first assay's
counts and `meta.data`, and does not execute R code. Its recorded peak on the authors'
host was about 2.2 GiB resident memory and 155 s. Three documented porting problems must
be fixed in a project copy before use: Linux `ru_maxrss` is reported in KiB, not bytes; the
built-in 35 GiB memory ceiling exceeds the analysis host's physical memory; and the
75 GiB disk-reserve rule needs a volume that can honor it.

## Measured resources on the native analysis host (SRV01)

| Resource | Measured 2026-10-26 | Consequence |
|---|---|---|
| Memory | 31 GiB physical, about 28 GiB available; one unrelated job running | Set the parser ceiling to 24 GiB; keep the project 12 GiB RSS guard for aggregation |
| Cores | 36 | Not limiting; the parser is single-threaded |
| Root volume | 42 GB free (91% used) | **Not usable**: a 39 GiB object would breach the project's 20 GiB reserve rule |
| Secondary volume | about 1.3 TB free | Use it for the object, the streamed payload and scratch; budget 100 GiB |
| Transfer rate | GEO, 2026-10-26: a 200 MiB range probe ran at 7.3 MB/s single-stream; the full object then transferred at 17.9 MB/s aggregate with four ranged streams (2,329 s), plus 662 s to concatenate the parts and 604 s to hash the result on the secondary volume | Four ranged, resumable streams are adequate; no multi-day budget is needed |

Disk budget: 39 GiB object + about 15 GiB streamed sparse payload (1.19 × 10⁹ nonzeros
as 32-bit indices plus 64-bit values, with pointers) + extracted metadata + reserve.
100 GiB is comfortable on the secondary volume and impossible on the root volume.

## Procedure, in order

1. **Pre-flight.** Record free space on the target volume (≥100 GiB required) and memory.
   Verify the archived parser scripts against the OSF member hashes. Apply the three
   fixes in a project copy under `scripts/` (archived originals untouched). Test the copy
   on a small R serialization fixture written by the project R runtime: a two-assay
   list with a sparse integer/fractional matrix and a data frame, checked for exact
   round-trip of values, barcode order and column names. No real object is opened
   before this passes.
2. **Download** to the secondary volume with a resumable ranged transfer, `Content-Length`
   checked against the recorded size before and after, then the full SHA-256 compared to
   the recorded value. Mismatch stops the step; the file is not "repaired". Record ETag,
   Last-Modified and wall time in a receipt.
3. **Stream-extract** the first assay's counts and `meta.data`. Check: 466,371 barcodes in
   metadata order, 34,727 features, nonzero total and fractional count equal to the
   author record, no negative or non-finite value, `orig.ident`/`sample` consistent with
   barcode prefixes.
4. **Aggregate** per donor for the 35 analyzed donors only (22 v1 proliferative + 13
   secretory) with the frozen lineage rule in
   [`config/discovery_lineage.json`](../config/discovery_lineage.json): `AnnotationMain`
   stromal, refined label in the 12-label list, mural and missing labels excluded. For
   the 22 v1 donors, per-donor cell membership must equal the author
   `DONOR_CELLTYPE_MANIFEST` counts exactly; for the 13 secretory donors, record cell
   counts and apply the ≥50 rule. Per-donor sums must equal the metadata `nCount_RNA`
   totals of the selected cells (relative tolerance 1e-12, absolute 1e-6, as in v1).
5. **Feature handling identical to v1.** Drop the 177 source features that form the 50
   irreversible multi-source canonical groups, so the singleton universe (34,550) and the
   pinned HGNC mapping (26,109 stable IDs) are the same as v1. Fractional estimates are
   kept unrounded.
6. **Identity check.** The 22 v1 donor columns must agree with the frozen discovery
   matrix within relative tolerance 1e-6 (the OSF export printed float64 sums; re-summing
   in float64 may differ in the last digits, so string identity is not expected). A larger
   difference stops the step and is reported.
7. **Outputs** go to a new tree, `data/processed/phase_positive_control/discovery/`, with
   a receipt at `receipts/phase_positive_control_discovery_counts.json` recording input
   and output hashes, the v1 matrix hash before and after (unchanged), peak memory and
   wall time. The v1 tree is never written.
8. **Clean-up.** Delete the 39 GiB object and the streamed payload after the receipt is
   final, recording their hashes and the deletion in the receipt. Per-donor sums are kept.

## Stop conditions and the declared reduced form

Hash mismatch; dimensions or totals that differ from the author record; any v1 column
beyond tolerance; parser memory above the ceiling; or fewer than 8 secretory donors with
≥50 primary-lineage cells after extraction. If the step cannot be completed within this
plan, the protocol's reduced form applies: HECA-only P1 and S1, with P2 and S2 reported
as not evaluable. The cohort roles are not swapped and no threshold is relaxed.
