# Phase positive-control feasibility — metadata only

Date: 2026-10-26. Protocol step 1 of
[PHASE_POSITIVE_CONTROL_PROTOCOL.md](../docs/PHASE_POSITIVE_CONTROL_PROTOCOL.md), executed
after the protocol freeze (receipt
[`phase_positive_control_freeze.json`](../receipts/phase_positive_control_freeze.json), git
head `b2670a1`). Script: [`scripts/phase_positive_control_feasibility.sh`](../scripts/phase_positive_control_feasibility.sh);
receipt: [`phase_positive_control_feasibility.json`](../receipts/phase_positive_control_feasibility.json).

**Result: both cohorts are estimable under the primary design, and the HECA secretory
arm is ready to run after a 1.38 GB download. The discovery secretory arm has no local
counts and requires the 41.7 GB public GEO object; that download has not started.** No
expression value was read, nothing was downloaded, and no model was fitted in this step.

## What was checked

| Check | Discovery (GSE266265) | HECA (E-MTAB-14039) |
|---|---|---|
| Frozen rule reproduces the frozen candidate IDs | yes, 35 = 22 + 13 | yes, 42 = 22 + 20 |
| Phase groups versus minimum | secretory 13 / proliferative 22, minimum 8: pass | secretory 20 / proliferative 22, minimum 5: pass |
| Disease within the secretory arm | 7 ENDO / 6 non-ENDO; 3 ENDO with adenomyosis | 16 ENDO / 4 non-ENDO; 3 ENDO with adenomyosis |
| Primary design | 13 columns, rank 13, residual df 22: **adopted** | 13 columns, rank 13, residual df 29: **adopted** |
| Fallback design (not used) | 6 columns, rank 6, df 29 | 5 columns, rank 5, df 37 |
| Age | complete, centered on 34.43 years | complete, centered on 32.5 years; H5 and clinical ages agree for all 20 |
| Secretory counts present locally | **0 of 13** (count matrix has the 22 v1 columns only) | **0 of 20** (count provenance lists the 22 v1 donors only) |
| Public source for the missing counts | one object, `GSE266265_EndoAtlas.rds`, HEAD 200, 41,729,517,596 bytes = recorded size | 20 archives, all HEAD 200, 1,378,729,112 bytes in total |

Design files are in the v1 format:
[discovery primary](phase_positive_control_design_discovery_primary.tsv),
[discovery fallback](phase_positive_control_design_discovery_fallback.tsv),
[HECA primary](phase_positive_control_design_heca_primary.tsv),
[HECA fallback](phase_positive_control_design_heca_fallback.tsv). Candidate ledgers:
[discovery](../manifests/phase_positive_control_discovery_candidates.tsv),
[HECA](../manifests/phase_positive_control_heca_candidates.tsv) (the HECA ledger also lists
the three secretory samples not selected and why). Download manifest with per-archive
sizes: [HECA archives](../manifests/phase_positive_control_heca_download_manifest.tsv).

## Decisions fixed by this step

- **Design level: primary in both cohorts.** Rank was computed on the actual numeric
  matrices by Gaussian elimination with partial pivoting (tolerance 1e-7 of the largest
  entry). Adenomyosis varies in both cohorts and stays in the model. Per the protocol,
  this choice may not be revisited after effects.
- **Donor sets.** Discovery: the 13 secretory candidates named in the protocol; sample87
  is excluded by the same minor-clinical-flag rule v1 applied (cervicitis). The 18
  periovulatory donors are not used. HECA: the 20 secretory donors named in the
  protocol; `UA_Endo11460962_FX1170` and `UA_Endo11460963_FX9006` are the smaller repeat
  samples of donors already represented, and `UA_Endo12865966_FX1223` has 49 primary
  nuclei and is below the 50 threshold; none is summed into another sample.
- **Batch structure, for interpretation.** Phase is not confounded with batch in
  either cohort, but it is unevenly distributed. Discovery capture batches 11 and 14
  (6 donors) are proliferative-only. HECA lanes 20 and 43103 (3 donors) are
  proliferative-only and lanes 44595 and 44999 (5 donors) are secretory-only; those 8
  HECA donors contribute no within-lane phase contrast. The phase coefficient is
  identified from the mixed batches and from the batch-adjusted between-batch
  comparison; this is recorded, not corrected.
- **Library quality notes carried, not applied.** As in v1, author library QC strings are
  recorded but are not exclusion criteria; eligibility uses the author per-nucleus
  `discard_lowQC` flag only. Four selected HECA secretory libraries carry a
  "degradation" note and the FX9006 library is annotated as not genotyped by the
  authors; v1 selected libraries with comparable notes.

## What remains before any fit

1. **Discovery secretory counts (protocol step 2).** The OSF compact export covers only
   the author's 23 proliferative donors, so the 13 secretory donors require the full GEO
   object. The server confirms the recorded size. Resource plan: the repository volume
   has about 330 GiB free; the object needs 39 GiB plus extraction space, so a 100 GiB
   reserve is adequate. The archived author stream parser reported a peak of about
   2.2 GiB resident memory on its original host; the local WSL2 guest is capped at 8 GB
   memory and 2 GB swap, so the parser's documented Linux fixes (resident-size units,
   memory ceiling, disk reserve) must be applied and tested on a small serialization
   fixture before the real object is opened. Per-donor count conservation against
   `nCount_RNA`, barcode/metadata order and author-label membership are required checks,
   as for v1. Download time depends on bandwidth and is not estimated here.
2. **HECA secretory counts (protocol step 4).** Twenty public archives, 1.38 GB, through
   the route v1 used; the v1 HECA count-preparation script applies unchanged. Downloading
   them does not compute any effect and is therefore permitted before the F1 list
   freeze; fitting the HECA phase model is not.
3. **Runtime.** The fits use the preserved container runtime described in
   [LOCAL_EXECUTION.md](../docs/LOCAL_EXECUTION.md); the guest was stopped at the time of
   this check and was not started.

If step 2 cannot be completed within the plan above, the protocol's declared reduced
form applies: HECA-only within-cohort detection (P1) and the marker panel (S1), with
the transfer endpoints P2 and S2 reported as not evaluable.

## Interpretation limits

Metadata eligibility is not power. The HECA secretory arm is 80% endometriosis cases,
so disease and phase are partly entangled there even though both are identifiable; the
design adjusts for disease additively and the contrast remains a phase association, not
a causal effect. Composition shifts among stromal states are part of the measured
contrast by design. Nothing in this step bears on the v1 disease result.
