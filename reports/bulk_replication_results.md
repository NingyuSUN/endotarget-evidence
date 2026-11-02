# Program-level replication in public bulk eutopic endometrium — results

**Labels under the frozen rules: E1 `no_detectable_difference` in both analyzed series and
`not_replicated` across series; E2 `does_not_transfer` for the frozen B0 top-20 in every
contrast; E3 descriptive.** In 144 natural-cycle whole-tissue samples of GSE51981 the secretory
program scored from the frozen single-cell phase signatures rises by about half the stromal
single-cell shift between proliferative and secretory tissue (gate AUC 0.96), and endometriosis
does not detectably change that rise: θ −0.09 program units [−0.21, +0.03], p 0.16, which
excludes an attenuation larger than 40% of the tissue-level shift and an enhancement larger
than 7%. In the 37 samples of GSE6364 the same estimate is −0.15 [−0.34, +0.04], p 0.13. The
frozen single-cell top-20 list's directions do not hold in whole tissue (12 / 20, 7 / 20 and
4 / 20 same direction; matched percentiles 81, 21 and 0.3). The per-gene proliferative contrast
of GSE51981 finds 9,593 genes at BH q < 0.05 against all non-endometriosis samples and none
against the other-pathology controls alone, reproducing the pattern reported by the series'
authors and leaving the v1 within-stroma estimand untouched.

Run on 2 November 2026 under the [protocol](../docs/BULK_REPLICATION_PROTOCOL.md) frozen the
same day before any expression value was read (freeze
[`bulk_replication_freeze.json`](../receipts/bulk_replication_freeze.json), git head 6ec376d).
Three runs are published: the first run and two corrections of defects found afterwards
(section [Runs and corrections](#runs-and-corrections)). The numbers in this report are those of
the second correction, [`bulk_replication_correction2/`](bulk_replication_correction2/), unless a
run is named. All runs executed once each on the analysis host in the v1 statistical runtime
(Python 3.10.22, NumPy 2.2.6, SciPy 1.15.2); the expression matrices were downloaded there after
the freeze with a [data receipt](../receipts/bulk_replication_data.json).

## Cohorts and samples

| Series | Platform | Authors' processing | Samples used | Proliferative ENDO / non-ENDO | Secretory ENDO / non-ENDO | Role |
|---|---|---|---:|---:|---:|---|
| GSE51981 (Tamaresis et al. 2014) | GPL570 | GCRMA log2 series matrix, range 2.2–16.5 | 144 (4 of indeterminate phase excluded) | 29 / 35 (15 other uterine or pelvic pathology, 20 no pathology) | 46 / 34 (early 18 / 12, mid 28 / 22) | primary |
| GSE6364 (Burney et al. 2007) | GPL570 | authors' series matrix, range 0.1–16.7 | 37 | 6 / 5 | 15 / 11 (early 6 / 3, mid 9 / 8) | secondary |
| GSE141549 (Gabriel et al. 2020) | GPL10558 + GPL13376 | batch-corrected workbook, downloaded, not parsed | — | — | — | `not_analyzed`, as frozen |

19,823 Ensembl genes per series after probe mapping through the GPL570 annotation (Entrez ID,
then symbol) and the pinned HGNC set, one probe per gene by highest mean expression; each
series is centered on its proliferative non-endometriosis samples. Sample table:
[`samples_used.tsv`](bulk_replication_correction2/samples_used.tsv).

## E1 — secretory program disease response at tissue level

Program score per sample = β-weighted projection on the frozen single-cell phase program
(HECA-derived primary, 3,986 genes, 3,574 in the bulk family; discovery-derived secondary,
7,367 / 6,537), in units where 1 = the signature cohort's full stromal secretory shift. Gate:
secretory-versus-proliferative coefficient positive with AUC ≥ 0.9. Primary θ = secretory ×
endometriosis interaction over all natural-cycle samples (OLS); label `attenuation_detected`
or `enhancement_detected` required two-sided p < 0.025.

| Series, program | Gate (shift, AUC) | Score means: P non-ENDO / P ENDO / S non-ENDO / S ENDO | θ [95% CI], p | MDE (80% power) | θ as share of the tissue shift (CI) | δ within secretory [95% CI], p | Label |
|---|---|---|---|---:|---|---|---|
| GSE51981, HECA program | 0.52, 0.96 | 0 / 0.10 / 0.57 / 0.58 | **−0.087** [−0.209, +0.035], 0.16 | 0.174 | −40% to +7% | +0.015 [−0.082, +0.113], 0.76 | `no_detectable_difference` |
| GSE51981, discovery program | 0.40, 0.96 | 0 / 0.05 / 0.42 / 0.43 | −0.045 [−0.139, +0.048], 0.34 | 0.133 | −35% to +12% | +0.004 [−0.073, +0.080], 0.92 | `no_detectable_difference` |
| GSE6364, HECA program | 0.49, 1.00 | 0 / 0.03 / 0.57 / 0.45 | **−0.148** [−0.341, +0.045], 0.13 | 0.273 | −70% to +9% | −0.115 [−0.239, +0.009], 0.068 | `no_detectable_difference` |
| GSE6364, discovery program | 0.37, 0.99 | 0 / 0.03 / 0.44 / 0.34 | −0.141 [−0.300, +0.019], 0.082 | 0.226 | −82% to +5% | −0.110 [−0.212, −0.007], 0.037 | `no_detectable_difference` |

Across series: `not_replicated` (no series detects an effect). The sub-phase-adjusted δ
(early versus mid secretory as a covariate; a prespecified sensitivity that never changes a
label) is +0.03 [−0.05, +0.10] in GSE51981 and −0.08 [−0.14, −0.02], p 0.008, in GSE6364 for
the HECA program; the GSE6364 secretory arm is small (15 / 11) and its early-secretory
subgroup (6 / 3) carries most of the negative estimate (see the first run below).
Per-sample scores: [`program_scores_GSE51981_heca_program.tsv`](bulk_replication_correction2/program_scores_GSE51981_heca_program.tsv)
and the three companion files.

The whole-tissue shift along the stromal program is about half the single-cell stromal shift
(0.52 and 0.49 for the HECA program), as expected when epithelial, immune and vascular
compartments dilute a stromal signature; the AUCs of 0.96–1.00 show that the frozen programs
separate phases in whole tissue as well as they did across single-cell cohorts.

## E2 — direction transfer of the frozen single-cell top-20

Per-gene disease coefficient in proliferative samples (moderated OLS, section E3), signs
compared with the frozen B0 discovery signs; `transfers` required ≥ 16 / 20 same direction and
a percentile ≥ 95 against 10,000 expression-quintile-matched random 20-gene sets (seed
20261102).

| Series, contrast | Same direction / evaluable | Matched reference mean | Matched percentile | Bulk q < 0.05 among the 20 (in the discovery direction) | Label |
|---|---:|---:|---:|---:|---|
| GSE51981, ENDO vs all non-ENDO | 12 / 20 | 9.6 | 80.8 | 13 (7) | `does_not_transfer` |
| GSE51981, ENDO vs other-pathology controls | 7 / 20 | 8.0 | 20.8 | 0 | `does_not_transfer` |
| GSE6364, ENDO vs all non-ENDO | 4 / 20 | 9.6 | 0.3 | 0 (smallest q 0.65) | `does_not_transfer` |

All 20 genes were measured in both series (missing-outcome bounds equal the point values).
In GSE51981 against all non-endometriosis samples, 13 of the 20 frozen genes reach q < 0.05 in
whole tissue, but only 7 of them in the single-cell direction: the list is sensitive to the
tissue-level contrast without carrying its directions, consistent with the global shift
described under E3 rather than with a replicated stromal effect. Gene-level values are in
`summary.json` (`E2.*.genes`).

## E3 — tissue-level per-gene proliferative contrast (descriptive)

Moderated t per gene (limma-style empirical Bayes on the gene-level log2 matrix; prior df
estimated by the method of moments with the corrected trigamma inverse).

| Series, contrast | n ENDO / non-ENDO | Residual df | Prior df | Genes at BH q < 0.05 | Median MDE (log2, fold) | Share of 95% CIs within ±1 log2 | Spearman of bulk t with pooled single-cell t (14,202 shared genes) | Sign agreement |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE51981, ENDO vs all non-ENDO | 29 / 35 | 62 | 0.18 | **9,593** (6,812 with \|β\| ≥ 0.585; 2,285 with \|β\| ≥ 1; 61.5% negative) | 0.51, 1.42 | 60.5% | 0.08 | 55.2% |
| GSE51981, ENDO vs other-pathology controls | 29 / 15 | 42 | 0.17 | **0** (smallest q 0.145; 226 genes with \|β\| ≥ 1) | 0.62, 1.54 | 78.7% | −0.22 | 43.2% |
| GSE6364, ENDO vs all non-ENDO | 6 / 5 | 9 | 2.3 | 1 | 0.79, 1.73 | 62.5% | −0.07 | 47.5% |

Per-gene tables: [`proliferative_de_GSE51981_all_non_ENDO.tsv.gz`](bulk_replication_correction2/proliferative_de_GSE51981_all_non_ENDO.tsv.gz),
[`proliferative_de_GSE51981_other_pathology_controls_only.tsv.gz`](bulk_replication_correction2/proliferative_de_GSE51981_other_pathology_controls_only.tsv.gz),
[`proliferative_de_GSE6364_all_non_ENDO.tsv.gz`](bulk_replication_correction2/proliferative_de_GSE6364_all_non_ENDO.tsv.gz).

Reading. The 9,593 genes against all non-endometriosis samples are driven by the 20
no-pathology controls: the series' authors reported that endometrium from women without
endometriosis and without other pelvic pathology differs broadly from both disease groups,
and a descriptive check on the first run found a mean log2 shift of −0.11 across all genes in
that contrast. Against the 15 other-pathology controls (the contrast closest to the v1 design,
whose non-ENDO donors had other gynecological indications) no gene reaches q < 0.05 and 79%
of the gene intervals lie within a twofold change. Bulk and pooled single-cell statistics are
essentially uncorrelated (Spearman 0.08 and −0.22; sign agreement at chance). The prior df
near 0.2 in GSE51981 reflects about 1,600–1,800 genes floored at the GCRMA minimum in nearly
every sample (residual variance near zero); moderation is therefore negligible, and only 8 of
the 9,593 significant genes are such floor genes.

## Runs and corrections

The protocol fixes the order: rules frozen → data downloaded with a receipt → one run. Defects
found after expression values were read do not alter a published run; corrections are
published separately with their own freeze and results receipts, and the earlier outputs stand
as reported.

| Run | Freeze (git head) | Outputs | Change | Effect on results |
|---|---|---|---|---|
| First run | [`bulk_replication_freeze.json`](../receipts/bulk_replication_freeze.json) (6ec376d); results [`bulk_replication_results.json`](../receipts/bulk_replication_results.json) | [`bulk_replication/`](bulk_replication/) | — | E3 invalid: the Newton update of the trigamma inverse in `squeeze_var` had the wrong sign, the prior df diverged to ~5·10¹⁴, every gene was tested against one pooled variance (14,311 / 11,209 / 1,135 genes at q < 0.05). E1 and E2 do not use the moderated variance. GSE6364 E1 used 20 samples (see the second defect): θ −0.24 [−0.37, −0.12], p 0.0008 (HECA) and −0.23 [−0.33, −0.13], p 0.0002 (discovery), label `attenuation_detected`; across series `not_replicated`. |
| First correction | [`bulk_replication_correction_freeze.json`](../receipts/bulk_replication_correction_freeze.json) (dd640c5); results [`bulk_replication_correction_results.json`](../receipts/bulk_replication_correction_results.json) | [`bulk_replication_correction/`](bulk_replication_correction/) | corrected trigamma inverse, unit-tested | E3 as tabulated above; E1 and E2 reproduced the first run exactly. |
| Second correction | [`bulk_replication_correction2_freeze.json`](../receipts/bulk_replication_correction2_freeze.json) (a606ade); results [`bulk_replication_correction2_results.json`](../receipts/bulk_replication_correction2_results.json) | [`bulk_replication_correction2/`](bulk_replication_correction2/) | GSE6364 phase regular expression accepts a space (`mid[- ]?secretory`); a unit test checks the configured expressions against every row of the committed GSE6364 sample manifest | The frozen expression `mid-?secretory` did not match the GEO text "Mid secretory", so the 17 mid-secretory samples (9 / 8) listed in the protocol's cohort table were silently dropped and the first two runs' GSE6364 secretory arm held 9 early-secretory samples (6 / 3). With all 37 samples, GSE6364 θ is −0.15 [−0.34, +0.04], p 0.13, label `no_detectable_difference`. GSE51981 (all estimands) and GSE6364 E2 and E3 are byte-identical to the first correction. |

The GSE6364 estimate therefore moved from a detected attenuation on 9 secretory samples to an
undetected one on 26, in the same direction. The first-run label rested on a parsing defect,
not on a prespecified sub-phase analysis, and is superseded; the direction and the
early-secretory sub-phase estimate are reported above as what they are, a sensitivity in a
small series.

## What this establishes and does not

- At tissue level, in the largest public natural-cycle cohort, the stromal secretory program
  responds to disease by less than 40% attenuation and less than 7% enhancement of its
  proliferative-to-secretory shift (HECA program; 35% and 12% for the discovery program). The
  single-cell bound (0.31–0.42 program units on 13 and 20 secretory donors) is sharpened in
  tissue units, with the estimand changed from stromal lineage to whole tissue.
- The frozen single-cell top-20 does not transfer in direction to whole tissue in either
  series or contrast. Whole-tissue expression mixes compartments and reflects composition as
  well as within-cell change, and GCRMA microarray values are not count-based estimates; a
  non-transfer here does not contradict the single-cell null and does not test it.
- The 9,593-gene contrast against all non-endometriosis samples is a property of the
  no-pathology control group of GSE51981, mirrors the series' published pathology-versus-normal
  pattern, vanishes against other-pathology controls and is uncorrelated with the single-cell
  disease statistics. It is reported as context; it nominates nothing.
- GSE141549 was downloaded but not analyzed under the frozen rules; GSE25628 offers no eutopic
  disease contrast. Disease ascertainment, control definitions, sample overlap between GSE6364
  and GSE51981 (both from the same group; overlap not documented in GEO) and the absence of age
  in GEO for both series limit the comparison. Nothing here nominates a target.

## Files

- Protocol: [`BULK_REPLICATION_PROTOCOL.md`](../docs/BULK_REPLICATION_PROTOCOL.md); configurations
  [`bulk_replication.json`](../config/bulk_replication.json),
  [`bulk_replication_correction.json`](../config/bulk_replication_correction.json),
  [`bulk_replication_correction2.json`](../config/bulk_replication_correction2.json) (each
  correction documents its defect under `correction_of`).
- Metadata audit before the freeze: [`bulk_replication_feasibility.json`](bulk_replication_feasibility.json);
  sample manifests [`GSE51981`](../manifests/bulk_replication_GSE51981_samples.tsv),
  [`GSE6364`](../manifests/bulk_replication_GSE6364_samples.tsv),
  [`GSE141549`](../manifests/bulk_replication_GSE141549_samples.tsv).
- Code: [`bulk_replication_metadata_audit.py`](../scripts/bulk_replication_metadata_audit.py),
  [`bulk_replication_pin.py`](../scripts/bulk_replication_pin.py),
  [`bulk_replication_fetch.py`](../scripts/bulk_replication_fetch.py),
  [`bulk_replication_analysis.py`](../scripts/bulk_replication_analysis.py),
  [`followon_freeze.py`](../scripts/followon_freeze.py); tests
  [`test_bulk_replication.py`](../tests/test_bulk_replication.py) (planted attenuation and full
  transfer recovered, null fixture, receipt and freeze enforcement, parsers and moderation,
  GSE6364 regexes against the manifest).
- Receipts: data [`bulk_replication_data.json`](../receipts/bulk_replication_data.json) (six files,
  URL, HTTP status, bytes, SHA-256, UTC); freeze and results receipts per run as tabulated above.
- Outputs per run: `summary.json`, `samples_used.tsv`, per-gene tables and per-sample program
  scores in [`bulk_replication/`](bulk_replication/), [`bulk_replication_correction/`](bulk_replication_correction/)
  and [`bulk_replication_correction2/`](bulk_replication_correction2/); the first run's folder also
  holds the descriptive [`posthoc_unmoderated_check.json`](bulk_replication/posthoc_unmoderated_check.json).
- Raw series matrices, annotation files and the GSE141549 workbooks are not committed.
