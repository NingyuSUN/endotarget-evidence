# Changelog

Versions are Git tags with GitHub releases. Scientific freezes and their receipts
are described in the [method card](docs/METHOD_CARD.md) and the
[publication history](docs/PUBLICATION_HISTORY.md); this file lists what each
tagged version contains.

## Unreleased

- **Published bulk signatures in stroma (2026-11-02).** Four fixed whole-tissue disease
  signatures (Poli-Neto 2020 proliferative stage I–II and III–IV, Bunis 2022 proliferative
  disease versus control, their consensus), built from publisher supplements downloaded with
  receipts, tested under a frozen protocol in the pooled stromal statistics (directional
  concordance with expression-matched references carrying the same direction labels; unsigned
  enrichment) and as donor-level signature scores on the v1 primary design with MDE. Two sets
  `concordant` (58% and 55% sign agreement versus 48% and 47% matched), none enriched, no
  detectable donor-level score difference (−0.01 [−0.08, +0.07] and +0.08 [−0.06, +0.23] log2
  units for the largest set). Supplement fetch script with receipts; standard-library xlsx reader.
- **Bulk program-level replication (2026-11-02).** The frozen single-cell phase programs
  and B0 top-20 taken into public bulk eutopic cohorts GSE51981 (144 natural-cycle samples)
  and GSE6364 (37) under a protocol frozen before any expression value was read; data
  downloaded on the analysis host with a receipt. Secretory program disease interaction
  θ −0.09 [−0.21, +0.03] and −0.15 [−0.34, +0.04], `no_detectable_difference` in both;
  top-20 `does_not_transfer`; 9,593 tissue-level genes against all controls, 0 against
  other-pathology controls. One run plus two published corrections (variance moderation
  sign error; GSE6364 phase regex that dropped 17 mid-secretory samples), each with its own
  freeze and results receipts; the earlier outputs stand as reported. Metadata audit, fetch,
  pin and analysis scripts with tests.
- **Stromal-state composition test (2026-11-02).** Prespecified per-cohort test of the
  decidualized-like share of primary stromal cells on the v1 design: odds ratio 0.75
  [0.48, 1.17] and 1.08 [0.64, 1.80], no author state at q < 0.05; labels
  `no_detectable_shift`, pair `no_replicated_shift`.
- **Genetic set-level test (2026-11-02).** Four fixed gene sets (Rahmioglu 2023
  nearest and SMR genes; Open Targets 26.09 genetic association > 0 and ≥ 0.5) show no
  enrichment of the pooled stromal disease statistics against expression-matched
  genes. Generic freeze/verify module `scripts/followon_freeze.py`.
- **Two-cohort combined bound (2026-11-02).** Prespecified fixed-effect combination
  of the two saved v1 disease fits over 15,435 shared genes: 0 genes at BH q < 0.05;
  median detectable effect 1.27-fold at α 0.05 and 1.5–1.8-fold at the multiplicity
  faced; twofold differences excluded family-wide for 60% of genes. Protocol, freeze
  and results receipts, per-gene table, figure and report; tests; CI installs the
  v1 runtime's NumPy and SciPy versions.
- Summary figure of the follow-on analyses with receipt and boundary tests; Markdown
  link checker with tests; changelog and citation sections (2026-10-30).

## 1.0.0 — 2026-10-28

First tagged release of the public repository.

- **v1 benchmark (frozen 2026-10-13 to 2026-10-14).** Donor-level
  proliferative-phase stromal disease contrast in GSE266265 (12 / 10 donors) and
  HECA E-MTAB-14039 (15 / 7); 0 genes at BH q < 0.05 in each cohort; frozen
  expression, degree and propagation lists evaluated externally
  (14 / 20 / 20, 3 / 13 / 20, 6 / 19 / 20); no supported target shortlist;
  one inconclusive PARP9 dossier, recommendation hold.
- **Phase positive control (2026-10-26 to 2026-10-27).** The same pipeline on a
  secretory-versus-proliferative contrast: 7,367 and 3,986 genes at q < 0.05;
  the frozen expression top-20 transfers 20 / 20; G1 − B0 −0.37 [−0.53, −0.21].
- **Minimum detectable effect (2026-10-27).** Median gene 1.5-fold at α 0.05,
  2.2–3.4-fold at the multiplicity faced; interval-exclusion counts; donor
  subsampling curves on the phase program.
- **Secretory program disease contrast (2026-10-27).** θ = +0.04 [−0.18, +0.25]
  and −0.11 [−0.40, +0.19] program units; detectable attenuation 0.31 and 0.42;
  pair label `no_detectable_difference`.
- **Cohort extension audit (2026-10-27).** GSE213216 and GSE179640 hold no
  eligible natural-cycle eutopic disease stratum; nothing downloaded.
- **Write-up.** Method card, research brief (two-page PDF re-rendered), English
  results summary and both README homepages integrate the four follow-on analyses.
- **Software.** `CITATION.cff`; GitHub Actions workflow running the synthetic
  unit tests (Python on the pinned probe lock, R on Bioconductor 3.20) and a
  Markdown link check.

Not included: large count matrices, `work/` checkpoints, runtime environments
and caches ([reproduction guide](docs/LOCAL_EXECUTION.md)).
