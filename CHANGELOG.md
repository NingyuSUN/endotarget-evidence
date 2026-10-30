# Changelog

Versions are Git tags with GitHub releases. Scientific freezes and their receipts
are described in the [method card](docs/METHOD_CARD.md) and the
[publication history](docs/PUBLICATION_HISTORY.md); this file lists what each
tagged version contains.

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
