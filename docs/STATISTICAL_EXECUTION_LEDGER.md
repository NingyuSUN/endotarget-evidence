# Statistical execution ledger

Scope: implement the approved donor-level statistics and ranking workflow on SRV01 under /home/analyst/te. Entry base: 9252a54. Commit identifiers in this ledger are historical execution identifiers; see the [publication history](PUBLICATION_HISTORY.md) for how they map to published commits.

- Input gate: M0_INPUT_FEASIBILITY_PASS_WITH_DOCUMENTED_LIMITS; no disease effects existed at entry.
- AI-assisted roles: gpt-6-astra (high) drafted the pre-effect method contract; gpt-6.1-sol prepared the isolated runtime (medium) and the synthetic network ranker (high); the lead session performed the statistical implementation and integration.
- Ruling: reuse the existing approved numeric primary design matrices verbatim, verifying exact donor order and disease column. This avoids a second implementation of categorical reference levels. Cost if incorrect: source design errors persist; independent full-rank, label and coefficient checks remain mandatory.
- Synthetic tests use invented counts and persons, never evidence of biological performance.
- No external disease fit before discovery/config/code/software/rank freeze. No real container execution is claimed while daemon account access is unavailable.

- Runtime committed 39f1053; exact R 4.4.3 / edgeR 4.4.0 / limma 3.62.1 / statmod 1.5.2 / jsonlite 2.0.0; isolated 117-package environment, no Docker access.
- Statistical contract committed 1515ba5; ranker 62504ea; null sampler 62938ad.
- Root R model: 24 synthetic checks PASS, including independent weighted-lm coefficients, sign reversal, group-preserving filter invariance, cell/donor/count/weight/rank guards. Failed initial missing implementation and one fixture syntax issue preserved in logs.
- Final bounded gpt-6-astra source review identified no current primary estimator blocker. Installed limma supports both adaptive.span and legacy; explicit FALSE/TRUE used. Additional finite-positive precision and per-gene df guards implemented and tested.
- Freeze module: six boundary tests PASS. It verifies real-versus-synthetic provenance, common U, exact top-list prefixes, null upstream hashes, explicit null failure, full-U directions, and immutable input hashes.
- Pipeline integration and actual discovery execution remain pending at this ledger entry. Do not infer scientific success from software checks.

- Final integration: 96 Python tests PASS in m0-h5. Statistical-only runtime lacked requests for the legacy ingestion test module; recorded failure, then used the existing correct data-preparation environment without mutating either lock. R 24 checks and actual three-stage Nextflow synthetic execution ran in the statistical runtime.
- Nextflow committed 34660b4: real synthetic 3/3 completed, 3/3 cached on resume, byte-identical outputs, four preflight negative gates. Each method 200 × 20 matched sets. Final receipt nextflow_synthetic_20261013T100618Z.json. Container remains unexecuted at this ledger entry.
- Review all current artifacts, commit implementation, and create pre-effect freeze before first real discovery invocation.

- Actual discovery ran after the 7ad353f pre-effect freeze: 3 Nextflow processes completed, no warnings, 17,733 tests, 0 BH q < 0.05. Same U: 11,749 ranks and 600 matched sets created. No outcome-driven parameter adjustment.
- Complete discovery/null freeze committed ebb5652 BEFORE any HECA effects. Negative/inconclusive interpretation preserved in reports/discovery_primary_results.md. External evaluation followed.
