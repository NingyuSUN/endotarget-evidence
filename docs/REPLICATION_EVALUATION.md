# External replication evaluation — pre-result implementation

Status: IMPLEMENTED_AND_SYNTHETICALLY_TESTED; no real HECA effects read or evaluated by this implementation round. Original discovery/statistical freezes and their artifacts remain unchanged. This addition consumes their saved predictions and never revises a discovery rank, sign, universe, or random membership.

## Interface and source verification

Production command, from /home/analyst/te:

```bash
envs/m0-h5/bin/python scripts/evaluate_replication.py \
  --discovery-freeze receipts/discovery_rank_freeze.json \
  --external-dir work/statistical/heca_primary/de \
  --output-dir work/statistical/heca_primary/evaluation
```

Production requires DISCOVERY_FROZEN_PRE_EXTERNAL, external_disease_effects_computed=false in the discovery checkpoint, public_cohort HECA provenance and scientific_effect_run=true. Verify every original freeze artifact first. Verify all three full rankings, exact top20/top50 prefixes, all-U saved coefficients/signs, exactly 200 original sets per estimable null method, null members' own signs, nonselected pool membership and no within-set duplicates. No new random sets are generated.

External de_all.tsv, de_fitted.tsv, expression_filter.tsv and provenance.json must exist. Check external input signatures against the original discovery checkpoint, including frozen count/design/clinical/config/model code, original implementation freeze and exact discovery-freeze file hash. Every external output in provenance is rehashed at its published basename, including the model/session artifacts. Verify the external all-input-gene filter ledger has the exact original counts IDs and order, and that count dimensions agree with provenance. This distinguishes absence from expression failure.

The de_all family must equal the expression-pass IDs. de_fitted must equal its valid finite-effect/t/p subset without changing rows. Recompute BH as a validation check over valid p-values using n equal to ALL expression-pass genes, including failed fits. Preserve original p/q/log2 effect/moderated-t/95% interval; the check does not substitute subset-adjusted values. Failed fits legitimately retain NA q; zero or finite q on failed fits is rejected.

Synthetic mode adds --synthetic-test, requires SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE and synthetic_test_only HECA provenance, and restricts freeze/external/output paths to work/synthetic/. Synthetic cannot enter production. Its provenance signatures are checked against the fixture freeze; real count-ID equality is validated in production. This is a distinct test mode, never relabelled public data.

## Fixed estimand and missingness

Each selected gene gets exactly one exclusion reason, in this precedence: absent_id, expression_filter_fail, invalid_fit, discovery_zero_sign, external_zero_sign. Valid direction means nonzero saved discovery coefficient sign and nonzero finite external coefficient/t/p with fit_estimable=true. Zero is unevaluable, never opposite; tiny finite nonzero effects have their actual sign, with no numerical epsilon.

Report agree/evaluable/selected separately and never refill. Concordance is agree/evaluable, NA if none evaluable. Missing bounds are agree/k through (agree+unevaluable)/k. G1-minus-B0 subtracts the opposite missing-bound endpoints. At least 16/20 direction-evaluable genes are required in both primary lists, and zero-containing or negative contrast bounds cannot support superiority. Even favorable point/bounds leave superiority_claim=false until the separately planned paired bootstrap is available. Top50 is secondary. Per-gene signed_external_t is the original valid-fit external moderated t multiplied by that gene's nonzero saved discovery sign. A valid external exact-zero effect retains signed t=0 while remaining direction-unevaluable; absent, failed, or discovery-zero predictions retain NA signed t.

## Null summaries

Evaluate all original 200 sets per estimable method using exactly the same exclusion precedence and each member's own frozen discovery sign. Save selected/evaluable/agree/opposite, each missing category, coverage and bounds for every set. Coverage adequate means at least 16/20. The observed descriptive percentile is 100 times the fraction of coverage-adequate sets whose concordance is <= observed, including ties. Save total/adequate/inadequate counts and the complete adequate distribution. No adequate sets or no evaluable observed list yields NA. Observed-list coverage remains explicitly available in summary.json; a percentile never bypasses its primary coverage guard. p_value is always null. Gene dependence, shared donors and discovery-informed directions preclude a calibrated permutation p-value.

## Predefined dossier pool and selection

The pool is the fixed union of B0 and G1 top20; in the actual frozen discovery lists it contains 40 IDs. The full pool is always written, including absence, filter failure, invalid fit and zero signs. Deterministic ordering is ascending full-discovery B0 rank, then stable ID, for every category:

1. Up to three supported_same_direction: direction-evaluable, same sign and external full-family BH q<0.05.
2. Up to one strongly_conflicting: direction-evaluable, opposite sign and the external **unadjusted** 95% interval excludes zero (CI_low>0 or CI_high<0). Report BH q separately. A q>=0.05 conflict is not multiplicity-adjusted significant evidence.
3. Up to one inconclusive: everything outside those two qualified classes. Its exact missing/zero reason remains visible.

Choose only qualifying genes. Do not force positive/conflicting categories or search beyond the fixed pool. No genetic annotations enter this transcriptomic selection. Dossier tables are selection preparation, not completed literature/genetic dossiers or therapeutic claims.

## Outputs and bootstrap handoff

- selected_genes.tsv: 210 rows = 3 methods times (20+50), no refill, original effects/p/q/intervals and exclusions.
- bootstrap_interface.tsv: identical LONG table with method, k, rank, gene_id, frozen discovery_sign and original_evaluable TRUE/FALSE. Resampling deduplicates the union of gene IDs while retaining method/k membership and fixed original masks. Originally missing genes never become newly evaluable.
- null_sets.tsv and null_summary.json: all frozen random-set coverage/missing bounds and descriptive distributions.
- summary.json: per-list counts, primary contrast guards and selection counts.
- dossier_pool.tsv and dossier_selection.tsv: complete fixed pool and deterministic qualified choices.
- provenance.json: exact code/freeze/external-source hashes, validated external signatures and output hashes.

Existing nonempty outputs are preserved and rejected for overwrite. All validation runs before writing evaluation artifacts.

## Actual synthetic validation and limits

Manual golden cases cover 10 agreeing + 5 opposite + one of each five exclusion reasons (15/20 evaluable, bounds 0.50–0.75), exact/tiny-zero behavior, missing-reason precedence, no refill, complete-family BH including failed fits, zero-q rejection, null coverage, zero-spanning contrast bounds, own random-gene signs and deterministic dossier eligibility. An artificial 90-gene CLI fixture checks all three full/prefix lists, 210 output rows, 600 membership sets, a 40-gene pool, 16/20 B0 coverage and no forced supported/conflicting dossiers. Negative file tests reject output changes, mismatched original input signatures, synthetic freeze in production, and wrong null-member signs.

These are engineering checks on invented data. They do not establish real replication, meaningful biological effects, sampling uncertainty, causal genes, treatment targets, independent annotation truth, actual Docker execution or AWS execution. The separate gated Nextflow HECA-model plus evaluator synthetic smoke remains the integrator's responsibility before real effects. A paired 500-draw external bootstrap remains a separate stage.

YELLOW — understanding not yet assessed. The user should be able to explain why a missing or zero-effect gene is not an opposite-direction observation, why BH uses the entire external expression family, why 200 matched sets are not 200 biological controls, and why an unadjusted interval conflict is different from full-family significant evidence.
