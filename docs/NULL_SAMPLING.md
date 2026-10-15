# Discovery matched null memberships

Implemented on SRV01 `/home/analyst/te`, tested only on wholly synthetic fixtures. No actual disease effects, target lists, real null memberships, external outcomes or scientific validation were produced in this implementation round. Adopted rules: `docs/STATISTICAL_IMPLEMENTATION_DECISIONS.md` section 6 and the exact `null` object of `config/statistical_analysis.json`.

## Question and interpretation

These sets ask whether an observed fixed-list result differs descriptively from gene sets with the same coarse expression/degree structure. They are not 200 biological controls, independent studies, permutations of disease labels, a calibrated null of no disease effect, or proof of causal specificity. Sets can reuse genes across draws, share donors and use disease-informed discovery expression directions. External evaluation and its missingness/coverage rules are a later stage; a sampled gene that is missing externally is never replaced.

## Inputs and execution

```bash
cd /home/analyst/te
source scripts/statistical_runtime_env.sh
python scripts/sample_matched_nulls.py \
  --de DISCOVERY_MODEL_DIR/de_fitted.tsv \
  --ranking-dir DISCOVERY_RANKING_DIR \
  --config config/statistical_analysis.json \
  --output-dir NEW_NULL_DIRECTORY
```

The discovery fitted table must be the same byte-hashed input used for ranking, with unique unversioned `ENSG` IDs, finite `logFC` and `null_matching_logCPM`. The expression field is computed upstream as donor-mean edgeR normalized log CPM with prior count 0.5; it is not `AveExpr` or the average voom expression. The sampler does not recompute counts/normalization or inspect external data. It validates the frozen rank provenance and recorded output hashes, sorted unique U of >=50 genes, full B0/B1/G1 tables on exactly U, rank/method coding, each top20's exact prefix, saved directions and weighted degrees. It recomputes weighted degree from saved induced edges and verifies agreement, retaining the original ranking float64 value for binning to avoid summation-order changes. Self loops, duplicate undirected edges, endpoints outside U, invalid edge weights, source changes, zero/nonfinite degrees and malformed TSVs fail instead of shrinking U.

`--synthetic-test` requires synthetic rank provenance and explicitly labels every result `SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE`. Production rejects synthetic ranks. There is no external-data interface and no seed, set-size, pool-expansion, quantile or retry override. The root orchestration owns the pre-effect implementation freeze, the complete pre-external freeze and authorization to calculate real cohort effects; successful sampling alone does not establish those gates.

## Exact bin and sampling contract

Take pooled U's expression and frozen weighted degree separately. Calculate internal quintile cutpoints at 0.2/0.4/0.6/0.8 with `numpy.quantile(method="linear")`, equivalent to R type 7. Collapse duplicated internal cutpoints with `numpy.unique`. Assign `1 + numpy.searchsorted(cutpoints, value, side="left")`, so equality belongs to the lower interval. Identical values are not split by stable gene-ID rank. Constant vectors may leave empty intervals; no rescue is attempted. Expression/degree crossing has at most 25 possible bin tuples.

For fixed method order B0/B1/G1, count its own top20's exact joint-bin occupancy. Candidate pool is U excluding only that method's top20; other methods' selected genes remain eligible. Stable-ID-sort candidates within bins. Before any draws, check all required-bin demand against available candidates. Any shortage makes the whole method `NOT_ESTIMABLE`, reason `INSUFFICIENT_OWN_EXCLUDED_JOINT_BIN_POOL`, with 0 completed sets and an empty membership table. Bin receipts show demand/availability even for empty insufficient pools. No bin merge, tolerance extension, replacement sampling or partial-set output occurs. Other methods retain their own independent estimability checks.

For each estimable method separately initialize `numpy.random.Generator(numpy.random.PCG64(42))`. Generate exactly 200 draws, traverse demanded bin tuples lexicographically, and call `choice` uniformly without replacement in each bin with its exact demand. Each resulting set contains 20 unique genes. Genes may recur between sets. Fresh method streams are reproducible but can be correlated; they are not independent scientific replicates. Each sampled gene retains its own exact discovery `logFC` sign `-1/0/+1`, including zero and arbitrarily small nonzero coefficients. Signs, t statistics, p values and external availability never form bins or choose memberships.

## Output artifacts and preservation

All output files are written into a new sibling `.partial` directory and renamed only after complete tables/status are written. Existing final or partial directories cause failure and are preserved. All source/config/implementation files are rehashed before output creation.

| File | Meaning |
|---|---|
| `universe_bins.tsv` | Every U gene's expression/degree, joint bins, coefficient and exact direction |
| `all_U_directions.tsv` | All U coefficients/directions for the aggregate pre-external freeze |
| `boundaries.json` | Exact float64 internal cutpoints, quantile method and equality convention |
| `bin_demands.tsv` | Each method/bin's exact requested and available counts and shortage status |
| `null_memberships.tsv` | Combined B0/B1/G1 actual members, method/draw/position, bins/covariates/signs |
| `B0_null_membership.tsv`, `B1_null_membership.tsv`, `G1_null_membership.tsv` | Same memberships separated per method; failed method has header only |
| `status.json` | Method estimability/reason, 200 requested/actual draws, own excluded genes, source/output hashes, actual software and timestamp |

`status.json` contains output hashes for every other artifact, but cannot hash itself. The root aggregate freeze must hash `status.json` together with all ranks, U, graph, bin boundaries, directions, actual memberships, donor/design/software/config inputs and evaluation rules. The status always declares `complete_external_freeze=false` and `scientific_validation=false`.

## Synthetic verification

```bash
source scripts/statistical_runtime_env.sh
python -m unittest discover -s tests -p test_matched_nulls.py -v
```

The executable tests cover independently hand-calculated type-7 boundaries, duplicate-cutpoint collapse and lower-boundary assignment; exact joint rather than just marginal occupancies; 200 sets of 20 unique genes and own-only exclusions; exact zero/tiny-sign retention; independently computed PCG64 first draw and stable-order reproducibility; reinitialization per method; invariance of membership to signs/t/p/external fields; complete-method failure on one scarce bin without rescue; malformed/nonfinite/duplicate inputs; CLI generated files/hash identity and preservation; changed discovery/rank files; synthetic production rejection and config drift. Tests use temporary fixtures under remote `work/` and remove them afterward. Receipt: `receipts/null_sampling_tests.json`.

An initial boundary assertion used literal 0.2; NumPy's type-7 interpolation produced 0.20000000000000018. The test assertion was corrected with a 1e-15 absolute floating-point comparison. The stored cutpoint and binning logic were not rounded or changed. All subsequent assertions pass; numerical agreement is software evidence, not external biological validation.

Interpretation notes: joint-bin occupancy determines which genes can be matched; identical values remain together; a scarce pool is not permission to rescue the null experiment; and these 200 sets are descriptive sensitivity evidence rather than 200 biological negative controls.
