# Prespecified clinical sensitivity execution

This new implementation executes adopted clinical subsets after both primary discovery and external results have been observed. It is not a preregistration,blind replication or a revision of the frozen primary analysis. All definitions,eligibility flags and numeric designs already existed; `work/statistical/clinical_sensitivities/input_snapshot.json` records their exact old bytes,source hashes,design last Git commits and selected donor order before the first sensitivity fit. Original primary/bootstrap/LODO scripts,configuration,counts,clinical ledgers,lists and outputs are preserved.

Run from `/home/analyst/te` with one process and one BLAS/OMP thread:

```bash
source scripts/statistical_runtime_env.sh
Rscript tests/test_clinical_sensitivities.R
python -m unittest discover -s tests -p test_clinical_sensitivity_reporting.py -v
Rscript scripts/run_clinical_sensitivities.R
python scripts/summarize_clinical_sensitivities.py
```

The final two commands have already run on real data. Do not repeat successful models. The wrapper refuses a nonempty output tree; reporting refuses existing report files and execution receipt. Fixes requiring new computation must preserve old artifacts and document the versioned reason.

## Exact scenarios

|Scenario|Discovery|HECA|
|---|---|---|
|Strict histology-conflict exclusion|20people11ENDO/9non-ENDO,rank12,df8;fitted|13people10/3;NOT_ESTIMABLE,no fit|
|Adenomyosis-free|18people8/10,rank11,df7;fitted|19people12/7,rank9,df10;fitted|
|Strict plus adenomyosis-free|16people7/9;NOT_ESTIMABLE,no fit|No extra combined scenario|

Subset the unchanged count matrix and primary clinical ledger to each existing design donor_id in exact order. Require exact membership equality with discovery strict_histology_sensitivity/adenomyosis_free_sensitivity and HECA histology_conflict/adenomyosis fields. Recheck group/cell minima,mixed batch support,rank,df,age centering and disease coding. For true viable subsets,use unchanged `fit_donor_model` from frozen donor_de_lib.R: rerun count filtering,retained-library sums,TMM,voom and full-family robust eBayes/BH. Adenomyosis columns absent from the existing constant-zero subset designs are not reintroduced. Failed minima are explicit outcomes; no covariate removal,threshold relaxation or alternate fitting is allowed.

## What the outputs answer

`reports/clinical_sensitivity_scenarios.tsv` lists all five outcomes,people,rank/df,filter families and BH counts. Discovery paired primary-gene effect tables reside under each executed discovery sensitivity directory. `reports/clinical_sensitivity_discovery_robustness.tsv` lists all-primary-family and frozen B0/B1/G1 top20/top50 direction agreement with explicit valid/missing denominators. This measures the stability of saved primary effects under clinical subsetting; no new rankings or predictions replace the original freeze. Missing/filter-failed/invalid/exact-zero genes remain separate from opposite effects. The unadjusted CI conflict flag is not multiplicity-adjusted significance.

For HECA adenomyosis-free,reuse unchanged evaluator classify/evaluate_set rules against frozen original discovery lists/signs. The separate fixed-list details and summaries include selected/evaluable/agree,coverage,missing bounds and BH-supported same-direction counts. They do not change original external bootstrap masks or imply a new bootstrap comparison. Observational associations and direction concordance do not prove mechanism or therapeutic action.

## Ascertainment and equivalent sensitivities

`clinical_sensitivity_ascertainment_donors.tsv` exports44selected people with unmodified source age,cycle day,progesterone where cached,histology,expression-assisted phase where separately recorded,author phase,hormone original fields/source-qualified eligibility and strict/free membership. Discovery original day/progesterone/histology are joined by sample_id to discovery_clinical_candidates.tsv. HECA values come from its frozen selected ledger. Discovery has no separately recorded per-person expression-phase field in this cached table; do not manufacture one from the final author label. Missing HECA hormone NA/nan remains missing,not FALSE. Groupwise missing/numeric counts and original ranges are in the ascertainment summary.

Strict means the predeclared exclusion of explicit histology conflicts. Missing histology is still missing,so this is not universal independent histology confirmation. The observed day28control stays visible without new exclusions. No post-result cycle-day boundary or imputation is introduced. All44ages are complete; no new age-unadjusted fit is warranted by the missing-age condition. The30/100cell thresholds retain identical original cohorts/libraries and therefore reuse primary results rather than supply independent refits.

## Verification and limits

Focused synthetic tests establish exact subset order,eligible membership,failed7/9minimum,no donor-order swaps,no invented phase values,unknown versus numeric day preservation,retaining day28,and missing genes not becoming opposite effects. Artificial end-to-end execution fits the three viable scenarios and records both failed cases without DE/RDS. The real receipt contains source/output hashes,the actual single-process PID/exit status and all scenario summaries. These engineering checks do not validate cycle-phase truth,remove clinical confounding,demonstrate power or turn0BH-significant genes into proof of no biological difference. Reduced subsets retain annotation/ascertainment dependence and residual pathology/batch uncertainty.

YELLOW — understanding not assessed: distinguish a preset subset robustness check from a newly selected model,explicit conflict exclusion from histological confirmation,missingness from negative evidence,and complete age availability from elimination of age confounding.
