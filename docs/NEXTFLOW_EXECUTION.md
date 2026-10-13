# Nextflow discovery execution — SRV01

This DSL2 workflow runs the adopted discovery model, the three frozen rankings and matched null-set generation. It does not run HECA effects, bootstrap, leave-one-donor-out, dossiers or reporting. Completing discovery still requires the separate aggregate discovery freeze before external effects.

## Inputs and safety gates

Production defaults are the accepted discovery counts, clinical ledger, numeric design, STRING gene-edge reference and statistical/ranking configurations. No production task starts without `--pre_effect_freeze` in phase `STATISTICAL_IMPLEMENTATION_FROZEN_PRE_EFFECT`.

The freeze must hash all selected inputs plus `main.nf`, `nextflow.config`, all three stage scripts, `donor_de_lib.R`, runtime wrapper, explicit runtime lock/YAML and runtime installation/R receipts. Additional artifacts in the manifest are also verified. Every path must be project-relative, normalized and resolve inside this project; arbitrary production-count substitutions are rejected. Each process checks the same freeze digest and every artifact again immediately before its command.

The real implementation freeze is created by the root integration step after code review and local commits. A synthetic manifest never authorizes production: it has phase `SYNTHETIC_TEST_ONLY_NOT_A_SCIENTIFIC_FREEZE`, requires `--synthetic_test true`, and all artificial biological inputs must remain under `work/synthetic/`. The R CLI separately enforces resolved synthetic paths.

## Process connections and outputs

1. `DONOR_DE`: donor-level edgeR normalization and limma-voom fitting, with exact frozen statistical settings. Retains the complete expression-pass family in `de/de_all.tsv`; finite estimable rows go to `de/de_fitted.tsv`.
2. `RANK_TARGETS`: uses both DE tables and the same discovery-induced network for B0 expression, B1 weighted degree and G1 propagation. Failed fits on the comparative universe stop ranking. Outputs `ranking/`.
3. `MATCHED_NULLS`: uses the saved discovery expression measure, weighted degree and each method's own top-20 exclusions to generate 200 matched sets per estimable method. Scarce bins produce method-level `NOT_ESTIMABLE`; no merging/rescue. Outputs `nulls/`.

All three stage directories are copied to the chosen output directory, preserving stage provenance. Task files remain under the configured Nextflow work directory. New output directories are mandatory even on `-resume`, so earlier copies cannot silently mask changed results.

## Actual native command

Run only on SRV01 under the project root. The profile is named `srv01` because `native` is a Groovy reserved word.

```bash
cd /home/analyst/te
source scripts/statistical_runtime_env.sh
nextflow -log logs/discovery_primary.nextflow.log run main.nf -profile srv01 \
  --pre_effect_freeze receipts/statistical_pre_effect_freeze.json \
  --outdir /home/analyst/te/work/statistical/discovery_primary \
  -work-dir /home/analyst/te/work/nextflow \
  -with-trace logs/discovery_primary.trace.tsv
```

This is an execution recipe, not evidence that production has run. The parent integration step controls when the prerequisite freeze exists and execution is authorized. `-work-dir` is a Nextflow option with one dash; `--work-dir` would be an unrelated workflow parameter.

The isolated runtime currently records Nextflow 25.10.2, Java 21 and R 4.4.3; exact builds/package versions are in the committed runtime receipts/explicit lock. The workflow pins Nextflow 25.10.2. It uses one CPU per process, a 4 GB declared memory budget, one local executor slot and terminating errors. Local-executor memory is a scheduling declaration, not an OS-enforced RSS limit. No dependency installation or downloads occur inside tasks.

## Cache validity and restart

`cache = 'deep'` hashes staged content. The manifest itself, its digest and every manifest artifact are process inputs, including source/config/runtime files. Code or input changes require a new authorized freeze and cannot reuse the old task cache. Upstream stage directories also participate in downstream cache hashes.

A valid unchanged freeze can use `-resume`; choose a new publication directory. A changed source with an unchanged manifest stops at the hash guard before launching tasks. The original work/cache/output directories remain preserved.

## Synthetic verification

```bash
cd /home/analyst/te
bash tests/smoke_nextflow.sh
```

The fixture generator creates 1,500 artificial genes, 20 artificial donors (10 per group), a connected synthetic association graph and a separately labeled test freeze. It reads configuration/code hashes but never opens actual biological count files.

The smoke runs all three actual scripts through Nextflow, then checks 3/3 cached processes and byte-identical outputs on `-resume`. Negative gates reject modified frozen data, missing freezes, production use of a synthetic freeze and missing required code hashes. A receipt is saved as `receipts/nextflow_synthetic_<UTC>.json`; artificial inputs, outputs, logs and task work remain under `work/synthetic/<id>/`.

`bash tests/smoke_nextflow.sh --verify-only /home/analyst/te/work/synthetic/<id>` rechecks preserved outputs and the omitted-code negative case without recomputing the stages. This is useful only when the pipeline artifacts have not changed; code changes require a fresh fixture/run.

The first full chain and cached rerun succeeded, but its final verifier initially expected the wrong synthetic label. The harness was corrected to the R provenance value `synthetic_test_only`; the preserved outputs were verified without recomputation. A later R path/weight/residual-df hardening and per-task guards justified a new complete smoke. The relevant receipts retain these boundaries.

## Docker design — not executed

The `docker` profile requires an immutable `--container_image IMAGE@sha256:DIGEST`. It mounts the project root read-only and the active task work directory read/write. The image must itself contain the pinned R/Python packages; the native Conda wrapper is not sourced inside containers.

Docker daemon access is currently denied to this account. No container image was built or run, and mount behavior is unverified on the daemon. A native smoke cannot be reported as container validation. AWS remains unexecuted.

## Interpretation

These checks establish workflow connections, frozen-input enforcement and cache mechanics on synthetic data. They do not establish real-cohort model suitability, cross-study replication, target accuracy, therapeutic actionability, causal effects or clean-container execution. Each null method's status must be inspected; successful process exit alone does not mean all null bins were estimable.
