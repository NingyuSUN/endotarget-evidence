# Replication execution and interpretation

Status before public external effects: implementation under final pre-effect review.
Source of scientific decisions: the frozen RESEARCH_PROTOCOL, PRE_EFFECT_METHOD_DECISIONS,
STATISTICAL_IMPLEMENTATION_DECISIONS and config/statistical_analysis.json.
Historical status wording inside those frozen documents is retained; docs/CURRENT_STATE.md
records later milestones without altering prior checkpoints.

## Native SRV01 workflow

All commands run from /home/analyst/te on SRV01, with scripts/statistical_runtime_env.sh sourced.
The two-process replication.nf reads the original discovery freeze and a second implementation
freeze made before any external disease effects. It fits the separately normalized HECA donors,
then applies the frozen discovery lists and precomputed matched null memberships.
No HECA effect enters discovery selection.

Planned first real command:

    source scripts/statistical_runtime_env.sh
    nextflow -C replication.config -log logs/heca_primary.log run replication.nf -profile srv01 -params-file config/replication_run.json -with-trace logs/heca_primary_trace.tsv

config/replication_run.json will point to receipts/replication_implementation_freeze.json
and work/statistical/heca_primary. It is a launch record, not a new statistical choice.
A new output directory is required. To use Nextflow cache, select another empty output
directory and add -resume; never overwrite the first published output.

## Evaluation

B0 is the expression-only ranking, B1 the network-degree-only comparator, and G1 the fixed
network-propagation ranking. Each method retains its original top20 and secondary top50.
Every gene is classified as evaluable or by one mutually exclusive missing/zero/invalid reason.
The reported count is agreeing / evaluable / originally selected. Zero external effects remain
unevaluable for direction but retain a finite signed t of zero when available.
Missingness bounds and >=80% coverage gate any possible comparative claim.

External multiple-testing adjustment uses the complete HECA testing family. Direction agreement
can reflect tiny noisy changes. A positive direction count is neither target accuracy nor
drug efficacy. Matched-set percentiles are descriptive, not empirical significance p-values.

The fixed 40-gene dossier pool is the union of B0/G1 top20. A strongly_conflicting category
requires opposite direction and a pointwise, unadjusted 95% external interval excluding zero.
That label is not a claim of multiplicity-adjusted conflict; always retain the full-family q.
Do not force any supported dossier if no same-direction q<0.05 gene exists.

## Resampling, after first external evaluation

    Rscript scripts/resample_donor_models.R --mode bootstrap --output-dir work/statistical/bootstrap --workers 4
    python scripts/run_discovery_lodo.py --output-dir work/statistical/discovery_lodo

All 500 within-disease donor draws are saved before any resampled fit. The same draw applies
to every method. Gene family, discovery lists/directions and original external evaluability
are fixed. Draws that fail fitting are retained as failures; >50 failures suppress the paired
interval. Resampled donor duplicates are not new people. Only constant nuisance terms and
unused factor levels may be dropped; nonconstant confounding is not repaired.

LODO removes one discovery person at a time, repeats filtering, model fitting and graph ranking,
and reports overlap/20 plus intersection/union. Failed folds stay failed. Both procedures bind
checkpoints to exact input/code hashes. Their uncertainty scopes differ: external bootstrap
conditions on discovery choices; LODO describes discovery rank sensitivity.

Four bootstrap workers and one serial LODO worker use one numerical thread per worker.
A fresh host read showed 36 CPUs, approximately 29 GiB available memory and ~60 GiB disk free.
Preserve at least 20 GiB free and reduce concurrency if actual measured resource use warrants it.
Do not restart completed analysis solely to refresh a report.

## Verification and remaining execution limits

An actual two-stage artificial-data Nextflow run completed. Its cached resume and negative gates
are recorded in receipts/replication_nextflow_tests.json. Component tests and resampling
artificial runs have separate receipts. These validate software behavior, not real replication.

The docker profile requires an immutable image digest and binds the project read-only except
task output. It is an unexecuted design until a real image build and container run succeed.
The installed Docker client is currently inaccessible to analyst at the daemon socket.
No host privilege change is authorized. Native execution must never be reported as container
validation. AWS is an optional, unexecuted deployment path.

## Learning target

YELLOW — understanding not yet assessed. Explain donor-level replication, fixed evaluation
denominators, why no list refill is allowed, full-family multiple testing, the two uncertainty
scopes, and observational association versus intervention causality.
