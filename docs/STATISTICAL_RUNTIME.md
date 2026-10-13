# Statistical/runtime environment — SRV01

VERIFIED: installed /home/analyst/te/envs/statistical-runtime with117 packages. No actual disease effects or target ranks were computed. No biological, Docker, AWS or full-pipeline validation is claimed.

## Versions and reproduction

R4.4.3; edgeR4.4.0; limma3.62.1; statmod1.5.2; jsonlite2.0.0; Python3.10.22; numpy2.2.6; scipy1.15.2; OpenJDK21.0.10-internal; Nextflow25.10.2 build10555.
R4.4/Bioconductor3.20-generation packages were solved together. Latest online releases are not implicitly substituted.

Exact117 versions/builds: requirements/statistical-runtime.yml.
Install URL/MD5 lock: requirements/statistical-runtime-planned.explicit.txt.
Independent installed URL/MD5 export: requirements/statistical-runtime.explicit.txt.
Receipt checks every installed version/build against plan and only conda-forge/bioconda channels. Planned download592,345,435bytes.

Fresh-prefix installer: bash scripts/setup_statistical_runtime.sh. It refuses to mutate an existing conda runtime and requires26GiB free (20GiB reserve+6GiB cache/unpack budget). It uses project cache/conda-pkgs and tmp/, with explicit channel override and no global conda configuration writes. Existing base and m0-h5 are not installation targets. Base conda metadata predates this resolution; m0-h5 is a Python venv whose configuration hash/imports are recorded. No before/after hash baseline was collected; do not claim exhaustive unchanged-file verification.

From /home/analyst/te, source scripts/statistical_runtime_env.sh before R/Python/Nextflow use. Wrapper defines project R user library, disables Python user-site, sets JAVA_HOME and project NXF_HOME/cache/tmp, and bounds BLAS/OpenMP threads to1.

## Validation and limits

Actual Rscript loaded all four packages, ran synthetic200genes/12samples through TMM, voom(span0.5), lmFit and robust eBayes: design rank2, residual df10, finite expression/t and finite positive weights PASS. This is explicitly a toy runtime test, not evidence of real biological performance. Actual Python numpy/scipy imports, java -version, nextflow -version PASS. Installer existing-prefix rejection also passed (expected exit1), protecting repeat execution from runtime mutation. Full pipeline execution remains unverified by this task.

Receipt files: receipts/runtime_statistics_R_smoke.json (sessionInfo/toy checks), receipts/runtime_statistics_install.json (exact package records, hashes, timings, memory, disk reserve). Raw logs/cache/environments are ignored by Git. No Docker access was retried.

Initial setup resource check failed before installation because shell quotation stripped tr argument; unnecessary tr removed and install passed. Receipt audit initially assumed m0-h5 was conda; corrected to its verified venv structure. Conda emitted implicit-defaults warning, but117 actual package channels were checked and wrapper now overrides channels. Java internal build suffix preserved honestly.

Official sources: [Bioconductor release compatibility](https://bioconductor.org/about/release-announcements/), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [Nextflow Java requirements](https://training.nextflow.io/latest/envsetup/02_local/).

YELLOW: understanding not assessed. User should distinguish toy software cooperation from count provenance, donor-design validity, scientific inference, external replication and container portability.
