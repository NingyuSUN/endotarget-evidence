# Resource gate — 2026-10-13

Status: OPEN / full-input path not demonstrated. This is a capacity calculation, not an observed out-of-memory failure.

SRV01 project filesystem had 69,935,509,504 free bytes (65.13 GiB) at the recorded probe. The host has approximately 31 GiB total RAM. Preserve at least 20 GiB free disk for the existing server.

- Discovery RDS: 41,729,517,596 bytes (38.86 GiB), gzip prefix verified via a 16-byte HTTP range read.
- Replication H5AD: 8,739,186,814 bytes (8.14 GiB), HDF5 prefix verified via a 16-byte HTTP range read.
- Together: 50,468,704,410 bytes. After both, only 18.13 GiB would remain BEFORE conversion, environments, caches and outputs.

Neither complete object has been downloaded. The archive size does not determine RDS in-memory size. A compressed 38.86 GiB Seurat object is not a defensible whole-object import plan for this host without additional measured evidence.

Next allowed work: inspect small public metadata/clinical supplements and public smaller-count representations, or investigate bounded HDF5 metadata access. Do not alter cohort, phase or donor thresholds to evade feasibility. Do not modify/delete existing server files to create room. If no suitable discovery representation exists, report the need for revised resources or a revised protocol before proceeding.

Evidence: receipts/input_head_probes.json, receipts/srv01_environment.json, receipts/resource_budget.json. Docker access is a separate runtime issue documented in docs/RUNTIME_STATUS.md.
