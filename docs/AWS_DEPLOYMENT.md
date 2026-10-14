# AWS deployment specification

**Status: documented; not executed.** Local native and Docker runs are recorded;
there is no AWS execution, cost or performance result for this project.

## Deployment scope

The proposed target is a Linux x86_64 EC2 host with encrypted storage for the
donor-level primary workflow. A single host is sufficient for the intended
orchestration design; this specification does not implement AWS Batch or repeat
raw-data acquisition. Instance size and storage must be selected from the actual
artifact inventory and measured workload, not inferred from cell counts alone.

## Required inputs

- A fixed source revision and the complete hashed artifact set referenced by the
  three scientific freezes, including frozen discovery outputs. A Git clone
  excludes required count matrices and intermediate products.
- The [locked scientific environment](../requirements/statistical-runtime.explicit.txt),
  Nextflow 25.10.2 and Java 21, or a rebuilt and verified Linux amd64 container.
  The recorded loopback registry address is local to the tested workstation and
  cannot be pulled from EC2 as a public image.
- The recorded `/home/analyst/te` path layout, or a separately tested portability
  adaptation. Frozen scripts and provenance records contain absolute paths.
- An authorized account, region, spending limit, scoped storage access and
  shutdown policy. No resources have been provisioned by this specification.

Source data and transfer bundles retain their original access and attribution
terms. Runtime credentials belong in the host identity mechanism, not in source,
containers, configuration committed to Git or evidence records.

## Proposed execution sequence

1. Provision the Linux host and required storage. Record the region, machine
   image, architecture, CPU, memory, storage and timestamps.
2. Restore the source and referenced artifacts. Recreate the locked environment
   or verify the rebuilt container, including its registry manifest and package
   inventory. Do not copy a machine-dependent environment prefix blindly.
3. From the restored project, verify all three scientific freezes:

   ```bash
   python scripts/freeze_analysis.py --verify receipts/statistical_pre_effect_freeze.json
   python scripts/freeze_analysis.py --verify receipts/discovery_rank_freeze.json
   python scripts/freeze_analysis.py --verify receipts/replication_implementation_freeze.json
   ```

4. Use new work/output directories for the designated discovery and external
   primary tasks. Adapt the tested commands in the [reproduction guide](LOCAL_EXECUTION.md)
   to the host and verify task mounts, user identity and input visibility before
   executing the scientific tasks. The existing workstation-specific checker is
   not an EC2 installer.
5. Compare outputs, preserve the execution record and shut down the allocated
   resources according to the retention policy.

## Acceptance criteria

Both primary workflows must complete. Selected gene IDs, directions, ranks,
filter decisions, tested-family sizes and agreement counts must match exactly.
Declare the numeric comparison tolerance before running; the tested local replay
used absolute tolerance 1e-8 and zero relative tolerance. Record any deviations.

Retain commands, input/container/configuration hashes, output comparisons,
execution traces, measured resources and costs, with estimates distinguished
from actual bills. Host paths and timestamps are expected to differ and do not
substitute for comparing scientific values.

A cloud replay would demonstrate portability of the same analysis, not independent
biological validation. Compute shutdown and remaining storage, snapshot and
registry resources require separate inventory checks; stopping an instance alone
does not establish that all charges have ended.
