# Local execution and reproduction

The user moved remaining execution from SRV01 to the local computer on 13 October 2026.
This instruction supersedes the old server-only location wording. The original scientific
analysis remains attributed to SRV01. Local files are at `G:\projects\te`
(`/mnt/g/projects/te`), on branch `codex/local-execution`, based on source commit
`4cc7580851cbfcb776385ac0ae858fe5f81e818b`. SRV01 originals are preserved.

## Current evidence

The full copied data/work/log/Nextflow tree contains 7,501 regular files and 495 links,
with approximately 3.48 GB of content. A complete rsync checksum dry-run found no content
difference; only the two destination parent-directory timestamps differ. The separate
initial migration check verified 1,607 declared hashes and all 500 bootstrap / 22 LODO
checkpoint records. These are integrity checks, not new scientific findings.

The local Docker build, real registry-manifest chain and actual workflow checks passed.
[The synthetic check](../receipts/local_container_nextflow.json) completed five fresh and
five cached tasks; [the real-primary replay](../receipts/local_primary_container_replay.json)
completed five fresh tasks and compared all 34 scientific TSVs plus eight JSON result
metadata files. Discrete choices matched exactly; the maximum continuous difference was
1.10e-13, below the declared absolute tolerance 1e-8. All 120 protected source/native-output
files retained their hashes. These checks establish reproduction, not new biological evidence.
The [final delivery verdict](../receipts/local_delivery_verification.json) also checks the
historical evidence, current documents, figures and PDF against their recorded hashes.

## Why these tools

R and Python perform the statistical analysis. Nextflow connects its steps, records task
execution and checks whether cached outputs may be reused. Docker supplies a fixed
software environment. The project uses them to meet the accepted reproducible-delivery
scope; they do not improve biological validity or replace donor-aware study design.
A standalone exploration can use ordinary scripts and a locked environment.

## Local runtime

`EndoTargetRuntime` is a project-owned WSL2 distribution with its virtual disk under
`work/local-runtime/distro` on G. It has a separate root filesystem but shares the WSL
kernel, networking and global 8 GiB memory / 2 GiB swap budget. The default swap file
can still use C. The existing Ubuntu filesystem and global `.wslconfig` were not changed.

The official Ubuntu Base 24.04.5 amd64 archive was checked against Canonical's HTTPS
SHA256 listing; its GPG signature was not checked. Setup prevented service autostart.
Docker runs as root **inside this new guest**, with a root-owned Unix socket only.
This is not rootless Docker and the guest is not an independent security VM.

The daemon uses no bridge, iptables/ip6tables editing, IP forwarding changes or masquerade.
Before/after network digests matched. Image/layer storage is Linux ext4 inside the G-hosted
virtual disk, not the Windows project mount. `/home/analyst/te` bind-mounts the same G
project to preserve old absolute evidence paths; this does not imply local runs used SRV01.

Scientific task containers use UID 1003 / GID 1004, no network, read-only root and project
inputs, a private tmpfs, and one writable task directory. Nextflow supplies the mounts.
Concurrency is one and RAM is capped at 2 GiB per task; the daemon reports **no per-container
swap-limit support**. Java uses a 512 MiB heap. Actual Docker inspect checks these settings.
Builtin seccomp and cgroup namespaces are recorded; AppArmor is not claimed.

## Reproduction entry points

After successful verification, the project daemon was [stopped with no active containers](../receipts/local_runtime_shutdown.json). All data, images, caches and the WSL distribution were preserved. Restart it only when another execution is needed.

Use a PowerShell window to start the existing guest daemon (leave it open while running):

~~~powershell
$env:PATHEXT='.COM;.EXE;.BAT;.CMD'
& 'C:\Program Files\WSL\wsl.exe' -d EndoTargetRuntime -u root -- /bin/bash /mnt/g/projects/te/scripts/local_runtime_start.sh
~~~

The PATHEXT setting affects only that PowerShell process. The launcher refuses to start
a duplicate daemon or to run in the wrong distribution. It restores the project bind mount
after a guest restart. It does not reinstall a completed guest.

From another PowerShell window, invoke the guest Python scripts using their
`/home/analyst/te/scripts/` paths. Relevant scripts are:

- `build_local_container.py`: build the new local recipe and inspect all 117 locked
  scientific packages plus 37 explicitly pinned OS packages supporting task metrics and Git provenance.
- `verify_local_container_manifest.py`: use a short-lived, non-root registry listening
  only on `127.0.0.1:5500`, verify raw manifest/config bytes and pull the true digest,
  then stop the registry. Nothing is pushed to an external service.
- `check_local_container_nextflow.py --image IMAGE@sha256:DIGEST`: two synthetic
  workflows, five fresh tasks and five cached tasks, with all 34 scientific TSVs compared.
- `check_local_primary_replay.py --image IMAGE@sha256:DIGEST`: one fresh discovery and
  one fresh HECA primary run; original frozen discovery outputs still drive HECA evaluation.

Use the exact reference in `receipts/local_container_manifest.json`. A Docker 29 image
ID can denote an OCI index; it must not be assumed to be an image-config digest.
The receipt records the verified index, amd64 manifest, configuration and pulled identity.

Completed receipts/output directories are deliberately preserved. Review and archive a
previous attempt under an explicit new name before an intentional repeat; do not erase
a receipt or silently overwrite results. The exact executed commands and paths are stored
in the receipts. Git alone excludes the required data/work artifact closure.

## Scope of the bounded primary replay

The original ROADMAP M5 / ETE-12 requires a clean environment to reproduce designated
outputs. Synthetic tests alone cannot establish this for the delivered scientific tables.
The local primary replay therefore repeats only three discovery tasks and two HECA tasks,
once each, under the frozen parameters with fresh output/work directories. Discrete
identities, ranks, counts and filter decisions must agree exactly; declared continuous
numbers use absolute tolerance 1e-8, with no relative tolerance. The real replay completed in
96.46 seconds wall time, including driver/check work; synthetic fresh/resume checks took
87.93 seconds. These are these particular local runs, not generalized performance benchmarks.

This is the necessary limited replay authorized by the original full-project request,
not a production-run authorization inferred from the earlier container-candidate guide.
It does not repeat raw import, 500 bootstrap fits, 22 donor-exclusion fits or clinical
sensitivities. Host paths/times and runtime serialization are handled separately from
scientific values. No new biological or independent-cohort evidence is claimed.

## Historical evidence and known fixes

The image recipe defaults to UID/GID 57439:57439; actual workflow tasks explicitly use
1003:1004 for the project files. Both are non-root. Build verification checks the default
against the recipe, and task verification separately checks each actual container. An
initial delivery-checker assertion conflated these values; its failed receipt is retained.

The three scientific freezes and all historical scientific outputs remain unchanged.
The two previous root READMEs are exact byte archives under `docs/archive/`; current
root READMEs describe local delivery. Use `verify_local_delivery.py` for the explicitly
versioned final state. The initial `verify_local_migration.py` checks the original root
README locations and is a historical migration check; after the documented homepage
update it is expected to reject those two changed root files.

Preserved failed attempts explain the fixes: the registry image's default debug listener
was caught and the service stopped before pushing project layers; an explicit configuration
now leaves only the loopback registry port. Nextflow's automatic mounts made the candidate's
manual mounts redundant; the new overlay removes the duplicates. The first image lacked
`ps`; the additive local recipe pins procps and its two added libraries. The original
117 scientific package lock is unchanged. The subsequent Git-enabled recipe pins all 37 selected
OS packages, including the metrics utilities and Git dependencies; it reads the actual Git
commit as UID 1003 with the project read-only.

A direct Java JAR launch initially omitted the official Nextflow 25.10.2 launcher's 17 JVM
module-opening arguments. Restoring that fixed-version launcher behavior allowed task
cache records to be written and reused; no particular underlying reflection exception
was captured. Synthetic checks use `work/synthetic/local-container` to meet the existing
synthetic-data path guard. Frozen workflow code and guards were not relaxed.

AWS remains a documented, unexecuted plan. No cloud resources, external image publication
or GitHub push is included in this local run.

YELLOW — understanding not yet assessed: statistical validity versus engineering
reproduction; task cache identity; donor replication; image index versus config identity.

## Observed resource snapshot

At 02:41:29 UTC, after the synthetic check and before the real replay, the [resource receipt](../receipts/local_container_resources.json) recorded 8.33 GB total guest RAM, 6.66 GB available, and 442.8 GB free on G. Guest ext4 used 7.70 GB. Docker reported images 5.528 GB and build cache 3.224 GB; shared layers mean these figures must not be added as independent physical usage. These are snapshot values, not peak memory, final disk usage or cloud cost. All failed-attempt artifacts and earlier images remain preserved.
