# Reproducing the primary analyses

The donor-level analyses were executed in a native SRV01 environment. Their five
primary Nextflow tasks were subsequently reproduced in a fixed Docker environment
on Windows/WSL2. This guide describes that tested environment, its required artifacts,
and the scope of the completed checks.

## Completed checks

| Check | Recorded outcome | Evidence |
|---|---|---|
| Fixed software environment | 117 locked scientific packages and 37 pinned OS packages checked in the built image | [Build receipt](../receipts/local_container_build.json) |
| Image identity | Actual registry index, amd64 manifest, configuration and pulled identity verified | [Manifest receipt](../receipts/local_container_manifest.json) |
| Synthetic workflow and resume | Five fresh tasks and five cached tasks completed; all 34 scientific TSVs compared | [Synthetic receipt](../receipts/local_container_nextflow.json) |
| Real primary reproduction | Five fresh tasks completed; 34 scientific TSVs and eight JSON result metadata files matched | [Primary replay receipt](../receipts/local_primary_container_replay.json) |
| Continuous software checks | On every push: synthetic unit tests (Python 3.10 on the pinned probe lock, with the checkout bind-mounted at the documented project path; R on the Bioconductor 3.20 release series) and a Markdown link check; software checks only, no scientific rerun | [Workflow](../.github/workflows/ci.yml), [test runner](../scripts/run_unit_tests.py) |

For the real replay, identifiers, counts, ranks, directions and filter decisions
matched exactly. Continuous scientific values used absolute tolerance 1e-8 and
zero relative tolerance; the largest difference was 1.10e-13. All 120 protected
source and native-output files retained their hashes. Runtime timestamps, paths,
host descriptions and RDS serialization were handled separately from scientific
values. These checks reproduce the same data and analysis; they do not add an
independent cohort or establish biological or therapeutic validity.

## Required artifacts and software

A Git clone contains source, configuration, reports and receipts. It excludes the
large inputs, environments and `work/` products needed to restore the full checkpoint.
The tested replay requires:

- The hashed discovery and HECA count matrices, STRING gene-edge reference, frozen
  designs, clinical ledgers and all artifacts required by the three scientific freezes.
- The original discovery and HECA outputs used as comparison references. HECA
  evaluation continues to use the original frozen discovery lists.
- The preserved synthetic fixtures referenced by
  [the container checker](../scripts/check_local_container_nextflow.py), plus the
  local runtime assets and fixed Nextflow 25.10.2 JAR under `work/local-runtime/`.
- The existing `EndoTargetRuntime` WSL2 distribution, its project bind mount and the
  verified Docker image. The recorded registry reference is local to that environment;
  it is not a publicly hosted image that a fresh clone can pull.

The native statistical environment used R 4.4.3, edgeR 4.4.0, limma 3.62.1,
statmod 1.5.2, Python 3.10.22, NumPy 2.2.6 and SciPy 1.15.2. The
[explicit Linux package lock](../requirements/statistical-runtime.explicit.txt)
records the complete scientific environment. Nextflow 25.10.2 and Java 21 run the
workflow. The scripts retain `/home/analyst/te` as their required project path;
the tested guest bind-mounts `G:\projects\te` (`/mnt/g/projects/te`) at that path.
Relocation is not a tested drop-in replacement.

R and Python perform the statistical analysis. Nextflow connects the steps and
records execution/cache identity. Docker fixes the software environment. None of
these tools replaces donor-aware design or improves the biological evidence.

## Starting the tested runtime

The completed environment was [stopped with no active containers](../receipts/local_runtime_shutdown.json).
From PowerShell, start the existing guest daemon and leave that window open during use:

```powershell
$env:PATHEXT='.COM;.EXE;.BAT;.CMD'
& 'C:\Program Files\WSL\wsl.exe' -d EndoTargetRuntime -u root -- /bin/bash /mnt/g/projects/te/scripts/local_runtime_start.sh
```

The `PATHEXT` setting affects only this PowerShell process. The launcher verifies
the distribution, refuses a duplicate daemon, and restores the project bind mount
after a guest restart. It does not install the runtime.

The following scripts execute as root inside `EndoTargetRuntime`, using their
`/home/analyst/te/scripts/` paths:

| Script | Operation |
|---|---|
| `build_local_container.py` | Build the local recipe and inspect the locked scientific and OS packages |
| `verify_local_container_manifest.py` | Check raw manifest/config bytes through a temporary non-root registry on `127.0.0.1:5500`, pull the real digest, then stop the registry |
| `check_local_container_nextflow.py --image IMAGE@sha256:DIGEST` | Execute the preserved synthetic discovery/replication fixtures and cached resumes |
| `check_local_primary_replay.py --image IMAGE@sha256:DIGEST` | Execute three discovery and two HECA tasks with fresh work/output directories |

Use the exact image reference recorded in
[`receipts/local_container_manifest.json`](../receipts/local_container_manifest.json).
`IMAGE@sha256:DIGEST` above is a placeholder, not a runnable reference. A Docker 29
image ID may denote an OCI index and must not be treated as an image-configuration
digest. Exact executed commands and paths are recorded in the receipts.

Completed receipts and output directories are deliberately retained. These are
checkpoint-specific reproduction tools, not an unattended fresh-clone installer.
Before a deliberate repeat, preserve the previous attempt under an explicit new
name; do not erase receipts or silently overwrite results.

## Execution controls and resource limits

`EndoTargetRuntime` has its virtual disk under `work/local-runtime/distro` on G.
It has a separate root filesystem but shares the WSL kernel, networking and global
8 GiB memory / 2 GiB swap budget. The default swap file may use C. The existing
Ubuntu filesystem and global `.wslconfig` were not changed. The Ubuntu Base
24.04.5 amd64 archive was checked against Canonical's HTTPS SHA256 listing;
its GPG signature was not checked.

Docker runs as root inside this guest with a root-owned Unix socket. It is not
rootless Docker, and the guest is not an independent security VM. The daemon uses
no bridge, iptables/ip6tables edits, IP-forwarding changes or masquerade; recorded
before/after network digests matched. Image/layer storage uses the guest's ext4
filesystem rather than the Windows project mount.

Actual workflow containers use UID 1003 / GID 1004, no network, a read-only root
and project inputs, a private tmpfs, and one writable task directory. Nextflow
supplies the mounts. The image's separate default user is UID/GID 57439:57439.
Concurrency is one, RAM is capped at 2 GiB per task, and Java uses a 512 MiB heap.
Docker reported no per-container swap-limit support. Builtin seccomp and cgroup
namespaces were recorded; AppArmor is not claimed. The receipts contain actual
container inspection rather than configuration-only assertions.

The real replay took 96.46 seconds including its driver/check work; synthetic
fresh/resume checks took 87.93 seconds. These are individual observed runs, not
performance benchmarks. The [resource snapshot](../receipts/local_container_resources.json)
between those checks recorded 8.33 GB total guest RAM, 6.66 GB available,
442.8 GB free on G and 7.70 GB guest ext4 usage. Docker image and build-cache
figures share layers and cannot be summed as independent physical storage.
These observations are not peak-memory, final-storage or cloud-cost measurements.

## Reproduction scope

The real replay repeats only the three discovery and two external primary tasks.
It does not repeat raw import, 500 bootstrap fits, 22 donor-exclusion fits or
clinical sensitivities. The original bootstrap and donor-exclusion checkpoints
were checked against recorded hashes. The three scientific freezes and original
scientific outputs remain unchanged.

[AWS deployment instructions](AWS_DEPLOYMENT.md) describe an unexecuted portability
plan. No cloud execution or external image publication was performed. The
[publication history](PUBLICATION_HISTORY.md) explains the scope of historical
documentation hashes and verification receipts after the public documentation update.
