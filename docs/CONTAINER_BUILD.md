# Statistical container candidate v1 — NOT_BUILT_NOT_RUN

Scope: frozen discovery and replication statistical tasks only. No ingestion tests,
scientific reruns, frozen-source edits, privilege changes, global installs or registry pushes.
Docker account access is blocked; daemon health remains unknown. Static checks cannot
establish that this image builds or runs.

## Inputs and immutable identity

The existing requirements/statistical-runtime.explicit.txt selects exact Linux-64 package
URLs/builds/MD5 fragments, installed at /opt/endo without fresh dependency solving.
The lock SHA256 is in receipts/container_recipe_v1.json; MD5 fragments are integrity
checks, not an independent supply-chain audit.

Pins: R 4.4.3, limma 3.62.1, edgeR 4.4.0, statmod 1.5.2, jsonlite 2.0.0,
Python 3.10.22, numpy 2.2.6, scipy 1.15.2, OpenJDK 21.0.10, Nextflow 25.10.2.
The Docker Hub registry API returned mambaorg/micromamba:2.3.2 index
sha256:ce2026639d612fbe837a35ca9d513292a8d776a24a5e6c82cccf9dddf45ad3d2,
whose header matched the SHA256 computed from 2364 bytes of index JSON.
Its Linux amd64 manifest is
sha256:955819619f303e2aa1dc1ba89beefe7f326a378b3d0782a4a0d28b1cf11b68b6.
No image/package layers were downloaded. This is a registry identity receipt, not image validation.

PATH/JAVA_HOME/R_LIBS_USER select /opt/endo. Do not source the frozen native
statistical_runtime_env.sh inside the image: it selects the host environment.
The frozen workflows already omit that wrapper in docker mode. The host Nextflow driver
still uses its pinned native environment and host cache/log directories.

## Build once authorized server access is available

Run on SRV01; these commands have NOT been executed:

~~~bash
cd /home/analyst/te
context=$(mktemp -d work/statistics_container_build_v1.XXXXXX)
cp containers/statistics/Dockerfile containers/statistics/.dockerignore "$context/"
cp requirements/statistical-runtime.explicit.txt "$context/"
docker build --platform linux/amd64 -t endotarget-statistics:v1 "$context"
~~~

The context contains recipe, ignore rules and explicit lock only. No source tree, data,
secrets, credentials, host environments or results enter the image. Record command,
lock/recipe hashes, logs, platform and image identity; failed builds remain failures.
Dependencies install only at build time, not inside statistical tasks.

## Local Docker smoke before registry publication

After build, inspect endotarget-statistics:v1 and record its local image Id.
That Id hashes image configuration; it is NOT a registry manifest digest and cannot
be relabeled as IMAGE@sha256. Frozen workflows reject docker mode without a real
repository manifest reference. Do not bypass their checks.

A direct Docker smoke may use the recorded local image Id. Create a new task directory
under work/synthetic/; bind /home/analyst/te to the identical container path ro and bind
only that task directory rw, pass --user "$(id -u):$(id -g)", --read-only, and a /tmp
tmpfs with the same uid/gid, mode 0700. Use project-root cwd for frozen scripts;
send outputs only to the writable task path. First check exact R/Python/Java/package
versions and read/write behavior. Frozen script tests then need a synthetic fixture/freeze
plan with review and integrity checks. The user has already
authorized container implementation/testing on SRV01: actual build and synthetic runs may
proceed once authorized server access is available, without renewed human approval.
This does not authorize production reruns.

Full Nextflow requires a real accessible IMAGE@sha256 registry manifest. An authenticated
ECR push or a reviewed local registry can provide it later. Privileged account changes
and external registry publication require separate authorization; neither is executed here.
Record the image manifest and overlay hash in a new
container-execution receipt. Old scientific freezes do not cover these runtime additions.

## Separate overlay, frozen files preserved

The overlay is conf/statistics_container.config; use both profiles. Illustrative discovery
command (replace ALL placeholders with reviewed synthetic paths/digest):

~~~bash
source scripts/statistical_runtime_env.sh
export ENDOTARGET_CONTAINER_UID="$(id -u)" ENDOTARGET_CONTAINER_GID="$(id -g)"
nextflow -c conf/statistics_container.config run main.nf -profile docker,statistics_container --container_image IMAGE@sha256:REAL_MANIFEST_DIGEST --pre_effect_freeze SYNTHETIC_FREEZE --synthetic_test true --counts SYNTHETIC_COUNTS --design SYNTHETIC_DESIGN --clinical SYNTHETIC_CLINICAL --edges SYNTHETIC_EDGES --outdir NEW_SYNTHETIC_OUTPUT -work-dir NEW_SYNTHETIC_WORK
~~~

For replication select replication.nf with -C replication.config,conf/statistics_container.config,
exported numeric UID/GID, the same profiles and an integrity-checked synthetic replication parameter file.
Do not edit main.nf/nextflow.config/replication.nf/replication.config or old manifests.

Containers use host numeric UID/GID and readonly root filesystem. Original project paths
remain available ro for hash guards. The exact active task directory is nested rw inside
that mount; upstream outputs remain readable. TMPDIR, HOME and cache use private writable
/tmp tmpfs. The host publishes output after tasks finish. No chmod 777 is used.

Before a future task runs, inspect generated .command.run for project ro, task rw,
correct UID/GID and tmpfs, and no broad writable bind overriding the intended mounts.
Nextflow 25.10.2 DockerBuilder explicitly handles writableInputMounts=false and appends
runOptions after input mounts; actual nested mount enforcement remains NOT VERIFIED.
The image arbitrary numeric user, /opt/endo prefix and package ABI compatibility also
remain NOT VERIFIED.

## Acceptance boundary

Current receipt proves only static package-pin references and preserved frozen file hashes.
It does not prove image build/run, mount enforcement, cache resume, clean-container execution,
daemon health, scientific model suitability or ingestion correctness. Actual goal completion
still requires a built image, isolated runtime/filesystem smoke, and integrity-checked artificial
Nextflow execution. Current Docker account access remains the blocker.

Primary references:
- https://github.com/mamba-org/micromamba-docker
- https://docs.docker.com/build/building/context/
- https://docs.docker.com/engine/storage/bind-mounts/
- https://github.com/nextflow-io/nextflow/blob/v25.10.2/modules/nextflow/src/main/groovy/nextflow/container/DockerBuilder.groovy

YELLOW — understanding not yet assessed: defend local image Id vs registry manifest digest,
software lock vs tested execution, and nested project ro/task rw mounts.
