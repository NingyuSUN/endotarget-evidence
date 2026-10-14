# Statistical container

The verified local recipe is
[`containers/statistics-local/Dockerfile`](../containers/statistics-local/Dockerfile).
It supports the three discovery and two external-evaluation Nextflow tasks.

## Environment and identity

- The [explicit scientific lock](../requirements/statistical-runtime.explicit.txt)
  fixes 117 packages, including R 4.4.3, edgeR 4.4.0, limma 3.62.1,
  Python 3.10.22, NumPy 2.2.6 and SciPy 1.15.2.
- The [OS package lock](../containers/statistics-local/os-packages.lock.tsv)
  fixes 37 packages supporting process metrics and Git provenance.
- The [build record](../receipts/local_container_build.json) contains the actual
  package inspection. The [manifest record](../receipts/local_container_manifest.json)
  binds the registry index, Linux amd64 manifest, image configuration and pulled
  identity. An image ID must not be substituted for a registry manifest digest.

The recorded image reference belongs to a local registry. No public container
image is distributed with this repository.

## Build and replay requirements

The [reproduction guide](LOCAL_EXECUTION.md) describes the tested WSL2 runtime,
fixed project paths, build driver, manifest checker and necessary saved artifacts.
The build driver constructs a bounded context from the recipe and package locks;
scientific inputs, credentials and project outputs are not image layers.

The actual task overlay is [`conf/local_container.config`](../conf/local_container.config).
Nextflow controls input and task-directory mounts. Verified task containers use
UID 1003 / GID 1004, no network, a read-only root and inputs, a private temporary
filesystem and one writable task directory. The image default UID/GID is
57439:57439; the actual task user is inspected separately.

The earlier [`containers/statistics/Dockerfile`](../containers/statistics/Dockerfile)
and `conf/statistics_container.config` are preserved source artifacts. They are
not interchangeable with the tested local recipe and overlay.

## Recorded verification

Synthetic checks completed five fresh tasks and five cached tasks. The real
primary replay completed five fresh tasks and reproduced 34 scientific tables
and eight result metadata files. Discrete selections were exact; continuous
values met absolute tolerance 1e-8, with maximum difference 1.10e-13.

See [synthetic checks](../receipts/local_container_nextflow.json) and
[primary replay](../receipts/local_primary_container_replay.json). These checks
establish same-data reproduction in the recorded environment, not new biological
evidence, fresh-machine portability or AWS execution.
