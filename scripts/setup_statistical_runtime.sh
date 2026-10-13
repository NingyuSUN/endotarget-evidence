#!/usr/bin/env bash
set -euo pipefail
cd /home/analyst/te
export CONDA_PKGS_DIRS="$PWD/cache/conda-pkgs"
export CONDA_ENVS_PATH="$PWD/envs"
export XDG_CACHE_HOME="$PWD/cache"
export NXF_HOME="$PWD/cache/nextflow"
export TMPDIR="$PWD/tmp"
mkdir -p "$CONDA_PKGS_DIRS" "$TMPDIR" logs/runtime
available=$(df -B1 --output=avail . | tail -n1 )
# Reserve 20 GiB plus a conservative 6 GiB unpack/cache budget.
(( available >= 26*1024*1024*1024 )) || { echo Insufficient resource reserve >&2; exit 1; }
[[ ! -e envs/statistical-runtime/conda-meta/history ]] || { echo Refusing to mutate an existing runtime >&2; exit 1; }
/usr/bin/time -v -o logs/runtime/conda-install.resources.txt /home/analyst/miniconda3/bin/conda create --override-channels -c conda-forge -c bioconda --yes --prefix "$PWD/envs/statistical-runtime" --file requirements/statistical-runtime-planned.explicit.txt --json > logs/runtime/conda-install.json 2> logs/runtime/conda-install.stderr
