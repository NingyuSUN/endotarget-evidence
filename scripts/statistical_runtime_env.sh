# Source from project scripts; all caches and R/Python library resolution isolated.
export ENDOTARGET_ROOT=/home/analyst/te
export PATH="$ENDOTARGET_ROOT/envs/statistical-runtime/bin:$PATH"
export JAVA_HOME="$ENDOTARGET_ROOT/envs/statistical-runtime/lib/jvm"
export NXF_HOME="$ENDOTARGET_ROOT/cache/nextflow"
export XDG_CACHE_HOME="$ENDOTARGET_ROOT/cache"
export CONDA_PKGS_DIRS="$ENDOTARGET_ROOT/cache/conda-pkgs"
export CONDA_ENVS_PATH="$ENDOTARGET_ROOT/envs"
export TMPDIR="$ENDOTARGET_ROOT/tmp"
export R_LIBS_USER="$ENDOTARGET_ROOT/envs/statistical-runtime/lib/R/library"
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
