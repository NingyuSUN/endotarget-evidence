# Candidate v1: NOT_BUILT_NOT_RUN; linux/amd64 manifest verified 2026-10-13.
FROM --platform=linux/amd64 mambaorg/micromamba@sha256:955819619f303e2aa1dc1ba89beefe7f326a378b3d0782a4a0d28b1cf11b68b6
USER root
COPY statistical-runtime.explicit.txt /tmp/statistical-runtime.explicit.txt
RUN micromamba create --yes --prefix /opt/endo --file /tmp/statistical-runtime.explicit.txt && micromamba clean --all --yes && rm /tmp/statistical-runtime.explicit.txt
ENV PATH=/opt/endo/bin:$PATH JAVA_HOME=/opt/endo/lib/jvm R_LIBS_USER=/opt/endo/lib/R/library PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 TMPDIR=/tmp HOME=/tmp XDG_CACHE_HOME=/tmp/.cache
# Skip base entrypoint activation of /opt/conda.
ENTRYPOINT []
USER 57439:57439
WORKDIR /tmp
CMD ["/bin/bash"]
