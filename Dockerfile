# Build stage image
FROM mambaorg/micromamba:2.0.5 AS builder
# Install environment
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
# Active environment in build
ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Install ete4 via git clone since branch is still active (as of 2024_01_10)
USER root
WORKDIR /venv
RUN mkdir -p ete4 && \
    git clone -b ete4 https://github.com/etetoolkit/ete ete4 && \
    cd ete4 && pip install -e .
# Install conda-pack:
RUN micromamba install -c conda-forge conda-pack
# Standalone environment
RUN conda-pack -p /opt/conda -o /tmp/env.tar && \
    tar xf /tmp/env.tar && \
    rm /tmp/env.tar
# Same path as final image 
RUN bin/conda-unpack

# Runtime image
FROM debian:12-slim
# Venv from previous stage
COPY --from=builder /venv /venv
# Copy scripts
WORKDIR /usr/local/bin
ARG MODE=555
COPY --chmod=${MODE} scripts/* Snakefile .
# Variables
ENV XDG_CACHE_HOME=/tmp/.cache
# When image runs, environment is activated
COPY --chmod=${MODE} _entrypoint.sh /usr/local/bin/_entrypoint.sh
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
# Default command for docker run
CMD ["/bin/bash"]
