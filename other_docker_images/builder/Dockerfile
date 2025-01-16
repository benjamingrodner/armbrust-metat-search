# Build stage image
FROM mambaorg/micromamba:2.0.5
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