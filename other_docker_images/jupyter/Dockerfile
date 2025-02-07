# Pull from builder
FROM benjamingrodner/metat-builder
# Update with new packages
# COPY --chown=$MAMBA_USER:$MAMBA_USER env_jupyter.yml /tmp/env_jupyter.yml
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install --yes \
    -c conda-forge -c bioconda -c anaconda \
    'jupyterhub-singleuser' \
    'jupyterlab' \
    'nbclassic' \
    # Sometimes, when the new version of `jupyterlab` is released, latest `notebook` might not support it for some time
    # Old versions of `notebook` (<v7) didn't have a restriction on the `jupyterlab` version, and old `notebook` is getting installed
    # That's why we have to pin the minimum notebook version
    # More info: https://github.com/jupyter/docker-stacks/pull/2167
    'notebook>=7.2.2' \
    'umap-learn' \
    'scipy' \
    'statsmodels' \
    'scikit-learn' \
    'matplotlib' \
    'seaborn' \
    && \
    jupyter server --generate-config && \
    micromamba clean --all -f -y && \
    jupyter lab clean

WORKDIR /app
ENV XDG_CACHE_HOME=/app/.cache
ENV XDG_CONFIG_HOME=/app/.config
ENV XDG_DATA_HOME=/app/.local/share
ENV XDG_BIN_HOME=/app/.local/bin
ENV IPYTHONDIR=/app/.local/var