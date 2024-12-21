FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

WORKDIR /usr/local/bin

ARG MODE=555
COPY --chmod=${MODE} scripts/sc_get_metat_dict.py Snakefile .


ENV XDG_CACHE_HOME=/tmp/.cache