# Build stage image
FROM benjamingrodner/metat-builder AS builder
ARG MAMBA_DOCKERFILE_ACTIVATE=1
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
ARG LBIN=/usr/local/bin
WORKDIR ${LBIN}
ARG MODE=555
COPY --chmod=${MODE} scripts/* functions/* Snakefile _entrypoint.sh ./
# Variables
ENV XDG_CACHE_HOME=/tmp/.
ENV PYTHONPATH=${LBIN}
# When image runs, environment is activated
# COPY --chmod=${MODE} _entrypoint.sh /usr/local/bin/_entrypoint.sh
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
# Default command for docker run
CMD ["/bin/bash"]
