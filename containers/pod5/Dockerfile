FROM mambaorg/micromamba:latest

USER root
RUN apt-get update
RUN apt-get install -y --no-install-recommends procps
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/bin:$PATH
