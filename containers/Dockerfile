# One Dockerfile, one image per process spec.
#
#   docker build -f containers/Dockerfile --build-arg SPEC=kraken2 \
#                -t biothreatfinder-kraken2:0.4.0 .
#
# The spec files under envs/process/ are GENERATED from the `conda` directive of the
# processes that use them (scripts/sync_container_specs.py), so an image and the conda
# route it mirrors cannot drift: they are built from the same pinned line. A pin moves in
# a module, the generator rewrites the spec, the suite fails until the image is rebuilt.
#
# WHY NOT ONE IMAGE WITH EVERYTHING. Measured, not assumed: installing kraken2 beside
# mmseqs2/hmmer/mafft makes the solver resolve kraken2 2.0.8-beta (2019), which has no
# --report-minimizer-data - so the index builds, works, and carries none of the evidence
# stage 1 scores on (D49). A fat image reproduces that trap. These stay small and
# single-purpose, which also means they build in parallel and a rebuild touches one.
#
# WHY micromamba AND NOT A BIOCONTAINER. Eight of the thirteen specs are one tool and
# would map to a single biocontainer. Five are two or three tools together
# (minimap2+samtools, samtools+python, mmseqs2+seqkit, hmmer+seqkit) and would need a
# "mulled" image whose name is a hash of the exact package set - which exists only if
# someone has already built that combination. Building from our own pins is fewer moving
# parts and a stronger guarantee.
FROM mambaorg/micromamba:2.0.5-debian12-slim

# The spec to install, by file name under envs/process/ without .yml.
ARG SPEC
# Repeated here because a FROM resets ARGs declared before it.
ARG VERSION=dev

LABEL org.opencontainers.image.title="bioThreatFinder ${SPEC}" \
      org.opencontainers.image.description="Pinned tools for the bioThreatFinder \
processes whose conda directive is envs/process/${SPEC}.yml. Research software; not a \
validated diagnostic assay - see DISCLAIMER.md." \
      org.opencontainers.image.source="https://github.com/duceppemo/bioThreatFinder" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${VERSION}"

USER root
# procps for `ps`, which Nextflow calls to collect per-task resource metrics: without it
# every task still runs and the trace file silently reports nothing, which is the kind of
# quiet degradation this project spends most of its guards on.
RUN apt-get update -qq \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/process/${SPEC}.yml /tmp/spec.yml
RUN micromamba install -y -n base -f /tmp/spec.yml \
 && micromamba clean --all --yes \
 && rm /tmp/spec.yml

# The disclaimer travels with the image for the same reason it travels inside the database
# archives: whoever runs the container is often not whoever read the README.
COPY --chown=$MAMBA_USER:$MAMBA_USER DISCLAIMER.md /opt/bioThreatFinder/DISCLAIMER.md

# Nextflow runs each task through `bash -ue`, not through an entrypoint, so the tools have
# to be on PATH for a plain non-login shell rather than activated by one.
ENV PATH=/opt/conda/bin:$PATH
