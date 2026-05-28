# Multi-target Dockerfile for Oncodrive3D.
#
# Targets (each builds on the previous; pick with `--target`):
#   light    (default)  - oncodrive3d only.
#                         Covers `run` and `plot` (including normal-tissue
#                         mutability runs). No external prerequisites to build.
#   chimerax            - light + ChimeraX. Adds `chimerax-plot`. Requires the
#                         ChimeraX installer staged in the build context as
#                         `chimerax.deb` (override via --build-arg CHIMERAX_DEB=).
#   full                - chimerax + bgdata caches (hg38, mm39) + PDB_Tool.
#                         Adds `build-datasets` and `build-annotations`.
#
# Examples:
#   docker build --target light    -t oncodrive3d:1.0.9         -t oncodrive3d:latest   .
#   docker build --target chimerax -t oncodrive3d:1.0.9-chimerax -t oncodrive3d:chimerax .
#   docker build --target full     -t oncodrive3d:1.0.9-full     -t oncodrive3d:full     .

# Stage 1: build the oncodrive3d wheel
FROM ghcr.io/astral-sh/uv:0.5-python3.10-bookworm AS build-stage

ENV UV_COMPILE_BYTECODE=1

WORKDIR /oncodrive3d

COPY uv.lock uv.lock
COPY pyproject.toml pyproject.toml
COPY scripts scripts
COPY README.md README.md

RUN mkdir -p /root/.cache/uv \
    && uv sync --frozen --no-dev \
    && uv build

# Stage 2: compile PDB_Tool from source (only consumed by the `full` target)
FROM python:3.10-bullseye AS pdb-tool-stage
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates git \
    && git clone --depth 1 https://github.com/realbigws/PDB_Tool.git /pdb-tool \
    && make -C /pdb-tool/source_code \
    && rm -rf /var/lib/apt/lists/*

# Stage 3 (target=light): just oncodrive3d. No ChimeraX, no PDB_Tool, no bgdata.
# Debian 11 base keeps layers shareable with the downstream targets that need
# libssl1.1 + libffi7 for the Ubuntu 20.04 ChimeraX .deb.
FROM python:3.10-bullseye AS light

COPY --from=build-stage /oncodrive3d/dist /oncodrive3d/dist
WORKDIR /oncodrive3d
RUN pip install --no-cache-dir dist/*.tar.gz

LABEL Author="stefano.pellegrini@irbbarcelona.org"

CMD ["oncodrive3D"]

# Stage 4 (target=chimerax): light + ChimeraX
FROM light AS chimerax

ENV DEBIAN_FRONTEND=noninteractive

ARG CHIMERAX_DEB=chimerax.deb
COPY ${CHIMERAX_DEB} /tmp/chimerax.deb

RUN apt-get update \
    && apt-get install -y --no-install-recommends /tmp/chimerax.deb \
    && rm /tmp/chimerax.deb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Stage 5 (target=full): chimerax + bgdata caches (hg38, mm39) + PDB_Tool
FROM chimerax AS full

ENV BGDATA_LOCAL="/bgdatacache" \
    BBGLAB_HOME="/home/user/.config/bbglab/"

# PDB_Tool binary (for `build-annotations`)
COPY --from=pdb-tool-stage /pdb-tool/PDB_Tool /usr/local/bin/PDB_Tool

RUN install -d -m 0755 "$BGDATA_LOCAL" "$BBGLAB_HOME"

# Write bgdata configuration
RUN echo "# Version of the bgdata config file\n\
version = 2\n\
\n\
# The default local folder to store the data packages\n\
local_repository = \"$BGDATA_LOCAL\"\n\
\n\
# The remote URL to download the data packages\n\
remote_repository = \"https://bbglab.irbbarcelona.org/bgdata\"\n\
\n\
# If you want to force bgdata to work only locally\n\
# offline = True\n\
\n\
# Cache repositories\n\
[cache_repositories]" > "$BBGLAB_HOME/bgdatav2.conf"

# Pre-fetch and prepare genome data (needed by `build-datasets`)
RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && pip install --no-cache-dir bgdata bgreference \
    && bgdata get datasets/genomereference/hg38 \
    && bgdata get datasets/genomereference/mm39 \
    && python3 -c "from bgreference import hg38; hg38(1, 1300000, 3000)" \
    && python3 -c "from bgreference import mm39; mm39(1, 1300000, 3000)" \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN chmod -R 0755 "$BGDATA_LOCAL"
