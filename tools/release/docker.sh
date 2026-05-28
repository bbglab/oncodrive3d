#!/usr/bin/env bash
#
# Build and push all three Oncodrive3D Docker image variants in one shot.
#
# Usage:
#   ./tools/release/docker.sh <version> <repo>
#
# Examples:
#   ./tools/release/docker.sh 1.0.9 spellegrini87/oncodrive3d
#   ./tools/release/docker.sh 1.0.9 bbglab/oncodrive3d
#
# Prerequisites:
#   - Logged in to Docker Hub for <repo> (`docker login`).
#   - `chimerax.deb` staged in the repo root (used by the chimerax and full
#     variants — see Dockerfile).

set -euo pipefail

VERSION="${1:?usage: docker.sh <version> <repo>}"
REPO="${2:?usage: docker.sh <version> <repo>}"

# Run from the repo root so `chimerax.deb` and the Dockerfile resolve.
cd "$(dirname "$0")/../.."

if [[ ! -f chimerax.deb ]]; then
    cat <<EOF >&2
Missing chimerax.deb in repo root ($(pwd)).

Download the Ubuntu 20.04 ChimeraX .deb from:
  https://www.cgl.ucsf.edu/chimerax/download.html

and save it as: $(pwd)/chimerax.deb
EOF
    exit 1
fi

docker build --target light    -f Dockerfile \
    -t "${REPO}:${VERSION}"          -t "${REPO}:latest" \
    -t "${REPO}:${VERSION}-light"    -t "${REPO}:light"    .

docker build --target chimerax -f Dockerfile \
    -t "${REPO}:${VERSION}-chimerax" -t "${REPO}:chimerax" .

docker build --target full     -f Dockerfile \
    -t "${REPO}:${VERSION}-full"     -t "${REPO}:full"     .

docker push "${REPO}" --all-tags

echo
echo "Released ${VERSION} to ${REPO}. Tags now on Docker Hub:"
docker images "${REPO}" --format "  {{.Tag}}"
