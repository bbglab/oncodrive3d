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

Save the Ubuntu 20.04 ChimeraX .deb installer as:
  $(pwd)/chimerax.deb

How to get it:

  External users:
    Download from https://www.cgl.ucsf.edu/chimerax/download.html

  BBGLab internal users:
    Download from the lab cluster path:
      /data/bbg/datasets/oncodrive3d/chimerax/chimerax-rc_2004.deb
EOF
    exit 1
fi

LIGHT_TAGS=(
    "${REPO}:${VERSION}"
    "${REPO}:latest"
    "${REPO}:${VERSION}-light"
    "${REPO}:light"
)
CHIMERAX_TAGS=(
    "${REPO}:${VERSION}-chimerax"
    "${REPO}:chimerax"
)
FULL_TAGS=(
    "${REPO}:${VERSION}-full"
    "${REPO}:full"
)

# Expand a tag array into `-t tag1 -t tag2 ...` for docker build.
build_tag_args() {
    local args=() tag
    for tag in "$@"; do args+=(-t "$tag"); done
    printf '%s\n' "${args[@]}"
}

mapfile -t LIGHT_BUILD_ARGS    < <(build_tag_args "${LIGHT_TAGS[@]}")
mapfile -t CHIMERAX_BUILD_ARGS < <(build_tag_args "${CHIMERAX_TAGS[@]}")
mapfile -t FULL_BUILD_ARGS     < <(build_tag_args "${FULL_TAGS[@]}")

docker build --target light    -f Dockerfile "${LIGHT_BUILD_ARGS[@]}"    .
docker build --target chimerax -f Dockerfile "${CHIMERAX_BUILD_ARGS[@]}" .
docker build --target full     -f Dockerfile "${FULL_BUILD_ARGS[@]}"     .

# Push only the tags created above; avoid `--all-tags` so stale local tags
# (e.g. from prior releases) are not silently re-published.
for tag in "${LIGHT_TAGS[@]}" "${CHIMERAX_TAGS[@]}" "${FULL_TAGS[@]}"; do
    docker push "${tag}"
done

echo
echo "Released ${VERSION} to ${REPO}. Tags pushed:"
printf '  %s\n' "${LIGHT_TAGS[@]}" "${CHIMERAX_TAGS[@]}" "${FULL_TAGS[@]}"
