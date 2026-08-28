#!/usr/bin/env bash
#
# Regenerate python/.build/env-{312,313,314}-linux.yml from clean, minimal specs.
#
# The committed env-*.yml files are fully-pinned `conda env export` dumps. They
# are not meant to be hand-edited: over time they drift (channel order differs
# per file, env-314-linux.yml was exported from a live dev environment and
# carried ~140 unpinned pip entries including magtense itself). This script is
# the supported way to change them - edit the spec below, run this on a Linux
# x86 box, commit the result.
#
# Usage:
#   bash python/.build/regen-envs-linux.sh [312|313|314 ...]
#
# Run from the repository root. Requires conda on PATH (or CONDA_BIN set).

set -euo pipefail

CONDA_BIN="${CONDA_BIN:-conda}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/python/.build"

# Latest CUDA 12.x on the nvidia channel. Bump this one line to move CUDA.
CUDA_LABEL="cuda-12.9.2"

VERSIONS=("$@")
[ ${#VERSIONS[@]} -eq 0 ] && VERSIONS=(312 313 314)

for V in "${VERSIONS[@]}"; do
    case "$V" in
        312) PY=3.12; ONEAPI="2025.2.*" ;;
        313) PY=3.13; ONEAPI="2025.2.*" ;;
        314) PY=3.14; ONEAPI="2025.3.*" ;;
        *)   echo "unknown python version: $V" >&2; exit 1 ;;
    esac

    ENV_NAME="magtense-regen-${V}"
    SPEC="$(mktemp -t "spec-${V}-XXXXXX.yml")"
    trap 'rm -f "$SPEC"' EXIT

    # Direct dependencies only - conda solves the rest, which is what keeps the
    # export clean. numpy is a conda dependency here rather than a pip one:
    # scipy/matplotlib/h5py pull conda's numpy in regardless, and mixing that
    # with a pip-installed numpy in the same prefix is how you get an f2py build
    # that links one ABI and imports another.
    cat > "$SPEC" <<SPEC_EOF
name: ${ENV_NAME}
channels:
  - https://software.repos.intel.com/python/conda
  - nvidia/label/${CUDA_LABEL}
  - conda-forge
dependencies:
  - python=${PY}
  - pip
  - cmake

  # Intel oneAPI toolchain (ifx / icx) and MKL
  - ifx_linux-64=${ONEAPI}
  - dpcpp_linux-64=${ONEAPI}
  - mkl
  - mkl-devel
  - mkl-include
  - mkl-static

  # GNU toolchain - nvcc needs a supported host compiler
  - gcc_linux-64=13.*
  - gxx_linux-64=13.*

  # CUDA
  - cuda-nvcc
  - cuda-cudart-dev
  - libcublas-dev
  - libcusparse-dev
  - libnvjitlink-dev

  # Numerics and visualisation, mirroring requirements-py3.txt
  - numpy=2.3.2
  - scipy>=1.16.0
  - h5py>=3.14.0
  - matplotlib>=3.10.5
  - ipympl>=0.9.7
  - pycairo>=1.28.0
  - notebook>=7.4.5
  - tqdm>=4.67.1
  - importlib_resources>=6.5.2

  # f2py and wheel building
  - pip:
      - build>=1.3.0
      - charset-normalizer>=3.4.3
      - meson>=1.8.3
      - ninja>=1.13.0
      - packaging>=25.0
      - pyproject-hooks>=1.2.0
SPEC_EOF

    echo "==> regenerating env-${V}-linux.yml (python ${PY}, ${CUDA_LABEL})"
    "$CONDA_BIN" env remove -n "$ENV_NAME" -y >/dev/null 2>&1 || true
    "$CONDA_BIN" env create -n "$ENV_NAME" -f "$SPEC"

    TARGET="${OUT_DIR}/env-${V}-linux.yml"
    "$CONDA_BIN" env export -n "$ENV_NAME" \
        | grep -v '^prefix: ' \
        | sed "s/^name: ${ENV_NAME}$/name: magtense-env/" \
        > "$TARGET"

    # The export must describe a *build* environment. magtense appearing in its
    # own build env means the export came from a polluted prefix - the exact
    # failure mode that produced the old env-314-linux.yml.
    if grep -qiE '^\s*-\s*magtense([=<> ]|$)' "$TARGET"; then
        echo "ERROR: ${TARGET} contains magtense itself - refusing to write." >&2
        exit 1
    fi
    if ! grep -q "cuda-version=" "$TARGET"; then
        echo "ERROR: ${TARGET} has no cuda-version pin." >&2
        exit 1
    fi

    echo "    wrote ${TARGET}"
    grep -E '^\s+- (cuda-version|cuda-nvcc|libcublas|libcusparse|libnvjitlink)=' "$TARGET" || true

    "$CONDA_BIN" env remove -n "$ENV_NAME" -y >/dev/null 2>&1 || true
    rm -f "$SPEC"
    trap - EXIT
done

echo
echo "Done. Review the diffs, then run python/.build/test-envs-linux.sh"
